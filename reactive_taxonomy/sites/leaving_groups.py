"""Leaving-group site detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_context
from ..models import SiteCandidate
from ..patterns import MatchIndex
from .common import bond_index, unique_indices


def _sulfonate_token(sulfur: Any) -> str:
    carbon_neighbors = [neighbor for neighbor in sulfur.GetNeighbors() if neighbor.GetSymbol() == "C"]
    if not carbon_neighbors:
        return "OSO2R"
    carbon = carbon_neighbors[0]
    if sum(neighbor.GetSymbol() == "F" for neighbor in carbon.GetNeighbors()) == 3:
        return "OTf"
    if carbon.GetIsAromatic():
        return "OTs"
    if carbon.GetTotalNumHs() == 3:
        return "OMs"
    return "OSO2R"


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    sites: List[SiteCandidate] = []
    candidate_handles = match_index.role_atoms("leaving_group", "handle")
    candidate_connectors = match_index.role_atoms("leaving_group", "connector")
    for handle in mol.GetAtoms():
        if handle.GetSymbol() in {"F", "Cl", "Br", "I"}:
            if handle.GetIdx() not in candidate_handles or handle.GetDegree() != 1:
                continue
            anchor = next(iter(handle.GetNeighbors()))
            if anchor.GetSymbol() != "C":
                continue
            if handle.GetSymbol() == "F" and sum(
                neighbor.GetSymbol() in {"F", "Cl", "Br", "I"} for neighbor in anchor.GetNeighbors()
            ) >= 2:
                continue
            context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()}, match_index=match_index)
            token = handle.GetSymbol()
            sites.append(SiteCandidate(
                site_type="leaving_group", topology="edge",
                atom_roles={"anchor": (anchor.GetIdx(),), "handle": (handle.GetIdx(),)},
                atom_indices=(anchor.GetIdx(), handle.GetIdx()),
                bond_indices=(bond_index(mol, anchor.GetIdx(), handle.GetIdx()),),
                canonical_signature=f"LG|{context.token}|{token}",
                render_kind="edge", render_data={"context": context.token, "handle": token},
                matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("leaving_group", "handle", handle.GetIdx())),
                details={"anchor_context": context.token, "handle_token": token},
                context_records=(context,), availability="available",
            ))
        if handle.GetSymbol() != "O" or handle.GetDegree() != 2 or handle.GetIdx() not in candidate_connectors:
            continue
        sulfur_neighbors = [n for n in handle.GetNeighbors() if n.GetSymbol() == "S"]
        carbon_neighbors = [n for n in handle.GetNeighbors() if n.GetSymbol() == "C"]
        if not sulfur_neighbors or not carbon_neighbors:
            continue
        sulfur, anchor = sulfur_neighbors[0], carbon_neighbors[0]
        oxo = [n for n in sulfur.GetNeighbors() if n.GetSymbol() == "O" and str(mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondType()) == "DOUBLE"]
        if len(oxo) < 2:
            continue
        token = _sulfonate_token(sulfur)
        context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()}, match_index=match_index)
        handle_atoms = tuple(unique_indices([handle.GetIdx(), sulfur.GetIdx(), *(n.GetIdx() for n in sulfur.GetNeighbors() if n.GetIdx() != handle.GetIdx())]))
        sites.append(SiteCandidate(
            site_type="leaving_group", topology="edge",
            atom_roles={"anchor": (anchor.GetIdx(),), "connector": (handle.GetIdx(),), "center": (sulfur.GetIdx(),), "handle": handle_atoms},
            atom_indices=tuple(unique_indices([anchor.GetIdx(), *handle_atoms])),
            bond_indices=(bond_index(mol, anchor.GetIdx(), handle.GetIdx()),),
            canonical_signature=f"LG|{context.token}|{token}",
            render_kind="edge", render_data={"context": context.token, "handle": token},
            matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("leaving_group", "connector", handle.GetIdx())),
            details={"anchor_context": context.token, "handle_token": token},
            context_records=(context,), availability="available",
        ))
    return sites
