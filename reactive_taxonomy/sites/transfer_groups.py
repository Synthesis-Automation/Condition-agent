"""Transfer-group site detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_context
from ..models import ReactiveSiteCandidate
from ..patterns import MatchIndex
from .common import bond_index, unique_indices


_METAL_TOKENS = {
    "Zn": "ZnX",
    "Mg": "MgX",
    "Li": "Li",
    "Cu": "Cu",
    "Al": "Al",
    "Sn": "SnR3",
    "Si": "SiR3",
}
_ANCHOR_CONTEXT_RANK = {"HeteroAr": 50, "Ar": 45, "Alkenyl": 40, "Alkynyl": 35, "Alkyl": 10}


def _boron_token(boron: Any) -> str:
    neighbors = list(boron.GetNeighbors())
    if sum(neighbor.GetSymbol() == "F" for neighbor in neighbors) >= 3:
        return "BF3K"
    oxygens = [neighbor for neighbor in neighbors if neighbor.GetSymbol() == "O"]
    if len(oxygens) >= 2:
        if sum(int(oxygen.GetTotalNumHs(includeNeighbors=True)) > 0 or oxygen.GetFormalCharge() < 0 for oxygen in oxygens) >= 2:
            return "B(OH)2"
        carbon_neighbors = [neighbor for oxygen in oxygens for neighbor in oxygen.GetNeighbors() if neighbor.GetSymbol() == "C"]
        if len({neighbor.GetIdx() for neighbor in carbon_neighbors}) == 2 and all(neighbor.GetDegree() >= 3 for neighbor in carbon_neighbors):
            return "BPin"
        return "B(OR)2"
    return "B"


def detect(mol: Any, match_index: MatchIndex) -> List[ReactiveSiteCandidate]:
    sites: List[ReactiveSiteCandidate] = []
    candidate_centers = match_index.role_atoms("transfer_group", "center")
    for handle in mol.GetAtoms():
        symbol = handle.GetSymbol()
        if (symbol != "B" and symbol not in _METAL_TOKENS) or handle.GetIdx() not in candidate_centers:
            continue
        token = _boron_token(handle) if symbol == "B" else _METAL_TOKENS[symbol]
        carbon_anchors = [neighbor for neighbor in handle.GetNeighbors() if neighbor.GetSymbol() == "C"]
        if symbol in {"Si", "Sn"}:
            if any(neighbor.GetSymbol() != "C" for neighbor in handle.GetNeighbors()):
                continue
            ranked = []
            for candidate in carbon_anchors:
                context = classify_context(mol, candidate.GetIdx(), {handle.GetIdx()}, match_index=match_index)
                bonus = 5 if candidate.GetTotalNumHs(includeNeighbors=True) < 3 else 0
                ranked.append((_ANCHOR_CONTEXT_RANK.get(context.token, 0) + bonus, candidate.GetIdx(), candidate))
            ranked.sort(key=lambda item: (-item[0], item[1]))
            if not ranked or ranked[0][0] <= 10 or (len(ranked) > 1 and ranked[0][0] == ranked[1][0]):
                continue
            carbon_anchors = [ranked[0][2]]
        for anchor in carbon_anchors:
            context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()}, match_index=match_index)
            handle_atoms = tuple(unique_indices([handle.GetIdx(), *(n.GetIdx() for n in handle.GetNeighbors() if n.GetIdx() != anchor.GetIdx())]))
            sites.append(ReactiveSiteCandidate(
                site_type="transfer_group", topology="edge",
                atom_roles={"anchor": (anchor.GetIdx(),), "center": (handle.GetIdx(),), "handle": handle_atoms},
                atom_indices=tuple(unique_indices([anchor.GetIdx(), *handle_atoms])),
                bond_indices=(bond_index(mol, anchor.GetIdx(), handle.GetIdx()),),
                canonical_signature=f"TM|{context.token}|{token}", render_kind="edge",
                render_data={"context": context.token, "handle": token},
                matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("transfer_group", "center", handle.GetIdx())),
                details={"anchor_context": context.token, "handle_token": token},
                context_records=(context,), availability="transferable",
            ))
    return sites
