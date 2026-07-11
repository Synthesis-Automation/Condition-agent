"""Leaving-group site detection."""

from __future__ import annotations

from typing import Any, Dict, List

from ..context import classify_context
from ..labels import render_edge
from ..patterns import matched_patterns_for_atom, matched_role_atoms
from .common import bond_index, unique_indices


def _sulfonate_token(mol: Any, sulfur: Any, connector_o: Any) -> str:
    carbon_neighbors = [n for n in sulfur.GetNeighbors() if n.GetSymbol() == "C"]
    if not carbon_neighbors:
        return "OSO2R"
    carbon = carbon_neighbors[0]
    fluorines = sum(1 for n in carbon.GetNeighbors() if n.GetSymbol() == "F")
    if fluorines == 3:
        return "OTf"
    if carbon.GetIsAromatic():
        return "OTs"
    if carbon.GetTotalNumHs() == 3:
        return "OMs"
    return "OSO2R"


def detect(mol: Any) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    candidate_handles = matched_role_atoms(mol, "leaving_group", "handle")
    candidate_connectors = matched_role_atoms(mol, "leaving_group", "connector")
    for handle in mol.GetAtoms():
        if handle.GetSymbol() in {"F", "Cl", "Br", "I"}:
            if handle.GetIdx() not in candidate_handles:
                continue
            # A leaving-group halogen is terminal. This excludes metal-halogen
            # bridges and malformed organometallic strings such as Mg-Br-Ar.
            if handle.GetDegree() != 1:
                continue
            for anchor in handle.GetNeighbors():
                if anchor.GetSymbol() != "C":
                    continue
                # Fluorines in retained polyfluoroalkyl groups (CF2/CF3) are
                # context, not individual alkyl-fluoride reactive handles.
                if handle.GetSymbol() == "F":
                    attached_halogen_count = sum(
                        neighbor.GetSymbol() in {"F", "Cl", "Br", "I"}
                        for neighbor in anchor.GetNeighbors()
                    )
                    if attached_halogen_count >= 2:
                        continue
                context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()})["token"]
                sites.append({
                    "topology": "edge", "atom_indices": [anchor.GetIdx(), handle.GetIdx()],
                    "bond_indices": [bond_index(mol, anchor.GetIdx(), handle.GetIdx())],
                    "signature": f"LG|{context}|{handle.GetSymbol()}",
                    "label": render_edge(context, handle.GetSymbol()),
                    "details": {"anchor_atom_index": anchor.GetIdx(), "handle_atom_indices": [handle.GetIdx()], "anchor_context": context, "handle_token": handle.GetSymbol()},
                    "matched_patterns": [item["id"] for item in matched_patterns_for_atom(mol, "leaving_group", "handle", handle.GetIdx())],
                })
        if handle.GetSymbol() != "O" or handle.GetDegree() != 2:
            continue
        if handle.GetIdx() not in candidate_connectors:
            continue
        sulfur_neighbors = [n for n in handle.GetNeighbors() if n.GetSymbol() == "S"]
        carbon_neighbors = [n for n in handle.GetNeighbors() if n.GetSymbol() == "C"]
        if not sulfur_neighbors or not carbon_neighbors:
            continue
        sulfur, anchor = sulfur_neighbors[0], carbon_neighbors[0]
        oxo = [n for n in sulfur.GetNeighbors() if n.GetSymbol() == "O" and str(mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondType()) == "DOUBLE"]
        if len(oxo) < 2:
            continue
        token = _sulfonate_token(mol, sulfur, handle)
        context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()})["token"]
        handle_atoms = unique_indices([handle.GetIdx(), sulfur.GetIdx(), *(n.GetIdx() for n in sulfur.GetNeighbors() if n.GetIdx() != handle.GetIdx())])
        sites.append({
            "topology": "edge", "atom_indices": unique_indices([anchor.GetIdx(), *handle_atoms]),
            "bond_indices": [bond_index(mol, anchor.GetIdx(), handle.GetIdx())],
            "signature": f"LG|{context}|{token}", "label": render_edge(context, token),
            "details": {"anchor_atom_index": anchor.GetIdx(), "handle_atom_indices": handle_atoms, "anchor_context": context, "handle_token": token},
            "matched_patterns": [item["id"] for item in matched_patterns_for_atom(mol, "leaving_group", "connector", handle.GetIdx())],
        })
    return sites
