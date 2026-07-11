"""Transfer-group site detection."""

from __future__ import annotations

from typing import Any, Dict, List

from ..context import classify_context
from ..labels import render_edge
from .common import bond_index, unique_indices


_METAL_TOKENS = {"Zn": "ZnX", "Mg": "MgX", "Sn": "SnR3", "Si": "SiR3"}

_ANCHOR_CONTEXT_RANK = {
    "HeteroAr": 50,
    "Ar": 45,
    "Alkenyl": 40,
    "Alkynyl": 35,
    "Alkyl": 10,
}


def _boron_token(boron: Any) -> str:
    neighbors = list(boron.GetNeighbors())
    f_count = sum(n.GetSymbol() == "F" for n in neighbors)
    oxygens = [n for n in neighbors if n.GetSymbol() == "O"]
    if f_count >= 3: return "BF3K"
    if len(oxygens) >= 2:
        oh_count = sum(int(o.GetTotalNumHs(includeNeighbors=True)) > 0 or o.GetFormalCharge() < 0 for o in oxygens)
        if oh_count >= 2: return "B(OH)2"
        carbon_neighbors = [n for o in oxygens for n in o.GetNeighbors() if n.GetSymbol() == "C"]
        if len({n.GetIdx() for n in carbon_neighbors}) == 2 and all(n.GetDegree() >= 3 for n in carbon_neighbors): return "BPin"
        return "B(OR)2"
    return "B"


def detect(mol: Any) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    for handle in mol.GetAtoms():
        symbol = handle.GetSymbol()
        if symbol != "B" and symbol not in _METAL_TOKENS:
            continue
        token = _boron_token(handle) if symbol == "B" else _METAL_TOKENS[symbol]
        carbon_anchors = [neighbor for neighbor in handle.GetNeighbors() if neighbor.GetSymbol() == "C"]
        if symbol in {"Si", "Sn"}:
            # V1 covers carbon-transfer silanes/stannanes, not silyl ethers or
            # other heteroatom-bound protecting groups.
            if any(neighbor.GetSymbol() not in {"C"} for neighbor in handle.GetNeighbors()):
                continue
            ranked = []
            for candidate in carbon_anchors:
                context = classify_context(mol, candidate.GetIdx(), {handle.GetIdx()})["token"]
                non_methyl_bonus = 5 if candidate.GetTotalNumHs(includeNeighbors=True) < 3 else 0
                ranked.append((_ANCHOR_CONTEXT_RANK.get(context, 0) + non_methyl_bonus, candidate.GetIdx(), candidate))
            ranked.sort(key=lambda item: (-item[0], item[1]))
            if not ranked or ranked[0][0] <= 10:
                continue
            # A tie between top non-identical ligands has no uniquely
            # identifiable transferable carbon in this compact v1 model.
            if len(ranked) > 1 and ranked[0][0] == ranked[1][0]:
                continue
            carbon_anchors = [ranked[0][2]]
        for anchor in carbon_anchors:
            context = classify_context(mol, anchor.GetIdx(), {handle.GetIdx()})["token"]
            handle_atoms = unique_indices([handle.GetIdx(), *(n.GetIdx() for n in handle.GetNeighbors() if n.GetIdx() != anchor.GetIdx())])
            sites.append({
                "topology": "edge", "atom_indices": unique_indices([anchor.GetIdx(), *handle_atoms]),
                "bond_indices": [bond_index(mol, anchor.GetIdx(), handle.GetIdx())],
                "signature": f"TM|{context}|{token}", "label": render_edge(context, token),
                "details": {"anchor_atom_index": anchor.GetIdx(), "handle_atom_indices": handle_atoms, "anchor_context": context, "handle_token": token},
            })
    return sites
