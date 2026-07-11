"""Pronucleophilic X-H site detection."""

from __future__ import annotations

from typing import Any, Dict, List

from ..context import classify_neighbor_contexts
from ..labels import render_xh


def _derived_family(center: str, h_count: int, contexts: List[str], aromatic: bool) -> str:
    if aromatic and center == "N": return "aromatic_nh"
    if center == "Csp": return "terminal_alkyne"
    if center == "O": return "phenol" if contexts == ["Ar"] else "alcohol"
    if center == "S": return "thiophenol" if contexts == ["Ar"] else "thiol"
    if "SO2R" in contexts: return "sulfonamide"
    if any(c.startswith("C(O)") for c in contexts): return "amide_like"
    if "N" in contexts: return "hydrazine"
    if center == "N":
        return "primary_amine" if h_count == 2 else "secondary_amine"
    return "pronucleophile"


def detect(mol: Any) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        h_count = int(atom.GetTotalNumHs(includeNeighbors=True))
        center = symbol
        supported = symbol in {"N", "O", "S"}
        if symbol == "C" and str(atom.GetHybridization()) == "SP" and h_count == 1:
            supported, center = True, "Csp"
        if not supported or h_count < 1 or atom.GetFormalCharge() != 0:
            continue
        # O-H atoms belonging to a composite transfer handle are internal to
        # that handle; they are not free alcohol/phenol pronucleophile sites.
        if symbol == "O" and any(
            neighbor.GetSymbol() in {"B"} for neighbor in atom.GetNeighbors()
        ):
            continue
        contexts_data = classify_neighbor_contexts(mol, atom.GetIdx())
        # ``classify_neighbor_contexts`` already applies the taxonomy's
        # reactivity-first precedence; preserve it in the canonical signature.
        contexts = [item["token"] for item in contexts_data]
        signature = f"XH|{center}|H{h_count}|{','.join(contexts)}"
        sites.append({
            "topology": "atom", "atom_indices": [atom.GetIdx()], "bond_indices": [],
            "signature": signature, "label": render_xh(center, h_count, contexts),
            "details": {"center_atom_index": atom.GetIdx(), "center_element": center, "formal_charge": atom.GetFormalCharge(), "aromatic": atom.GetIsAromatic(), "h_count": h_count, "contexts": contexts, "derived_family": _derived_family(center, h_count, contexts, atom.GetIsAromatic())},
        })
    return sites
