"""Acyl and sulfonyl electrophilic-center detection."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from ..context import classify_context
from ..labels import render_edge
from .common import bond_index, unique_indices


def _single_bond_neighbors(mol: Any, atom: Any) -> List[Any]:
    return [b.GetOtherAtom(atom) for b in atom.GetBonds() if str(b.GetBondType()) == "SINGLE"]


def _acyl_group(neighbors: List[Any]) -> Optional[Tuple[Any, str, str]]:
    for n in neighbors:
        if n.GetSymbol() in {"Cl", "Br", "F"}: return n, n.GetSymbol(), "activated"
    for n in neighbors:
        if n.GetSymbol() == "O":
            if n.GetFormalCharge() < 0: return n, "O-", "ionic"
            if n.GetTotalNumHs(includeNeighbors=True) > 0: return n, "OH", "latent"
            return n, "OR", "ester"
    return None


def detect(mol: Any) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    for center in mol.GetAtoms():
        symbol = center.GetSymbol()
        double_oxygens = [b.GetOtherAtom(center) for b in center.GetBonds() if str(b.GetBondType()) == "DOUBLE" and b.GetOtherAtom(center).GetSymbol() == "O"]
        singles = _single_bond_neighbors(mol, center)
        if symbol == "C" and double_oxygens:
            group = _acyl_group(singles)
            if group is None: continue
            leaving, token, state = group
            retained = next((n for n in singles if n.GetIdx() != leaving.GetIdx()), None)
            context = classify_context(mol, retained.GetIdx(), {center.GetIdx()})["token"] if retained is not None else "Other"
            label = f"{render_edge(context, 'C(O)' + token)}"
            atoms = unique_indices([center.GetIdx(), double_oxygens[0].GetIdx(), leaving.GetIdx(), *([retained.GetIdx()] if retained else [])])
            sites.append({"topology": "center", "atom_indices": atoms, "bond_indices": [bond_index(mol, center.GetIdx(), leaving.GetIdx())], "signature": f"EC|Acyl|{context}|{token}|{state}", "label": label, "details": {"center_family": "Acyl", "center_atom_index": center.GetIdx(), "multiple_bond_atom_indices": [double_oxygens[0].GetIdx()], "retained_context": context, "leaving_or_activatable_group": token, "activation_state": state}})
        elif symbol == "S" and len(double_oxygens) >= 2:
            leaving = next((n for n in singles if n.GetSymbol() in {"Cl", "Br", "F"}), None)
            if leaving is None: continue
            retained = next((n for n in singles if n.GetIdx() != leaving.GetIdx()), None)
            context = classify_context(mol, retained.GetIdx(), {center.GetIdx()})["token"] if retained is not None else "Other"
            token = leaving.GetSymbol()
            sites.append({"topology": "center", "atom_indices": unique_indices([center.GetIdx(), leaving.GetIdx(), *(o.GetIdx() for o in double_oxygens), *([retained.GetIdx()] if retained else [])]), "bond_indices": [bond_index(mol, center.GetIdx(), leaving.GetIdx())], "signature": f"EC|Sulfonyl|{context}|{token}|activated", "label": f"{render_edge(context, 'S(O)2' + token)}", "details": {"center_family": "Sulfonyl", "center_atom_index": center.GetIdx(), "multiple_bond_atom_indices": [o.GetIdx() for o in double_oxygens], "retained_context": context, "leaving_or_activatable_group": token, "activation_state": "activated"}})
    return sites
