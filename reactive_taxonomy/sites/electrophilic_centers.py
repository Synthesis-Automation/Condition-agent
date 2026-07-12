"""Acyl and sulfonyl electrophilic-center detection."""

from __future__ import annotations

from typing import Any, List, Optional, Tuple

from ..context import classify_context
from ..models import SiteCandidate
from ..patterns import MatchIndex
from .common import bond_index, unique_indices


def _single_bond_neighbors(atom: Any) -> List[Any]:
    return [bond.GetOtherAtom(atom) for bond in atom.GetBonds() if str(bond.GetBondType()) == "SINGLE"]


def _acyl_group(neighbors: List[Any]) -> Optional[Tuple[Any, str, str]]:
    for neighbor in neighbors:
        if neighbor.GetSymbol() in {"Cl", "Br", "F"}: return neighbor, neighbor.GetSymbol(), "activated"
    for neighbor in neighbors:
        if neighbor.GetSymbol() == "O":
            if neighbor.GetFormalCharge() < 0: return neighbor, "O-", "ionic"
            if neighbor.GetTotalNumHs(includeNeighbors=True) > 0: return neighbor, "OH", "latent"
            return neighbor, "OR", "ester"
    return None


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    sites: List[SiteCandidate] = []
    candidate_centers = match_index.role_atoms("electrophilic_center", "center")
    for center in mol.GetAtoms():
        if center.GetIdx() not in candidate_centers:
            continue
        symbol = center.GetSymbol()
        double_oxygens = [bond.GetOtherAtom(center) for bond in center.GetBonds() if str(bond.GetBondType()) == "DOUBLE" and bond.GetOtherAtom(center).GetSymbol() == "O"]
        singles = _single_bond_neighbors(center)
        if symbol == "C" and double_oxygens:
            group = _acyl_group(singles)
            if group is None: continue
            leaving, token, state = group
            retained = next((neighbor for neighbor in singles if neighbor.GetIdx() != leaving.GetIdx()), None)
            context = classify_context(mol, retained.GetIdx(), {center.GetIdx()}, match_index=match_index) if retained is not None else None
            context_token = context.token if context is not None else "Other"
            atoms = tuple(unique_indices([center.GetIdx(), double_oxygens[0].GetIdx(), leaving.GetIdx(), *([retained.GetIdx()] if retained else [])]))
            sites.append(SiteCandidate(
                site_type="electrophilic_center", topology="center",
                atom_roles={"center": (center.GetIdx(),), "oxo": (double_oxygens[0].GetIdx(),), "leaving_or_activatable": (leaving.GetIdx(),), **({"retained": (retained.GetIdx(),)} if retained else {})},
                atom_indices=atoms, bond_indices=(bond_index(mol, center.GetIdx(), leaving.GetIdx()),),
                canonical_signature=f"EC|Acyl|{context_token}|{token}|{state}", render_kind="edge",
                render_data={"context": context_token, "handle": f"C(O){token}"},
                matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("electrophilic_center", "center", center.GetIdx())),
                details={"center_family": "Acyl", "retained_context": context_token, "leaving_or_activatable_group": token, "activation_state": state},
                context_records=(context,) if context is not None else (),
                availability="activated" if state == "activated" else ("latent" if state == "latent" else "conditional"),
            ))
        elif symbol == "S" and len(double_oxygens) >= 2:
            leaving = next((neighbor for neighbor in singles if neighbor.GetSymbol() in {"Cl", "Br", "F"}), None)
            if leaving is None: continue
            retained = next((neighbor for neighbor in singles if neighbor.GetIdx() != leaving.GetIdx()), None)
            context = classify_context(mol, retained.GetIdx(), {center.GetIdx()}, match_index=match_index) if retained is not None else None
            context_token = context.token if context is not None else "Other"
            token = leaving.GetSymbol()
            sites.append(SiteCandidate(
                site_type="electrophilic_center", topology="center",
                atom_roles={"center": (center.GetIdx(),), "oxo": tuple(oxygen.GetIdx() for oxygen in double_oxygens), "leaving_or_activatable": (leaving.GetIdx(),), **({"retained": (retained.GetIdx(),)} if retained else {})},
                atom_indices=tuple(unique_indices([center.GetIdx(), leaving.GetIdx(), *(oxygen.GetIdx() for oxygen in double_oxygens), *([retained.GetIdx()] if retained else [])])),
                bond_indices=(bond_index(mol, center.GetIdx(), leaving.GetIdx()),),
                canonical_signature=f"EC|Sulfonyl|{context_token}|{token}|activated", render_kind="edge",
                render_data={"context": context_token, "handle": f"S(O)2{token}"},
                matched_patterns=tuple(item["id"] for item in match_index.patterns_for_atom("electrophilic_center", "center", center.GetIdx())),
                details={"center_family": "Sulfonyl", "retained_context": context_token, "leaving_or_activatable_group": token, "activation_state": "activated"},
                context_records=(context,) if context is not None else (), availability="activated",
            ))
    return sites
