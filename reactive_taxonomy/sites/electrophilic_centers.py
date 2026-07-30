"""Acyl and sulfonyl electrophilic-center detection."""

from __future__ import annotations

from typing import Any, List, Optional, Tuple

from ..context import classify_context, classify_neighbor_contexts
from ..models import SiteCandidate
from ..patterns import MatchIndex
from .common import bond_index, unique_indices


def _single_bond_neighbors(atom: Any) -> List[Any]:
    return [bond.GetOtherAtom(atom) for bond in atom.GetBonds() if str(bond.GetBondType()) == "SINGLE"]


def _acyl_group(neighbors: List[Any]) -> Optional[Tuple[Any, str, str]]:
    for neighbor in neighbors:
        if neighbor.GetSymbol() in {"Cl", "Br", "F"}:
            return neighbor, neighbor.GetSymbol(), "activated"
    for neighbor in neighbors:
        if neighbor.GetSymbol() == "O":
            if neighbor.GetFormalCharge() < 0:
                return neighbor, "O-", "ionic"
            if neighbor.GetTotalNumHs(includeNeighbors=True) > 0:
                return neighbor, "OH", "latent"
            return neighbor, "OR", "ester"
    return None


def _strained_ring_sites(
    mol: Any,
    match_index: MatchIndex,
) -> List[SiteCandidate]:
    """Expose each carbon of an epoxide or aziridine as a ring-opening center."""
    sites: List[SiteCandidate] = []
    candidate_centers = match_index.role_atoms(
        "electrophilic_center",
        "center",
    )
    candidate_heteroatoms = match_index.role_atoms(
        "electrophilic_center",
        "heteroatom",
    )
    for center_index in sorted(candidate_centers):
        definitions = [
            definition
            for definition in match_index.patterns_for_atom(
                "electrophilic_center",
                "center",
                center_index,
            )
            if definition.get("center_family") == "StrainedRing"
        ]
        if not definitions:
            continue
        center = mol.GetAtomWithIdx(center_index)
        heteroatoms = [
            neighbor
            for neighbor in center.GetNeighbors()
            if neighbor.GetIdx() in candidate_heteroatoms
            and neighbor.IsInRingSize(3)
            and neighbor.GetSymbol() in {"N", "O"}
        ]
        for heteroatom in heteroatoms:
            partners = [
                neighbor
                for neighbor in heteroatom.GetNeighbors()
                if neighbor.GetIdx() != center_index
                and neighbor.GetIdx() in candidate_centers
                and neighbor.GetSymbol() == "C"
                and neighbor.IsInRingSize(3)
                and mol.GetBondBetweenAtoms(
                    center_index,
                    neighbor.GetIdx(),
                )
                is not None
            ]
            if len(partners) != 1:
                continue
            partner = partners[0]
            subtype = "epoxide" if heteroatom.GetSymbol() == "O" else "aziridine"
            pattern_id = f"{subtype}_ring"
            if not any(
                definition.get("id") == pattern_id
                for definition in definitions
            ):
                continue
            context = classify_context(
                mol,
                center_index,
                {heteroatom.GetIdx()},
                match_index=match_index,
            )
            sites.append(
                SiteCandidate(
                    site_type="electrophilic_center",
                    topology="center",
                    atom_roles={
                        "center": (center_index,),
                        "heteroatom": (heteroatom.GetIdx(),),
                        "leaving_or_activatable": (heteroatom.GetIdx(),),
                        "ring_partner": (partner.GetIdx(),),
                    },
                    atom_indices=tuple(
                        unique_indices(
                            [
                                center_index,
                                heteroatom.GetIdx(),
                                partner.GetIdx(),
                            ]
                        )
                    ),
                    bond_indices=(
                        bond_index(
                            mol,
                            center_index,
                            heteroatom.GetIdx(),
                        ),
                    ),
                    canonical_signature=(
                        f"EC|StrainedRing|{subtype}|"
                        f"{context.token}|ring_opening"
                    ),
                    render_kind="named_handle",
                    render_data={"template_id": f"{subtype}_carbon"},
                    matched_patterns=(pattern_id,),
                    details={
                        "center_family": "StrainedRing",
                        "strained_ring_type": subtype,
                        "reaction_mode": "ring_opening",
                        "anchor_context": context.token,
                        "activation_state": "strain_activated",
                    },
                    context_records=(context,),
                    availability="activated",
                )
            )
    return sites


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    sites: List[SiteCandidate] = _strained_ring_sites(mol, match_index)
    candidate_centers = match_index.role_atoms("electrophilic_center", "center")
    for center in mol.GetAtoms():
        if center.GetIdx() not in candidate_centers:
            continue
        symbol = center.GetSymbol()
        double_oxygens = [bond.GetOtherAtom(center) for bond in center.GetBonds() if str(bond.GetBondType()) == "DOUBLE" and bond.GetOtherAtom(center).GetSymbol() == "O"]
        singles = _single_bond_neighbors(center)
        if symbol == "C" and double_oxygens:
            group = _acyl_group(singles)
            if group is not None:
                leaving, token, state = group
                center_pattern_ids = {
                    str(item["id"])
                    for item in match_index.patterns_for_atom("electrophilic_center", "center", center.GetIdx())
                }
                if "acyl_anhydride" in center_pattern_ids:
                    token, state = "OC(O)R", "activated"
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
                    matched_patterns=tuple(sorted(center_pattern_ids & {"acyl_anhydride", "acyl_oxygen_or_halide"})),
                    details={
                        "center_family": "Acyl",
                        "acyl_subtype": "acid_anhydride" if "acyl_anhydride" in center_pattern_ids else state,
                        "reaction_mode": "substitution", "retained_context": context_token,
                        "leaving_or_activatable_group": token, "activation_state": state,
                    },
                    context_records=(context,) if context is not None else (),
                    availability="activated" if state == "activated" else ("latent" if state == "latent" else "conditional"),
                ))
                continue

            carbonyl_patterns = {
                str(item["id"])
                for item in match_index.patterns_for_atom("electrophilic_center", "center", center.GetIdx())
                if item.get("tokens") == ["Carbonyl"]
            }
            if not carbonyl_patterns:
                continue
            h_count = int(center.GetTotalNumHs(includeNeighbors=True))
            subtype = "ketone" if h_count == 0 else ("formaldehyde" if h_count == 2 else "aldehyde")
            contexts = classify_neighbor_contexts(mol, center.GetIdx(), {double_oxygens[0].GetIdx()}, match_index=match_index)
            context_tokens = [record.token for record in contexts]
            context_label = ",".join(context_tokens) if context_tokens else "Other"
            substituents = tuple(neighbor.GetIdx() for neighbor in singles if neighbor.GetAtomicNum() == 6)
            sites.append(SiteCandidate(
                site_type="electrophilic_center", topology="center",
                atom_roles={
                    "center": (center.GetIdx(),),
                    "heteroatom": (double_oxygens[0].GetIdx(),),
                    **({"substituents": substituents} if substituents else {}),
                },
                atom_indices=tuple(unique_indices([center.GetIdx(), double_oxygens[0].GetIdx(), *substituents])),
                bond_indices=(bond_index(mol, center.GetIdx(), double_oxygens[0].GetIdx()),),
                canonical_signature=f"EC|Carbonyl|{subtype}|{context_label}|addition",
                render_kind="named_handle", render_data={"template_id": f"carbonyl_{subtype}"},
                matched_patterns=tuple(sorted(carbonyl_patterns)),
                details={
                    "center_family": "Carbonyl", "carbonyl_subtype": subtype,
                    "reaction_mode": "addition", "contexts": context_tokens,
                    "activation_state": "intrinsic",
                },
                context_records=tuple(contexts), availability="available",
            ))
        elif symbol == "S" and len(double_oxygens) >= 2:
            leaving = next((neighbor for neighbor in singles if neighbor.GetSymbol() in {"Cl", "Br", "F"}), None)
            if leaving is None:
                continue
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
                matched_patterns=("sulfonyl_halide",),
                details={"center_family": "Sulfonyl", "reaction_mode": "substitution", "retained_context": context_token, "leaving_or_activatable_group": token, "activation_state": "activated"},
                context_records=(context,) if context is not None else (), availability="activated",
            ))
    return sites
