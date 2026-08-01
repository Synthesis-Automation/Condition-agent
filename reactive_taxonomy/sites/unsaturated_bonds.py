"""Carbon-carbon unsaturated-bond reactive-site detection."""

from __future__ import annotations

from typing import Any, List

from ..models import ReactiveSiteCandidate
from ..patterns import MatchIndex


def _carbon_substitution_degree(atom: Any, other_index: int) -> int:
    return sum(
        neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != other_index
        for neighbor in atom.GetNeighbors()
    )


def _heavy_substitution_degree(atom: Any, other_index: int) -> int:
    return sum(
        neighbor.GetAtomicNum() > 1 and neighbor.GetIdx() != other_index
        for neighbor in atom.GetNeighbors()
    )


def _stereochemistry(bond: Any) -> str:
    value = str(bond.GetStereo()).upper()
    if value.endswith("STEREOE"):
        return "E"
    if value.endswith("STEREOZ"):
        return "Z"
    return ""


def detect(mol: Any, match_index: MatchIndex) -> List[ReactiveSiteCandidate]:
    """Return localized C=C, C#C, C=N, and C#N bond-capacity sites."""
    sites: List[ReactiveSiteCandidate] = []
    candidate_endpoint_a = match_index.role_atoms("unsaturated_bond", "endpoint_a")
    candidate_endpoint_b = match_index.role_atoms("unsaturated_bond", "endpoint_b")
    for bond in mol.GetBonds():
        order = float(bond.GetBondTypeAsDouble())
        if order not in {2.0, 3.0} or bond.GetIsAromatic():
            continue
        left, right = bond.GetBeginAtom(), bond.GetEndAtom()
        matched_endpoints = (
            left.GetIdx() in candidate_endpoint_a
            and right.GetIdx() in candidate_endpoint_b
        ) or (
            right.GetIdx() in candidate_endpoint_a
            and left.GetIdx() in candidate_endpoint_b
        )
        if not matched_endpoints:
            continue
        if order == 3.0 and {left.GetAtomicNum(), right.GetAtomicNum()} == {6, 7}:
            carbon = left if left.GetAtomicNum() == 6 else right
            nitrogen = right if carbon is left else left
            if (
                carbon.GetIdx() not in candidate_endpoint_a
                or nitrogen.GetIdx() not in candidate_endpoint_b
            ):
                continue
            attachment = next(
                (
                    neighbor
                    for neighbor in carbon.GetNeighbors()
                    if neighbor.GetIdx() != nitrogen.GetIdx()
                    and neighbor.GetAtomicNum() > 1
                ),
                None,
            )
            if attachment is None:
                continue
            patterns = tuple(
                definition["id"]
                for definition in match_index.patterns_for_atom(
                    "unsaturated_bond", "carbon_endpoint", carbon.GetIdx()
                )
            )
            sites.append(
                ReactiveSiteCandidate(
                    site_type="unsaturated_bond",
                    topology="bond",
                    atom_roles={
                        "endpoint_a": (carbon.GetIdx(),),
                        "endpoint_b": (nitrogen.GetIdx(),),
                        "carbon_endpoint": (carbon.GetIdx(),),
                        "nitrogen_endpoint": (nitrogen.GetIdx(),),
                        "attachment": (attachment.GetIdx(),),
                    },
                    atom_indices=(carbon.GetIdx(), nitrogen.GetIdx()),
                    bond_indices=(bond.GetIdx(),),
                    canonical_signature="PI|Nitrile",
                    render_kind="named_handle",
                    render_data={"template_id": "nitrile"},
                    matched_patterns=patterns,
                    details={
                        "handle_token": "Nitrile",
                        "bond_order": 3,
                        "endpoint_elements": ["C", "N"],
                        "attachment_element": attachment.GetSymbol(),
                        "electrophilic_endpoint_atom_index": carbon.GetIdx(),
                        "reaction_modes": ["addition", "reduction"],
                    },
                    availability="available",
                )
            )
            continue
        if order == 2.0 and {left.GetAtomicNum(), right.GetAtomicNum()} == {6, 7}:
            carbon = left if left.GetAtomicNum() == 6 else right
            nitrogen = right if carbon is left else left
            if (
                carbon.GetIdx() not in candidate_endpoint_a
                or nitrogen.GetIdx() not in candidate_endpoint_b
            ):
                continue
            patterns = tuple(
                definition["id"]
                for definition in match_index.patterns_for_atom(
                    "unsaturated_bond",
                    "carbon_endpoint",
                    carbon.GetIdx(),
                )
                if definition.get("id")
                == "carbon_nitrogen_polarized_double_bond"
            )
            if not patterns:
                continue
            sites.append(
                ReactiveSiteCandidate(
                    site_type="unsaturated_bond",
                    topology="bond",
                    atom_roles={
                        "endpoint_a": (carbon.GetIdx(),),
                        "endpoint_b": (nitrogen.GetIdx(),),
                        "carbon_endpoint": (carbon.GetIdx(),),
                        "nitrogen_endpoint": (nitrogen.GetIdx(),),
                    },
                    atom_indices=(carbon.GetIdx(), nitrogen.GetIdx()),
                    bond_indices=(bond.GetIdx(),),
                    canonical_signature="PI|PolarizedC=N",
                    render_kind="named_handle",
                    render_data={"template_id": "polarized_c_n"},
                    matched_patterns=patterns,
                    details={
                        "handle_token": "PolarizedC=N",
                        "bond_order": 2,
                        "endpoint_elements": ["C", "N"],
                        "electrophilic_endpoint_atom_index": carbon.GetIdx(),
                        "nucleophilic_endpoint_atom_index": nitrogen.GetIdx(),
                        "reaction_modes": [
                            "addition",
                            "reduction",
                            "hydrolysis",
                        ],
                    },
                    availability="available",
                )
            )
            continue
        if left.GetAtomicNum() != 6 or right.GetAtomicNum() != 6:
            continue
        token = "Alkene" if order == 2.0 else "Alkyne"
        pattern_id = "carbon_carbon_alkene" if order == 2.0 else "carbon_carbon_alkyne"
        endpoint_records = [
            {
                "atom": atom,
                "h_count": int(atom.GetTotalNumHs(includeNeighbors=True)),
                "heavy_substituents": _heavy_substitution_degree(
                    atom, other.GetIdx()
                ),
                "carbon_substituents": _carbon_substitution_degree(
                    atom, other.GetIdx()
                ),
            }
            for atom, other in ((left, right), (right, left))
        ]
        endpoint_records.sort(
            key=lambda record: (
                -record["h_count"] if order == 2.0 else record["h_count"],
                record["heavy_substituents"],
                record["atom"].GetIdx(),
            )
        )
        endpoint_a = int(endpoint_records[0]["atom"].GetIdx())
        endpoint_b = int(endpoint_records[1]["atom"].GetIdx())
        endpoint_h_counts = [
            int(record["h_count"]) for record in endpoint_records
        ]
        endpoint_substituent_counts = [
            int(record["heavy_substituents"]) for record in endpoint_records
        ]
        substitution = sum(
            int(record["carbon_substituents"]) for record in endpoint_records
        )
        heavy_substitution = sum(endpoint_substituent_counts)
        stereo = _stereochemistry(bond)
        sites.append(ReactiveSiteCandidate(
            site_type="unsaturated_bond", topology="bond",
            atom_roles={"endpoint_a": (endpoint_a,), "endpoint_b": (endpoint_b,)},
            atom_indices=(endpoint_a, endpoint_b), bond_indices=(bond.GetIdx(),),
            canonical_signature=f"PI|{token}",
            render_kind="unsaturated_bond",
            render_data={
                "bond_order": int(order),
                "endpoint_h_counts": endpoint_h_counts,
                "endpoint_substituent_counts": endpoint_substituent_counts,
                "stereochemistry": stereo,
            },
            matched_patterns=(pattern_id,),
            details={
                "handle_token": token, "bond_order": int(order),
                "substitution_degree": substitution,
                "heavy_atom_substitution_degree": heavy_substitution,
                "endpoint_h_counts": endpoint_h_counts,
                "endpoint_substituent_counts": endpoint_substituent_counts,
                "stereochemistry": stereo or None,
            },
            availability="available",
        ))
    return sites


__all__ = ["detect"]
