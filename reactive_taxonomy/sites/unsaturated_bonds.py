"""Carbon-carbon unsaturated-bond reactive-site detection."""

from __future__ import annotations

from typing import Any, List

from ..models import SiteCandidate
from ..patterns import MatchIndex


def _substitution_degree(atom: Any, other_index: int) -> int:
    return sum(
        neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != other_index
        for neighbor in atom.GetNeighbors()
    )


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Return non-aromatic C=C and C#C bonds as two-endpoint sites."""
    sites: List[SiteCandidate] = []
    candidate_atoms = (
        match_index.role_atoms("unsaturated_bond", "endpoint_a")
        | match_index.role_atoms("unsaturated_bond", "endpoint_b")
    )
    for bond in mol.GetBonds():
        order = float(bond.GetBondTypeAsDouble())
        if order not in {2.0, 3.0} or bond.GetIsAromatic():
            continue
        left, right = bond.GetBeginAtom(), bond.GetEndAtom()
        if left.GetAtomicNum() != 6 or right.GetAtomicNum() != 6:
            continue
        if left.GetIdx() not in candidate_atoms or right.GetIdx() not in candidate_atoms:
            continue
        endpoint_a, endpoint_b = sorted((left.GetIdx(), right.GetIdx()))
        token = "Alkene" if order == 2.0 else "Alkyne"
        pattern_id = "carbon_carbon_alkene" if order == 2.0 else "carbon_carbon_alkyne"
        substitution = (
            _substitution_degree(left, right.GetIdx())
            + _substitution_degree(right, left.GetIdx())
        )
        sites.append(SiteCandidate(
            site_type="unsaturated_bond", topology="bond",
            atom_roles={"endpoint_a": (endpoint_a,), "endpoint_b": (endpoint_b,)},
            atom_indices=(endpoint_a, endpoint_b), bond_indices=(bond.GetIdx(),),
            canonical_signature=f"PI|{token}",
            render_kind="named_handle",
            render_data={"template_id": token.lower()},
            matched_patterns=(pattern_id,),
            details={
                "handle_token": token, "bond_order": int(order),
                "substitution_degree": substitution,
            },
            availability="available",
        ))
    return sites


__all__ = ["detect"]
