"""Bond-localized organic heteroatom-pair handle detection."""

from __future__ import annotations

from typing import Any, List

from ..models import SiteCandidate
from ..patterns import MatchIndex

_KINDS = {
    ("N", "N", 2): (
        "Azo",
        "HB|Azo",
        "azo_bond",
        ["isomerization", "reduction"],
    ),
    ("S", "S", 1): (
        "Disulfide",
        "HB|Disulfide",
        "disulfide_bond",
        ["exchange", "reduction"],
    ),
    ("O", "O", 1): (
        "Peroxide",
        "HB|Peroxide",
        "peroxide_bond",
        ["homolysis", "reduction"],
    ),
}


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Return organic R–N=N–R, R–S–S–R, and R–O–O–R bonds."""
    endpoint_a_candidates = match_index.role_atoms("heteroatom_bond", "endpoint_a")
    endpoint_b_candidates = match_index.role_atoms("heteroatom_bond", "endpoint_b")
    attachment_a_candidates = match_index.role_atoms(
        "heteroatom_bond", "attachment_a"
    )
    attachment_b_candidates = match_index.role_atoms(
        "heteroatom_bond", "attachment_b"
    )
    sites: List[SiteCandidate] = []
    for bond in mol.GetBonds():
        left = bond.GetBeginAtom()
        right = bond.GetEndAtom()
        key = (
            left.GetSymbol(),
            right.GetSymbol(),
            int(round(float(bond.GetBondTypeAsDouble()))),
        )
        kind = _KINDS.get(key)
        if kind is None:
            continue
        if (
            left.GetIdx() in endpoint_a_candidates
            and right.GetIdx() in endpoint_b_candidates
        ):
            endpoint_a, endpoint_b = left, right
        elif (
            right.GetIdx() in endpoint_a_candidates
            and left.GetIdx() in endpoint_b_candidates
        ):
            endpoint_a, endpoint_b = right, left
        else:
            continue
        attachments_a = [
            neighbor
            for neighbor in endpoint_a.GetNeighbors()
            if neighbor.GetIdx() != endpoint_b.GetIdx()
            and neighbor.GetIdx() in attachment_a_candidates
        ]
        attachments_b = [
            neighbor
            for neighbor in endpoint_b.GetNeighbors()
            if neighbor.GetIdx() != endpoint_a.GetIdx()
            and neighbor.GetIdx() in attachment_b_candidates
        ]
        if len(attachments_a) != 1 or len(attachments_b) != 1:
            continue
        attachment_a = attachments_a[0]
        attachment_b = attachments_b[0]
        token, signature, template_id, reaction_modes = kind
        atoms = (
            attachment_a.GetIdx(),
            endpoint_a.GetIdx(),
            endpoint_b.GetIdx(),
            attachment_b.GetIdx(),
        )
        patterns = tuple(
            definition["id"]
            for definition in match_index.patterns_for_atom(
                "heteroatom_bond", "endpoint_a", endpoint_a.GetIdx()
            )
        )
        sites.append(
            SiteCandidate(
                site_type="heteroatom_bond",
                topology="bond",
                atom_roles={
                    "attachment_a": (attachment_a.GetIdx(),),
                    "endpoint_a": (endpoint_a.GetIdx(),),
                    "endpoint_b": (endpoint_b.GetIdx(),),
                    "attachment_b": (attachment_b.GetIdx(),),
                },
                atom_indices=atoms,
                bond_indices=(
                    mol.GetBondBetweenAtoms(
                        attachment_a.GetIdx(), endpoint_a.GetIdx()
                    ).GetIdx(),
                    bond.GetIdx(),
                    mol.GetBondBetweenAtoms(
                        endpoint_b.GetIdx(), attachment_b.GetIdx()
                    ).GetIdx(),
                ),
                canonical_signature=signature,
                render_kind="named_handle",
                render_data={"template_id": template_id},
                matched_patterns=patterns,
                details={
                    "handle_token": token,
                    "central_bond_order": int(
                        round(float(bond.GetBondTypeAsDouble()))
                    ),
                    "endpoint_elements": [
                        endpoint_a.GetSymbol(), endpoint_b.GetSymbol()
                    ],
                    "attachment_elements": [
                        attachment_a.GetSymbol(), attachment_b.GetSymbol()
                    ],
                    "reaction_modes": reaction_modes,
                },
                availability="available",
            )
        )
    return sites


__all__ = ["detect"]
