"""Aromatic carbon-hydrogen reactive-site detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_context
from ..models import ReactiveSiteCandidate
from ..patterns import MatchIndex


def detect(mol: Any, match_index: MatchIndex) -> List[ReactiveSiteCandidate]:
    """Return each aromatic carbon bearing an implicit or explicit hydrogen."""
    sites: List[ReactiveSiteCandidate] = []
    candidate_centers = match_index.role_atoms("aromatic_CH", "center")
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in candidate_centers:
            continue
        context = classify_context(mol, atom.GetIdx(), match_index=match_index)
        token = "HetArH" if context.token == "HeteroAr" else "ArH"
        h_count = int(atom.GetTotalNumHs(includeNeighbors=True))
        sites.append(ReactiveSiteCandidate(
            site_type="aromatic_CH", topology="atom",
            atom_roles={"center": (atom.GetIdx(),)},
            atom_indices=(atom.GetIdx(),), bond_indices=(),
            canonical_signature=f"CH|{token}",
            render_kind="named_handle",
            render_data={"template_id": "aromatic_ch", "context": context.token},
            matched_patterns=("aromatic_carbon_hydrogen",),
            details={
                "handle_token": token, "ring_context": context.token,
                "h_count": h_count, "aromatic": True,
            },
            context_records=(context,), availability="available",
        ))
    return sites


__all__ = ["detect"]
