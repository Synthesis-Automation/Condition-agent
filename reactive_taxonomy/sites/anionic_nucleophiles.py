"""Anionic heteroatom nucleophile detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_neighbor_contexts
from ..models import SiteCandidate
from ..patterns import MatchIndex


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Detect explicitly charged sulfur nucleophiles without inferring salts."""
    sites: List[SiteCandidate] = []
    candidate_centers = match_index.role_atoms("nucleophile_anion", "center")
    for atom_index in sorted(candidate_centers):
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() != "S" or atom.GetFormalCharge() >= 0:
            continue
        context_records = classify_neighbor_contexts(
            mol,
            atom_index,
            match_index=match_index,
        )
        contexts = [record.token for record in context_records]
        pattern_matches = match_index.patterns_for_atom(
            "nucleophile_anion",
            "center",
            atom_index,
        )
        sites.append(
            SiteCandidate(
                site_type="nucleophile_anion",
                topology="atom",
                atom_roles={"center": (atom_index,)},
                atom_indices=(atom_index,),
                bond_indices=(),
                canonical_signature=(
                    f"NU-|S|{atom.GetFormalCharge()}|{','.join(contexts)}"
                ),
                render_kind="anion",
                render_data={
                    "center": "S",
                    "contexts": contexts,
                    "charge": atom.GetFormalCharge(),
                },
                matched_patterns=tuple(
                    sorted(str(pattern["id"]) for pattern in pattern_matches)
                ),
                details={
                    "center_element": "S",
                    "center_token": "S",
                    "formal_charge": atom.GetFormalCharge(),
                    "h_count": int(atom.GetTotalNumHs(includeNeighbors=True)),
                    "contexts": contexts,
                    "derived_family": (
                        "thioacetate"
                        if "C(O)R" in contexts
                        else "thiolate"
                    ),
                },
                context_records=tuple(context_records),
                availability="ionic",
            )
        )
    return sites


__all__ = ["detect"]
