"""Explicit anionic nucleophile detection."""

from __future__ import annotations

from typing import Any, List

from ..context import classify_neighbor_contexts
from ..models import SiteCandidate
from ..patterns import MatchIndex


_SUPPORTED_ELEMENTS = {"C", "N", "O", "S"}


def _derived_family(element: str, contexts: List[str]) -> str:
    if element == "S":
        return "thioacetate" if "C(O)R" in contexts else "thiolate"
    if element == "O":
        return (
            "carboxylate"
            if any(value.startswith("C(O)") for value in contexts)
            else "oxide_or_alkoxide"
        )
    if element == "N":
        return (
            "amide_anion"
            if any(value.startswith("C(O)") for value in contexts)
            else "nitrogen_anion"
        )
    if element == "C":
        return (
            "stabilized_carbanion"
            if contexts
            else "carbanion"
        )
    return "anionic_nucleophile"


def detect(mol: Any, match_index: MatchIndex) -> List[SiteCandidate]:
    """Detect explicit C/N/O/S anions without inferring deprotonation or salts."""
    sites: List[SiteCandidate] = []
    candidate_centers = match_index.role_atoms("nucleophile_anion", "center")
    for atom_index in sorted(candidate_centers):
        atom = mol.GetAtomWithIdx(atom_index)
        element = atom.GetSymbol()
        if element not in _SUPPORTED_ELEMENTS or atom.GetFormalCharge() >= 0:
            continue
        if element == "N" and any(
            neighbor.GetSymbol() == "N"
            and neighbor.GetFormalCharge() > 0
            for neighbor in atom.GetNeighbors()
        ):
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
                    f"NU-|{element}|{atom.GetFormalCharge()}|"
                    f"{','.join(contexts)}"
                ),
                render_kind="anion",
                render_data={
                    "center": element,
                    "contexts": contexts,
                    "charge": atom.GetFormalCharge(),
                },
                matched_patterns=tuple(
                    sorted(str(pattern["id"]) for pattern in pattern_matches)
                ),
                details={
                    "center_element": element,
                    "center_token": element,
                    "formal_charge": atom.GetFormalCharge(),
                    "h_count": int(atom.GetTotalNumHs(includeNeighbors=True)),
                    "contexts": contexts,
                    "derived_family": _derived_family(element, contexts),
                },
                context_records=tuple(context_records),
                availability="ionic",
            )
        )
    return sites


__all__ = ["detect"]
