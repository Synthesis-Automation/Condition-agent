"""Vicinal carbon pairs supporting conservative beta elimination."""

from __future__ import annotations

from typing import Any, List

from ..models import ReactiveSiteCandidate
from ..patterns import MatchIndex


def detect(mol: Any, match_index: MatchIndex) -> List[ReactiveSiteCandidate]:
    """Return C-X/C-H vicinal motifs with explicit atom provenance."""
    sites: List[ReactiveSiteCandidate] = []
    seen: set[tuple[int, int, int]] = set()
    candidates = match_index.role_atoms("eliminable_pair", "endpoint_a")
    for atom_index in sorted(candidates):
        for definition, roles in match_index.role_matches_for_atom(
            "eliminable_pair", "endpoint_a", atom_index
        ):
            endpoint_a_values = roles.get("endpoint_a") or ()
            endpoint_b_values = roles.get("endpoint_b") or ()
            departing_values = roles.get("departing_a") or ()
            if not (
                len(endpoint_a_values) == 1
                and len(endpoint_b_values) == 1
                and len(departing_values) == 1
            ):
                continue
            endpoint_a = int(endpoint_a_values[0])
            endpoint_b = int(endpoint_b_values[0])
            departing_a = int(departing_values[0])
            key = (endpoint_a, endpoint_b, departing_a)
            if key in seen:
                continue
            seen.add(key)
            backbone = mol.GetBondBetweenAtoms(endpoint_a, endpoint_b)
            leaving = mol.GetBondBetweenAtoms(endpoint_a, departing_a)
            if backbone is None or leaving is None:
                continue
            hydrogen_count = int(
                mol.GetAtomWithIdx(endpoint_b).GetTotalNumHs(includeNeighbors=True)
            )
            if hydrogen_count < 1:
                continue
            departing_element = mol.GetAtomWithIdx(departing_a).GetSymbol()
            sites.append(
                ReactiveSiteCandidate(
                    site_type="eliminable_pair",
                    topology="bond",
                    atom_roles={
                        "endpoint_a": (endpoint_a,),
                        "endpoint_b": (endpoint_b,),
                        "departing_a": (departing_a,),
                        "hydrogen_carrier_b": (endpoint_b,),
                    },
                    atom_indices=(endpoint_a, endpoint_b, departing_a),
                    bond_indices=(int(backbone.GetIdx()), int(leaving.GetIdx())),
                    canonical_signature=(
                        f"EP|C-C|{departing_element}-H|"
                        f"H{hydrogen_count}"
                    ),
                    render_kind="named_handle",
                    render_data={"template_id": "eliminable_pair"},
                    matched_patterns=(str(definition["id"]),),
                    details={
                        "handle_token": "BetaHalideH",
                        "departing_element": departing_element,
                        "backbone_order": "SINGLE",
                        "hydrogen_count_b": hydrogen_count,
                    },
                    availability="available",
                )
            )
    return sites


__all__ = ["detect"]
