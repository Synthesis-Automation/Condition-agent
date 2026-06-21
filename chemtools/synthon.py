"""Taxonomy-driven synthon classification and role pairing utilities."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .taxonomy import loader as taxonomy_loader


@dataclass(frozen=True)
class SynthonDefinition:
    """A synthon class defined in taxonomy."""

    id: str
    name: str
    role: str
    motif_sets: Tuple[str, ...]
    members: Tuple[str, ...]
    priority: int


@dataclass(frozen=True)
class SynthonAssignment:
    """Synthon assignment for a single reactant."""

    synthon_id: str
    role: str
    matched_motifs: Tuple[str, ...]
    priority: int


def _coerce_str_list(values: Any) -> Tuple[str, ...]:
    if not values:
        return tuple()
    if isinstance(values, str):
        value = values.strip()
        return (value,) if value else tuple()
    if not isinstance(values, list):
        return tuple()
    out: List[str] = []
    for value in values:
        if not isinstance(value, str):
            continue
        token = value.strip()
        if token:
            out.append(token)
    return tuple(out)


def _to_priority(raw: Any, default: int = 0) -> int:
    try:
        value = int(raw)
    except Exception:
        value = default
    return max(0, value)


def _expand_motif_sets(
    motif_sets: Iterable[str],
    all_sets: Mapping[str, Any],
) -> Tuple[str, ...]:
    expanded: List[str] = []
    for set_id in motif_sets:
        entry = all_sets.get(set_id)
        if not isinstance(entry, dict):
            continue
        members = _coerce_str_list(entry.get("members"))
        expanded.extend(members)
    return tuple(sorted(set(expanded)))


@lru_cache(maxsize=1)
def load_synthon_definitions() -> Dict[str, SynthonDefinition]:
    """Load and normalize synthon definitions from taxonomy."""
    payload = taxonomy_loader.load_synthons()
    if not isinstance(payload, dict):
        return {}

    compound_logic = taxonomy_loader.load_compound_logic()
    motif_sets = {}
    if isinstance(compound_logic, dict):
        raw_sets = compound_logic.get("motif_sets")
        if isinstance(raw_sets, dict):
            motif_sets = raw_sets

    out: Dict[str, SynthonDefinition] = {}
    for entry in payload.get("synthons", []) or []:
        if not isinstance(entry, dict):
            continue
        synthon_id = str(entry.get("id") or "").strip()
        role = str(entry.get("role") or "").strip().lower()
        if not synthon_id or role not in {"electrophile", "nucleophile"}:
            continue
        motif_set_ids = _coerce_str_list(entry.get("motif_sets"))
        members = _expand_motif_sets(motif_set_ids, motif_sets)
        if not members:
            continue
        out[synthon_id] = SynthonDefinition(
            id=synthon_id,
            name=str(entry.get("name") or synthon_id),
            role=role,
            motif_sets=motif_set_ids,
            members=members,
            priority=_to_priority(entry.get("priority"), default=0),
        )
    return out


def _get_reactant_category_matches(smiles: str) -> Tuple[str, ...]:
    """
    Return taxonomy reactant categories for a SMILES string.

    Uses the calculable feature engine directly to avoid eager imports that can
    create circular dependency paths during package initialization.
    """
    try:
        from .featurizers import calculable as _calc
    except Exception:
        return tuple()
    try:
        reactant_features = _calc.get_reactant_type_features(smiles)
    except Exception:
        return tuple()
    if not isinstance(reactant_features, dict):
        return tuple()
    categories = reactant_features.get("categories", []) or []
    if not isinstance(categories, list):
        return tuple()
    out = sorted({str(category) for category in categories if str(category).strip()})
    return tuple(out)


@lru_cache(maxsize=8192)
def classify_reactant_synthons(smiles: str) -> Tuple[SynthonAssignment, ...]:
    """Return all matched synthons for a reactant SMILES."""
    text = str(smiles or "").strip()
    if not text:
        return tuple()

    definitions = load_synthon_definitions()
    if not definitions:
        return tuple()

    categories = set(_get_reactant_category_matches(text))
    if not categories:
        return tuple()

    assignments: List[SynthonAssignment] = []
    for definition in definitions.values():
        matched = sorted(categories.intersection(definition.members))
        if not matched:
            continue
        assignments.append(
            SynthonAssignment(
                synthon_id=definition.id,
                role=definition.role,
                matched_motifs=tuple(matched),
                priority=definition.priority,
            )
        )

    assignments.sort(
        key=lambda item: (
            -item.priority,
            -len(item.matched_motifs),
            item.synthon_id,
        )
    )
    return tuple(assignments)


def _best_role_hit(assignments: Sequence[SynthonAssignment], role: str) -> Optional[SynthonAssignment]:
    for hit in assignments:
        if hit.role == role:
            return hit
    return None


_ROLE_FALLBACK_GROUPS: Dict[str, Tuple[str, ...]] = {
    "electrophile": (
        "acyl_halide",
        "aryl_halide",
        "vinyl_halide",
        "alkyl_halide",
        "triflate",
        "mesylate",
        "tosylate",
        "aryl_sulfonate",
    ),
    "nucleophile": (
        "boron",
        "nucleophile_n",
        "nucleophile_o",
        "nucleophile_s",
    ),
}


def _fallback_role_matches(smiles: str, role: str) -> Tuple[str, ...]:
    """Use shared functional-group detection as the chemistry-first fallback."""
    groups = _ROLE_FALLBACK_GROUPS.get(role)
    if not groups:
        return tuple()
    try:
        from chemtools.molecule import functional_groups as _functional_groups
    except Exception:
        return tuple()

    detections = _functional_groups.detect_any([smiles], group_subset=groups)
    matches: List[str] = []
    for group_name in groups:
        definition = _functional_groups.get_group_definition(group_name)
        if not definition:
            continue
        if detections.get(definition.token):
            matches.append(group_name)
    return tuple(matches)


def _best_role_candidates(
    reactants: Sequence[str],
    role: str,
) -> List[Tuple[int, int, int]]:
    candidates: List[Tuple[int, int, int]] = []
    for idx, smiles in enumerate(reactants):
        assignments = classify_reactant_synthons(smiles)
        best = _best_role_hit(assignments, role)
        if best is not None:
            # Taxonomy synthons remain authoritative. Boost their score above
            # any functional-group fallback evidence.
            candidates.append((1000 + best.priority, len(best.matched_motifs), -idx))
            continue

        fallback_matches = _fallback_role_matches(smiles, role)
        if fallback_matches:
            candidates.append((len(fallback_matches), len(fallback_matches), -idx))
    candidates.sort(reverse=True)
    return candidates


def select_electrophile_nucleophile(reactants: Sequence[str]) -> Tuple[str, str]:
    """Select electrophile and nucleophile from a reactant list.

    Selection is taxonomy-driven via synthon classes with centralized
    functional-group detection as the chemistry-first fallback.
    """
    normalized = [str(value).strip() for value in reactants if str(value).strip()]
    if not normalized:
        return "", ""
    if len(normalized) == 1:
        return normalized[0], ""

    elec_candidates = _best_role_candidates(normalized, "electrophile")
    nuc_candidates = _best_role_candidates(normalized, "nucleophile")

    if elec_candidates and nuc_candidates:
        for e_pri, e_hits, e_neg_idx in elec_candidates:
            e_idx = -e_neg_idx
            for n_pri, n_hits, n_neg_idx in nuc_candidates:
                n_idx = -n_neg_idx
                if e_idx == n_idx:
                    continue
                return normalized[e_idx], normalized[n_idx]

    if elec_candidates:
        e_idx = -elec_candidates[0][2]
        n_idx = next((idx for idx in range(len(normalized)) if idx != e_idx), None)
        return normalized[e_idx], (normalized[n_idx] if n_idx is not None else "")

    if nuc_candidates:
        n_idx = -nuc_candidates[0][2]
        e_idx = next((idx for idx in range(len(normalized)) if idx != n_idx), None)
        return (normalized[e_idx] if e_idx is not None else normalized[0]), normalized[n_idx]

    return normalized[0], normalized[1]


__all__ = [
    "SynthonDefinition",
    "SynthonAssignment",
    "load_synthon_definitions",
    "classify_reactant_synthons",
    "select_electrophile_nucleophile",
]
