"""Type-agnostic hierarchical retrieval over generic reaction signatures."""

from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple

from .compatibility import CompatibilityAssessment, filter_compatible_precedents
from .generic_indexing import GenericIndexedReaction, GenericReactionIndex

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_retrieval.v1.json"


@lru_cache(maxsize=1)
def load_generic_retrieval_rules() -> Dict[str, Any]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


def _positions(mapping: Mapping[str, Tuple[int, ...]], key: Any) -> set[int]:
    return set(mapping.get(str(key or ""), ()))


def _compatible_edit_positions(
    signature: Mapping[str, Any], index: GenericReactionIndex
) -> set[int]:
    """Enforce the net bond-edit compatibility gate before similarity."""
    return _positions(index.bond_edits, signature.get("bond_edit_signature_key"))


def _candidate_levels(
    signature: Mapping[str, Any], index: GenericReactionIndex
) -> list[tuple[str, set[int]]]:
    compatible = _compatible_edit_positions(signature, index)
    if not compatible:
        return []
    rules = load_generic_retrieval_rules()
    levels: list[tuple[str, set[int]]] = [
        (
            "exact_signature",
            _positions(index.exact, signature.get("exact_signature_key")) & compatible,
        ),
        (
            "handle_signature",
            _positions(index.handles, signature.get("handle_signature_key"))
            & compatible,
        ),
    ]
    family = str(signature.get("named_family") or "")
    family_confidence = float(signature.get("family_confidence") or 0.0)
    threshold = float(rules["high_confidence_family_threshold"])
    family_positions = (
        _positions(index.families, family) & compatible
        if family and family_confidence >= threshold
        else set()
    )
    levels.extend(
        (
            ("named_family", family_positions),
            (
                "transformation_signature",
                _positions(
                    index.transformations,
                    signature.get("transformation_signature_key"),
                )
                & compatible,
            ),
            ("bond_edit_signature", compatible),
            (
                "environment_signature",
                _positions(
                    index.environments,
                    signature.get("environment_signature_key"),
                )
                & compatible,
            ),
        )
    )
    return levels


def _minimum_support(minimum_pool_size: int | None) -> int:
    rules = load_generic_retrieval_rules()
    minimum = (
        int(minimum_pool_size)
        if minimum_pool_size is not None
        else int(rules["minimum_pool_size"])
    )
    if minimum < 1:
        raise ValueError("minimum_pool_size must be positive")
    return minimum


def retrieve_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> Tuple[str, Tuple[GenericIndexedReaction, ...]]:
    """Select the first adequately supported chemistry-compatible tier."""
    minimum = _minimum_support(minimum_pool_size)
    compatible = _compatible_edit_positions(signature, index)
    if not compatible:
        return "no_compatible_bond_edit", ()
    levels = _candidate_levels(signature, index)
    fallback: tuple[str, set[int]] | None = None
    for level, positions in levels:
        if not positions:
            continue
        if fallback is None or len(positions) > len(fallback[1]):
            fallback = (level, positions)
        if len(positions) >= minimum:
            return level, index.select(sorted(positions))
    if fallback is None:
        return "no_compatible_precedent", ()
    level, positions = fallback
    return f"{level}_limited_support", index.select(sorted(positions))


def retrieve_compatible_generic_pool(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    *,
    minimum_pool_size: int | None = None,
) -> tuple[
    str,
    tuple[tuple[GenericIndexedReaction, CompatibilityAssessment], ...],
    int,
    int,
]:
    """Apply hard recipe compatibility at every level before support checks."""
    minimum = _minimum_support(minimum_pool_size)
    levels = _candidate_levels(signature, index)
    if not levels:
        return "no_compatible_bond_edit", (), 0, 0
    fallback = None
    for level, positions in levels:
        if not positions:
            continue
        rows = index.select(sorted(positions))
        accepted, excluded = filter_compatible_precedents(signature, rows)
        if not accepted:
            continue
        candidate = (level, accepted, len(rows), len(excluded))
        if fallback is None or len(accepted) > len(fallback[1]):
            fallback = candidate
        if len(accepted) >= minimum:
            return candidate
    if fallback is None:
        bond_rows = index.select(sorted(_compatible_edit_positions(signature, index)))
        _, excluded = filter_compatible_precedents(signature, bond_rows)
        return "no_compatible_condition_precedent", (), len(bond_rows), len(excluded)
    level, accepted, raw_count, excluded_count = fallback
    return f"{level}_limited_support", accepted, raw_count, excluded_count


def _jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a, b = set(left), set(right)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _multiset_jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_counts = Counter(left)
    right_counts = Counter(right)
    if not left_counts or not right_counts:
        return 0.0
    keys = set(left_counts).union(right_counts)
    intersection = sum(min(left_counts[key], right_counts[key]) for key in keys)
    union = sum(max(left_counts[key], right_counts[key]) for key in keys)
    return intersection / union if union else 0.0


def _partner_tokens(signature: Mapping[str, Any], field: str) -> Tuple[str, ...]:
    return tuple(
        sorted(
            str(value)
            for partner in signature.get("partners") or ()
            for value in partner.get(field) or ()
        )
    )


def _environment_tokens(signature: Mapping[str, Any]) -> Tuple[str, ...]:
    tokens = []
    for partner in signature.get("partners") or ():
        for category in ("steric", "electronic"):
            value = (partner.get(category) or {}).get("class")
            if value:
                tokens.append(f"{category}:{value}")
        for group in partner.get("nearby_groups") or ():
            group_id = group.get("group_id")
            if group_id:
                tokens.append(f"nearby:{group_id}")
    return tuple(sorted(tokens))


def _spectator_tokens(signature: Mapping[str, Any]) -> Tuple[str, ...]:
    return tuple(
        sorted(
            str(group.get("group_id"))
            for group in signature.get("spectator_groups") or ()
            if group.get("group_id")
        )
    )


def _event_tokens(signature: Mapping[str, Any]) -> Tuple[str, ...]:
    """Return multiplicity-preserving, environment-neutral event tokens."""
    tokens = []
    for event in signature.get("events") or ():
        if not isinstance(event, Mapping):
            continue
        token = {
            "formed": tuple(event.get("formed_bond_types") or ()),
            "broken": tuple(event.get("broken_bond_types") or ()),
            "order_changes": tuple(event.get("order_changes") or ()),
            "hydrogen_changes": tuple(event.get("hydrogen_changes") or ()),
            "transformation_class": str(event.get("transformation_class") or ""),
            "reaction_scope": str(
                (event.get("topology") or {}).get("reaction_scope") or ""
            ),
        }
        tokens.append(json.dumps(token, sort_keys=True, separators=(",", ":")))
    return tuple(sorted(tokens))


def reaction_scope(signature: Mapping[str, Any]) -> str:
    """Return the structurally observed reaction scope, when available."""
    topology = signature.get("topology") or {}
    if not isinstance(topology, Mapping):
        return ""
    return str(topology.get("reaction_scope") or "")


def _reaction_topology_similarity(
    query: Mapping[str, Any], precedent: Mapping[str, Any]
) -> float:
    query_topology = query.get("topology") or {}
    precedent_topology = precedent.get("topology") or {}
    if not isinstance(query_topology, Mapping) or not isinstance(
        precedent_topology, Mapping
    ):
        return 0.0
    query_scope = reaction_scope(query)
    precedent_scope = reaction_scope(precedent)
    if not query_scope or not precedent_scope or query_scope != precedent_scope:
        return 0.0

    score = 0.6
    query_bond_scopes = tuple(query_topology.get("formed_bond_scopes") or ())
    precedent_bond_scopes = tuple(precedent_topology.get("formed_bond_scopes") or ())
    if query_bond_scopes or precedent_bond_scopes:
        score += 0.15 * _jaccard(query_bond_scopes, precedent_bond_scopes)

    query_delta = query_topology.get("ring_count_delta")
    precedent_delta = precedent_topology.get("ring_count_delta")
    if query_delta is not None and precedent_delta is not None:
        score += 0.1 * float(int(query_delta) == int(precedent_delta))

    query_rings = tuple(
        int(value) for value in query_topology.get("formed_ring_sizes") or ()
    )
    precedent_rings = tuple(
        int(value) for value in precedent_topology.get("formed_ring_sizes") or ()
    )
    if not query_rings and not precedent_rings:
        score += 0.15
    elif query_rings and precedent_rings:
        distance = abs(min(query_rings) - min(precedent_rings))
        score += 0.15 * max(0.0, 1.0 - distance / 4.0)
    return round(min(score, 1.0), 6)


def generic_signature_similarity(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
) -> Tuple[float, Dict[str, float]]:
    """Return an interpretable score; absent features never count as matches."""
    edit_fields = (
        "formed_bond_types",
        "broken_bond_types",
        "order_changes",
        "hydrogen_changes",
    )
    query_edits = tuple(
        f"{field}:{value}" for field in edit_fields for value in query.get(field) or ()
    )
    precedent_edits = tuple(
        f"{field}:{value}"
        for field in edit_fields
        for value in precedent.get(field) or ()
    )
    query_transformation = str(query.get("transformation_class") or "")
    precedent_transformation = str(precedent.get("transformation_class") or "")
    query_family = str(query.get("named_family") or "")
    precedent_family = str(precedent.get("named_family") or "")
    components = {
        "edit_topology": _jaccard(query_edits, precedent_edits),
        "reaction_events": _multiset_jaccard(
            _event_tokens(query), _event_tokens(precedent)
        ),
        "handles": _jaccard(
            _partner_tokens(query, "handle_tokens"),
            _partner_tokens(precedent, "handle_tokens"),
        ),
        "contexts": _jaccard(
            _partner_tokens(query, "anchor_contexts"),
            _partner_tokens(precedent, "anchor_contexts"),
        ),
        "environment": _jaccard(
            _environment_tokens(query), _environment_tokens(precedent)
        ),
        "spectators": _jaccard(_spectator_tokens(query), _spectator_tokens(precedent)),
        "reaction_topology": _reaction_topology_similarity(query, precedent),
        "transformation": float(
            bool(query_transformation)
            and query_transformation == precedent_transformation
        ),
        "family": float(bool(query_family) and query_family == precedent_family),
    }
    weights = load_generic_retrieval_rules()["similarity_weights"]
    score = sum(float(weights[name]) * components[name] for name in weights)
    return round(score, 6), components


__all__ = [
    "generic_signature_similarity",
    "load_generic_retrieval_rules",
    "reaction_scope",
    "retrieve_compatible_generic_pool",
    "retrieve_generic_pool",
]
