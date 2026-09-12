"""Interpretable, versioned similarity for generic reaction signatures."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Tuple

from reactive_taxonomy.descriptors.tokens import reactivity_profile_is_usable

from .signature_features import (
    environment_profile_similarity,
    substituent_profile_similarity,
)

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_similarity.v1.json"
_DEFINITION_ID = "generic_similarity.v1"
_SCHEMA_VERSION = "1.0"
_FEATURES = (
    "edit_topology",
    "reaction_events",
    "reaction_topology",
    "handles",
    "contexts",
    "environment",
    "spectators",
    "transformation",
    "family",
)


@dataclass(frozen=True)
class SimilarityAssessment:
    """Total similarity with component values and weighted contributions."""

    score: float
    components: Dict[str, float]
    contributions: Dict[str, float]
    definition_id: str
    definition_version: str


@lru_cache(maxsize=1)
def load_generic_similarity_rules() -> Dict[str, Any]:
    """Load and validate the declarative generic-similarity definition."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    rules["weights"] = _validated_similarity_weights(rules)
    _validate_environment_weights(rules)
    return rules


def _validate_environment_weights(rules: Mapping[str, Any]) -> None:
    environment_weights = rules.get("environment_component_weights") or {}
    if (
        not isinstance(environment_weights, Mapping)
        or set(environment_weights) != {"molecular", "substituent"}
        or any(not 0 <= float(value) <= 1 for value in environment_weights.values())
        or abs(sum(float(value) for value in environment_weights.values()) - 1) > 1e-9
    ):
        raise ValueError("environment component weights must sum to one")


def _validated_similarity_weights(
    rules: Mapping[str, Any],
) -> Dict[str, float]:
    if str(rules.get("schema_version") or "") != _SCHEMA_VERSION:
        raise ValueError("unsupported generic similarity definition schema")
    if str(rules.get("definition_id") or "") != _DEFINITION_ID:
        raise ValueError("unexpected generic similarity definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("generic similarity definition requires calibration status")
    weights = rules.get("weights")
    if not isinstance(weights, Mapping) or set(weights) != set(_FEATURES):
        raise ValueError("generic similarity weights do not match feature vocabulary")
    normalized = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in normalized.values()):
        raise ValueError("generic similarity weights must be non-negative")
    if abs(sum(normalized.values()) - 1.0) > 1e-9:
        raise ValueError("generic similarity weights must sum to one")
    return normalized


def validate_generic_similarity_rules(rules: Mapping[str, Any]) -> None:
    """Validate a generic-similarity definition without loading global state."""
    _validated_similarity_weights(rules)
    _validate_environment_weights(rules)


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    """Set Jaccard where absent features never count as matches."""
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
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
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
        score += 0.15 * jaccard(query_bond_scopes, precedent_bond_scopes)

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


def assess_signature_similarity(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
    *,
    query_reaction_core: Mapping[str, Any] | None = None,
    precedent_reaction_core: Mapping[str, Any] | None = None,
    query_molecular_features: Mapping[str, Any] | None = None,
    precedent_molecular_features: Mapping[str, Any] | None = None,
) -> SimilarityAssessment:
    """Calculate a complete, auditable structural similarity assessment."""
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
    profile_similarity = substituent_profile_similarity(
        query_reaction_core,
        precedent_reaction_core,
    )
    query_features = (
        query if query_molecular_features is None else query_molecular_features
    )
    precedent_features = (
        precedent
        if precedent_molecular_features is None
        else precedent_molecular_features
    )
    molecular_similarity = environment_profile_similarity(
        query_features, precedent_features
    )
    available_profiles = all(
        any(
            reactivity_profile_is_usable(partner.get("reactivity_profile"))
            for partner in features.get("partners") or ()
            if isinstance(partner, Mapping)
        )
        for features in (query_features, precedent_features)
    )
    if profile_similarity is None:
        environment_similarity = molecular_similarity
    elif available_profiles:
        weights = load_generic_similarity_rules()["environment_component_weights"]
        environment_similarity = (
            float(weights["substituent"]) * profile_similarity
            + float(weights["molecular"]) * molecular_similarity
        )
    else:
        environment_similarity = profile_similarity
    components = {
        "edit_topology": jaccard(query_edits, precedent_edits),
        "reaction_events": _multiset_jaccard(
            _event_tokens(query), _event_tokens(precedent)
        ),
        "handles": jaccard(
            _partner_tokens(query, "handle_tokens"),
            _partner_tokens(precedent, "handle_tokens"),
        ),
        "contexts": jaccard(
            _partner_tokens(query, "anchor_contexts"),
            _partner_tokens(precedent, "anchor_contexts"),
        ),
        "environment": environment_similarity,
        "spectators": jaccard(_spectator_tokens(query), _spectator_tokens(precedent)),
        "reaction_topology": _reaction_topology_similarity(query, precedent),
        "transformation": float(
            bool(query_transformation)
            and query_transformation == precedent_transformation
        ),
        "family": float(bool(query_family) and query_family == precedent_family),
    }
    rules = load_generic_similarity_rules()
    contributions = {
        name: float(rules["weights"][name]) * components[name]
        for name in _FEATURES
    }
    return SimilarityAssessment(
        score=round(sum(contributions.values()), 6),
        components={name: round(components[name], 6) for name in _FEATURES},
        contributions={name: round(contributions[name], 6) for name in _FEATURES},
        definition_id=str(rules["definition_id"]),
        definition_version=str(
            rules.get("definition_version", rules["schema_version"])
        ),
    )


def generic_signature_similarity(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
) -> Tuple[float, Dict[str, float]]:
    """Compatibility wrapper returning the historical score and components."""
    assessment = assess_signature_similarity(query, precedent)
    return assessment.score, assessment.components


__all__ = [
    "SimilarityAssessment",
    "assess_signature_similarity",
    "generic_signature_similarity",
    "jaccard",
    "load_generic_similarity_rules",
    "reaction_scope",
    "validate_generic_similarity_rules",
]
