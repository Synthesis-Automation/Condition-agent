"""Interpretable similarity for unresolved-reaction fallback descriptors."""

from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Tuple

from .similarity import SimilarityAssessment


_RULES_PATH = Path(__file__).with_name("definitions") / "fallback_retrieval.v1.json"
_FEATURE_FIELDS = {
    "reactant_components": ("reactant_component_tokens",),
    "product_components": ("product_component_tokens",),
    "reactant_sites": ("reactant_site_tokens",),
    "product_sites": ("product_site_tokens",),
    "functional_groups": (
        "reactant_group_tokens",
        "product_group_tokens",
    ),
    "contexts": ("context_tokens",),
    "candidate_grammars": ("candidate_grammar_tokens",),
    "candidate_transformations": ("candidate_transformation_tokens",),
    "candidate_handles": ("candidate_handle_tokens",),
    "candidate_edits": ("candidate_edit_tokens",),
    "observed_or_partial_edits": ("verified_edit_tokens",),
    "bond_inventory_delta": ("bond_inventory_delta_tokens",),
    "element_delta": ("element_delta_tokens",),
}
_HIGH_SIGNAL_FIELDS = (
    "candidate_grammar_tokens",
    "candidate_transformation_tokens",
    "candidate_edit_tokens",
    "verified_edit_tokens",
    "bond_inventory_delta_tokens",
    "reactant_site_tokens",
    "product_group_tokens",
)


@lru_cache(maxsize=1)
def load_fallback_retrieval_rules() -> dict[str, Any]:
    """Load and validate the conservative fallback policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.2":
        raise ValueError("unsupported fallback retrieval definition schema")
    if str(rules.get("definition_id") or "") != "fallback_retrieval.v1":
        raise ValueError("unexpected fallback retrieval definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("fallback retrieval policy requires calibration status")
    weights = rules.get("weights")
    if not isinstance(weights, Mapping) or set(weights) != set(_FEATURE_FIELDS):
        raise ValueError("fallback weights do not match feature vocabulary")
    normalized = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in normalized.values()):
        raise ValueError("fallback weights must be non-negative")
    if abs(sum(normalized.values()) - 1.0) > 1e-9:
        raise ValueError("fallback weights must sum to one")
    rules["weights"] = normalized
    minimum = float(rules["minimum_similarity"])
    limited = float(rules["limited_support_minimum_similarity"])
    if not 0.0 < minimum <= limited <= 1.0:
        raise ValueError("invalid fallback similarity thresholds")
    if int(rules["minimum_independent_support"]) < 1:
        raise ValueError("fallback independent support must be positive")
    if int(rules["minimum_shared_high_signal_features"]) < 1:
        raise ValueError("fallback high-signal requirement must be positive")
    if int(rules["candidate_limit"]) < 1:
        raise ValueError("fallback candidate limit must be positive")
    source_requirements = rules.get("condition_source_requirements")
    if not isinstance(source_requirements, Mapping):
        raise ValueError("fallback condition source requirements must be an object")
    for requirement_id, requirement in source_requirements.items():
        if not str(requirement_id) or not isinstance(requirement, Mapping):
            raise ValueError("invalid fallback condition source requirement")
        if not tuple(requirement.get("elements") or ()):
            raise ValueError("condition source requirement needs elements")
        if not (
            tuple(requirement.get("family_ids") or ())
            or tuple(requirement.get("substance_ids") or ())
        ):
            raise ValueError("condition source requirement needs registry identities")
    return rules


def _tokens(
    descriptor: Mapping[str, Any],
    fields: Iterable[str],
) -> Tuple[str, ...]:
    return tuple(
        f"{field}:{value}" for field in fields for value in descriptor.get(field) or ()
    )


def _multiset_jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_counts = Counter(left)
    right_counts = Counter(right)
    if not left_counts or not right_counts:
        return 0.0
    keys = set(left_counts).union(right_counts)
    intersection = sum(min(left_counts[key], right_counts[key]) for key in keys)
    union = sum(max(left_counts[key], right_counts[key]) for key in keys)
    return intersection / union if union else 0.0


def fallback_index_tokens(descriptor: Mapping[str, Any]) -> Tuple[str, ...]:
    """Return high-signal tokens used only to narrow fallback comparisons."""
    if not bool(descriptor.get("retrieval_eligible")):
        return ()
    return tuple(sorted(set(_tokens(descriptor, _HIGH_SIGNAL_FIELDS))))


def shared_high_signal_count(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
) -> int:
    """Count distinct shared feature groups, not raw correlated tokens."""
    return sum(
        bool(
            set(str(value) for value in query.get(field) or ())
            & set(str(value) for value in precedent.get(field) or ())
        )
        for field in _HIGH_SIGNAL_FIELDS
    )


def assess_fallback_similarity(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
) -> SimilarityAssessment:
    """Score two descriptors while exposing every chemistry feature group."""
    rules = load_fallback_retrieval_rules()
    feature_tokens = {
        name: (
            _tokens(query, fields),
            _tokens(precedent, fields),
        )
        for name, fields in _FEATURE_FIELDS.items()
    }
    components = {
        name: _multiset_jaccard(left, right)
        for name, (left, right) in feature_tokens.items()
    }
    available = {
        name for name, (left, right) in feature_tokens.items() if left or right
    }
    denominator = sum(float(rules["weights"][name]) for name in available)
    contributions = {
        name: (
            float(rules["weights"][name]) / denominator * components[name]
            if denominator and name in available
            else 0.0
        )
        for name in _FEATURE_FIELDS
    }
    return SimilarityAssessment(
        score=round(sum(contributions.values()), 6),
        components={name: round(value, 6) for name, value in components.items()},
        contributions={name: round(value, 6) for name, value in contributions.items()},
        definition_id=str(rules["definition_id"]),
        definition_version=str(rules["schema_version"]),
    )


def compatibility_signature_from_fallback(
    descriptor: Mapping[str, Any],
) -> dict[str, Any]:
    """Expose all observed group tags conservatively to recipe compatibility."""
    tags = tuple(
        sorted(set(str(tag) for tag in descriptor.get("compatibility_tags") or ()))
    )
    return {
        "spectator_groups": (
            {
                "group_id": "fallback_all_observed_groups",
                "tags": tags,
            },
        )
        if tags
        else (),
        "events": (),
        "named_family": None,
        "family_confidence": 0.0,
    }


__all__ = [
    "assess_fallback_similarity",
    "compatibility_signature_from_fallback",
    "fallback_index_tokens",
    "load_fallback_retrieval_rules",
    "shared_high_signal_count",
]
