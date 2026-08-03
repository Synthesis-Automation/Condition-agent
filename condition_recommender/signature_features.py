"""Shared interpretable feature extraction from serialized reaction signatures."""

from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Tuple

from reactive_taxonomy.descriptors import reactivity_profile_tokens


_SUBSTITUENT_RULES_PATH = (
    Path(__file__).with_name("definitions")
    / "substituent_profile_similarity.v1.json"
)


@lru_cache(maxsize=1)
def load_substituent_profile_similarity_rules() -> Mapping[str, Any]:
    """Load and validate hierarchical R-group similarity weights."""
    with _SUBSTITUENT_RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if rules.get("schema_version") != "1.0" or rules.get(
        "definition_id"
    ) != "substituent_profile_similarity.v1":
        raise ValueError("unsupported substituent-profile similarity definition")
    weights = rules.get("level_weights") or {}
    if set(weights) != {"L0", "L1", "L2", "L3"}:
        raise ValueError("substituent-profile similarity levels are incomplete")
    normalized = {str(key): float(value) for key, value in weights.items()}
    if any(value < 0.0 for value in normalized.values()) or abs(
        sum(normalized.values()) - 1.0
    ) > 1e-9:
        raise ValueError("substituent-profile similarity weights must sum to one")
    rules["level_weights"] = normalized
    return rules


def environment_tokens(signature: Mapping[str, Any]) -> Tuple[str, ...]:
    """Return stable context-aware profile and nearby-group tokens."""
    tokens = []
    for partner in signature.get("partners") or ():
        role = str(partner.get("role") or "unassigned")
        for value in reactivity_profile_tokens(
            partner.get("reactivity_profile")
        ):
            tokens.append(f"{role}:{value}")
        for group in partner.get("nearby_groups") or ():
            group_id = group.get("group_id")
            if group_id:
                distance = group.get("distance")
                distance_token = (
                    f":d{int(distance)}" if distance is not None else ""
                )
                tokens.append(
                    f"{role}:nearby:{group_id}{distance_token}"
                )
    return tuple(sorted(tokens))


def _set_similarity(left: Tuple[str, ...], right: Tuple[str, ...]) -> float:
    a, b = set(left), set(right)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _multiset_similarity(left: Tuple[str, ...], right: Tuple[str, ...]) -> float:
    left_counts = Counter(left)
    right_counts = Counter(right)
    if not left_counts or not right_counts:
        return 0.0
    keys = set(left_counts).union(right_counts)
    intersection = sum(min(left_counts[key], right_counts[key]) for key in keys)
    union = sum(max(left_counts[key], right_counts[key]) for key in keys)
    return intersection / union if union else 0.0


def substituent_profile_tokens(
    reaction_core: Mapping[str, Any] | None,
    *,
    level: str | None = None,
) -> Tuple[str, ...]:
    """Return side- and continuity-aware profile tokens from a serialized core."""
    if not reaction_core:
        return ()
    tokens = []
    for remote in reaction_core.get("remote_subgraphs") or ():
        if not isinstance(remote, Mapping):
            continue
        side = str(remote.get("side") or "unknown")
        continuity = str(remote.get("continuity") or "unresolved")
        for port in remote.get("attachment_ports") or ():
            if not isinstance(port, Mapping):
                continue
            profile = port.get("substituent_profile") or {}
            if not isinstance(profile, Mapping):
                continue
            for value in profile.get("feature_tokens") or ():
                token = str(value)
                token_level = token.split(":", 1)[0]
                if level is not None and token_level != level:
                    continue
                tokens.append(f"{side}:{continuity}:{token}")
    return tuple(sorted(tokens))


def substituent_profile_similarity(
    query_core: Mapping[str, Any] | None,
    precedent_core: Mapping[str, Any] | None,
) -> float | None:
    """Compare port profiles hierarchically, or return ``None`` when absent."""
    if not substituent_profile_tokens(query_core) or not substituent_profile_tokens(
        precedent_core
    ):
        return None
    rules = load_substituent_profile_similarity_rules()
    return sum(
        float(weight)
        * _multiset_similarity(
            substituent_profile_tokens(query_core, level=level),
            substituent_profile_tokens(precedent_core, level=level),
        )
        for level, weight in rules["level_weights"].items()
    )


def _numeric_similarity(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
) -> float:
    left_steric = left.get("steric") or {}
    right_steric = right.get("steric") or {}
    steric = 0.0
    if isinstance(left_steric, Mapping) and isinstance(
        right_steric, Mapping
    ):
        left_score = left_steric.get("accessibility_score")
        right_score = right_steric.get("accessibility_score")
        if isinstance(left_score, (int, float)) and isinstance(
            right_score, (int, float)
        ):
            steric = max(0.0, 1.0 - abs(float(left_score) - float(right_score)))
    left_electronic = left.get("electronic") or {}
    right_electronic = right.get("electronic") or {}
    electronic = 0.0
    if isinstance(left_electronic, Mapping) and isinstance(
        right_electronic, Mapping
    ):
        left_axis = str(left_electronic.get("activation_axis") or "")
        right_axis = str(right_electronic.get("activation_axis") or "")
        left_score = left_electronic.get("activation_score")
        right_score = right_electronic.get("activation_score")
        if (
            left_axis
            and left_axis == right_axis
            and isinstance(left_score, (int, float))
            and isinstance(right_score, (int, float))
        ):
            electronic = max(
                0.0,
                1.0 - abs(float(left_score) - float(right_score)) / 2.0,
            )
    return 0.6 * steric + 0.4 * electronic


def environment_profile_similarity(
    query: Mapping[str, Any],
    precedent: Mapping[str, Any],
) -> float:
    """Compare role-aligned categorical profiles and compatible numeric axes."""
    query_partners = tuple(
        partner
        for partner in query.get("partners") or ()
        if isinstance(partner, Mapping)
        and isinstance(partner.get("reactivity_profile"), Mapping)
    )
    precedent_partners = tuple(
        partner
        for partner in precedent.get("partners") or ()
        if isinstance(partner, Mapping)
        and isinstance(partner.get("reactivity_profile"), Mapping)
    )
    if not query_partners or not precedent_partners:
        return 0.0
    used: set[int] = set()
    scores = []
    for query_partner in query_partners:
        query_profile = query_partner["reactivity_profile"]
        query_role = str(query_partner.get("role") or "")
        query_kind = str(query_profile.get("context_kind") or "")
        candidates = [
            (index, partner)
            for index, partner in enumerate(precedent_partners)
            if index not in used
            and (
                (
                    query_role
                    and str(partner.get("role") or "") == query_role
                )
                or (
                    not query_role
                    and str(
                        partner["reactivity_profile"].get("context_kind") or ""
                    )
                    == query_kind
                )
            )
        ]
        if not candidates:
            scores.append(0.0)
            continue
        query_tokens = reactivity_profile_tokens(query_profile)
        ranked = []
        for index, partner in candidates:
            profile = partner["reactivity_profile"]
            categorical = _set_similarity(
                query_tokens, reactivity_profile_tokens(profile)
            )
            numeric = _numeric_similarity(query_profile, profile)
            ranked.append((0.75 * categorical + 0.25 * numeric, index))
        score, selected = max(ranked, key=lambda item: (item[0], -item[1]))
        used.add(selected)
        scores.append(score)
    denominator = max(len(query_partners), len(precedent_partners))
    return sum(scores) / denominator if denominator else 0.0


__all__ = [
    "environment_profile_similarity",
    "environment_tokens",
    "load_substituent_profile_similarity_rules",
    "substituent_profile_similarity",
    "substituent_profile_tokens",
]
