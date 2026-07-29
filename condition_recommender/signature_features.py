"""Shared interpretable feature extraction from serialized reaction signatures."""

from __future__ import annotations

from typing import Any, Mapping, Tuple

from reactive_taxonomy.descriptors import reactivity_profile_tokens


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


__all__ = ["environment_profile_similarity", "environment_tokens"]
