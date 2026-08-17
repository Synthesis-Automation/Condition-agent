"""Weak-label condition fallback and diverse screening-array generation."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple

from condition_registry import CONDITION_RECIPE_COMPONENT_BUCKETS
from reactive_taxonomy import (
    featurize_reaction,
    load_reaction_type_hint_definitions,
)

from .compatibility import assess_recipe_compatibility
from .models import (
    WeakLabelConditionRecommendation,
    WeakLabelRecommendationResult,
    WeakLabelSourceMatch,
)
from .weak_label_indexing import (
    WeakLabelIndexedObservation,
    WeakLabelParticipant,
    load_weak_label_index,
)


_RULES_PATH = (
    Path(__file__).with_name("definitions") / "weak_label_retrieval.v1.json"
)
DEFAULT_WEAK_LABEL_RECORDS_PATH = (
    Path(__file__).resolve().parents[1]
    / "raw_dataset"
    / "weak_label"
    / "v2.1_cleaned.csv"
)


@lru_cache(maxsize=1)
def load_weak_label_retrieval_rules() -> Dict[str, Any]:
    """Load and validate weak-label retrieval and screening policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if rules.get("schema_version") != "1.0" or rules.get(
        "definition_id"
    ) != "weak_label_retrieval.v1":
        raise ValueError("unsupported weak-label retrieval definition")
    known_type_ids = {
        str(value["type_id"])
        for value in load_reaction_type_hint_definitions().values()
    }
    configured = set((rules.get("reaction_type_sources") or {}).keys())
    if not configured or not configured <= known_type_ids:
        raise ValueError(
            "weak-label retrieval contains unknown reaction type hints: "
            f"{sorted(configured - known_type_ids)}"
        )
    if set(rules.get("supported_reaction_scopes") or ()) - {
        "intramolecular",
        "intermolecular",
        "mixed",
        "unimolecular",
        "unresolved",
    }:
        raise ValueError("invalid weak-label supported reaction scope")
    for key in ("participant_similarity_weights", "ranking_weights"):
        weights = {name: float(value) for name, value in (rules.get(key) or {}).items()}
        if not weights or abs(sum(weights.values()) - 1.0) > 1e-9:
            raise ValueError(f"weak-label weights must sum to one: {key}")
        rules[key] = weights
    diversity_weight = float(rules.get("screening_diversity_weight") or 0.0)
    if not 0.0 <= diversity_weight <= 1.0:
        raise ValueError("screening diversity weight must be in [0, 1]")
    if rules.get("require_complete_participant_signatures") is not True:
        raise ValueError("weak-label retrieval requires complete participant labels")
    return rules


def validate_weak_label_retrieval_rules() -> list[str]:
    """Return retrieval-definition validation errors."""
    try:
        load_weak_label_retrieval_rules()
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return [str(exc)]
    return []


def _normalize_signature(value: str) -> str:
    return str(value or "").strip().replace("HeteroAr", "HetAr")


def _signature_similarity(query: str, precedent: str) -> float:
    query = _normalize_signature(query)
    precedent = _normalize_signature(precedent)
    if not query or not precedent:
        return 0.0
    if query == precedent:
        return 1.0
    query_parts = query.split("|")
    precedent_parts = precedent.split("|")
    if query_parts[0] != precedent_parts[0]:
        return 0.0
    if (
        query_parts[0] == "XH"
        and len(query_parts) > 1
        and len(precedent_parts) > 1
        and query_parts[1] != precedent_parts[1]
    ):
        return 0.0
    comparable = max(len(query_parts), len(precedent_parts)) - 1
    if comparable <= 0:
        return 0.35
    matches = sum(
        left == right
        for left, right in zip(query_parts[1:], precedent_parts[1:])
    )
    return round(0.35 + 0.65 * matches / comparable, 6)


def _known_incompatible(query: str, precedent: str) -> bool:
    query = _normalize_signature(query)
    precedent = _normalize_signature(precedent)
    if not query or not precedent:
        return False
    if _signature_similarity(query, precedent) == 0.0:
        return True
    query_parts = query.split("|")
    precedent_parts = precedent.split("|")
    if query_parts[0] != "TM" or precedent_parts[0] != "TM":
        return False
    families = load_weak_label_retrieval_rules()[
        "transfer_handle_token_families"
    ]
    return families.get(query_parts[-1], query_parts[-1]) != families.get(
        precedent_parts[-1], precedent_parts[-1]
    )


def _qualifier_similarity(query: Any, precedent: WeakLabelParticipant) -> float:
    comparisons = []
    if precedent.center_class:
        comparisons.append(
            query.center_substitution_class == precedent.center_class
        )
    if precedent.attachment_class:
        comparisons.append(
            precedent.attachment_class in query.attachment_carbon_classes
        )
    if precedent.alpha_branched is not None:
        comparisons.append(query.alpha_branched is precedent.alpha_branched)
    if not comparisons:
        return 0.5
    signed = sum(1.0 if value else -1.0 for value in comparisons)
    return max(0.0, min(1.0, 0.5 + 0.5 * signed / 3.0))


def _pair_similarity(
    query: Sequence[Any],
    precedent: Tuple[WeakLabelParticipant, WeakLabelParticipant],
) -> tuple[float, float, float, tuple[int, int]]:
    weights = load_weak_label_retrieval_rules()["participant_similarity_weights"]
    alternatives = []
    for order in ((0, 1), (1, 0)):
        ordered = (precedent[order[0]], precedent[order[1]])
        if any(not item.signature for item in ordered):
            alternatives.append((0.0, 0.0, 0.0, order))
            continue
        if any(
            _known_incompatible(query[index].canonical_signature, ordered[index].signature)
            for index in range(2)
        ):
            alternatives.append((0.0, 0.0, 0.0, order))
            continue
        signatures = tuple(
            _signature_similarity(
                query[index].canonical_signature,
                ordered[index].signature,
            )
            for index in range(2)
        )
        qualifiers = tuple(
            _qualifier_similarity(query[index], ordered[index])
            for index in range(2)
        )
        signature_score = sum(signatures) / 2.0
        qualifier_score = sum(qualifiers) / 2.0
        alternatives.append(
            (
                float(weights["signature"]) * signature_score
                + float(weights["qualifier"]) * qualifier_score,
                signature_score,
                qualifier_score,
                order,
            )
        )
    return max(alternatives, key=lambda value: value[:3])


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(-20.0, min(20.0, value))))


def _mean(values: Iterable[float]) -> float:
    items = tuple(values)
    return sum(items) / len(items) if items else 0.0


def _recipe_tokens(recipe: Mapping[str, Any]) -> frozenset[str]:
    tokens = set()
    for bucket in CONDITION_RECIPE_COMPONENT_BUCKETS:
        for component in recipe.get(bucket) or ():
            identity = (
                component.get("substance_id")
                or component.get("canonical_name")
                or component.get("raw_identifier")
                or "unresolved"
            )
            tokens.add(f"{bucket}:{identity}")
    temperature = recipe.get("temperature_c")
    if isinstance(temperature, (int, float)):
        tokens.add(f"temperature_bin:{int(float(temperature) // 25)}")
    return frozenset(tokens)


def _jaccard_distance(left: frozenset[str], right: frozenset[str]) -> float:
    if not left and not right:
        return 0.0
    return 1.0 - len(left & right) / len(left | right)


def _diverse_selection(
    recommendations: Sequence[WeakLabelConditionRecommendation],
    *,
    size: int,
) -> tuple[WeakLabelConditionRecommendation, ...]:
    rules = load_weak_label_retrieval_rules()
    limit = int(rules.get("screening_candidate_limit") or 250)
    remaining = list(recommendations[:limit])
    selected: list[WeakLabelConditionRecommendation] = []
    diversity_weight = float(rules["screening_diversity_weight"])
    tokens = {item.recipe_id: _recipe_tokens(item.resolved_recipe) for item in remaining}
    while remaining and len(selected) < size:
        if not selected:
            choice = remaining[0]
        else:
            choice = max(
                remaining,
                key=lambda item: (
                    (1.0 - diversity_weight) * item.score
                    + diversity_weight
                    * min(
                        _jaccard_distance(
                            tokens[item.recipe_id], tokens[value.recipe_id]
                        )
                        for value in selected
                    ),
                    item.score,
                    item.expected_yield_pct
                    if item.expected_yield_pct is not None
                    else -1.0,
                    item.recipe_id,
                ),
            )
        selected.append(choice)
        remaining.remove(choice)
    return tuple(
        replace(
            item,
            rank=index,
            explanation=tuple(
                dict.fromkeys(
                    (*item.explanation, "Selected for condition-space diversity")
                )
            ),
        )
        for index, item in enumerate(selected, start=1)
    )


def _query_compatibility_payload(analysis: Any) -> Dict[str, Any]:
    payload = (
        asdict(analysis.reaction_signature)
        if analysis.reaction_signature is not None
        else {}
    )
    payload.update(
        {
            "named_family": analysis.named_family,
            "family_confidence": 1.0 if analysis.named_family else 0.0,
            "spectator_groups": tuple(
                asdict(group) for group in analysis.spectator_groups
            ),
        }
    )
    return payload


def _rank_recipes(
    query_participants: Sequence[Any],
    rows: Sequence[WeakLabelIndexedObservation],
    *,
    query_compatibility: Mapping[str, Any],
) -> tuple[
    tuple[WeakLabelConditionRecommendation, ...], int, int
]:
    rules = load_weak_label_retrieval_rules()
    assessed = []
    excluded = 0
    compatibility_cache = {}
    for row in rows:
        label_score, signature_score, qualifier_score, participant_order = _pair_similarity(
            query_participants, row.participants
        )
        if signature_score <= 0.0:
            excluded += 1
            continue
        compatibility = compatibility_cache.get(row.recipe_id)
        if compatibility is None:
            compatibility = assess_recipe_compatibility(
                query_compatibility,
                row.resolved_recipe,
            )
            compatibility_cache[row.recipe_id] = compatibility
        if not compatibility.compatible:
            excluded += 1
            continue
        assessed.append(
            (
                label_score,
                signature_score,
                qualifier_score,
                participant_order,
                compatibility,
                row,
            )
        )

    groups = defaultdict(list)
    for value in assessed:
        groups[value[-1].recipe_id].append(value)
    weights = rules["ranking_weights"]
    maximum = int(rules.get("maximum_precedents_per_recipe") or 5)
    ranked = []
    for recipe_id, members in groups.items():
        members.sort(
            key=lambda item: (
                -item[0],
                -item[-1].yield_pct,
                item[-1].source_row_number,
            )
        )
        similarity_weights = tuple(max(item[0], 0.05) for item in members)
        weight_sum = sum(similarity_weights)
        expected_yield = sum(
            weight * item[-1].yield_pct
            for weight, item in zip(similarity_weights, members)
        ) / weight_sum
        z_values = tuple(
            (weight, item[-1].z_score)
            for weight, item in zip(similarity_weights, members)
            if item[-1].z_score is not None
        )
        mean_z = (
            sum(weight * float(value) for weight, value in z_values)
            / sum(weight for weight, _ in z_values)
            if z_values
            else None
        )
        (
            best_label,
            best_signature,
            best_qualifier,
            _best_order,
            compatibility,
            best_row,
        ) = members[0]
        support_score = min(1.0, math.log1p(len(members)) / math.log1p(10))
        score = (
            float(weights["label_similarity"]) * best_label
            + float(weights["compatibility"]) * compatibility.score
            + float(weights["yield"]) * expected_yield / 100.0
            + float(weights["z_score"]) * _sigmoid(mean_z or 0.0)
            + float(weights["support"]) * support_score
        )
        ranked.append(
            WeakLabelConditionRecommendation(
                rank=0,
                recipe_id=recipe_id,
                resolved_recipe=dict(best_row.resolved_recipe),
                score=round(score, 6),
                label_similarity=round(best_label, 6),
                signature_similarity=round(best_signature, 6),
                qualifier_similarity=round(best_qualifier, 6),
                compatibility_score=compatibility.score,
                expected_yield_pct=round(expected_yield, 2),
                mean_z_score=round(mean_z, 3) if mean_z is not None else None,
                support=len(members),
                source_reaction_types=tuple(
                    sorted({item[-1].reaction_type for item in members})
                ),
                source_row_numbers=tuple(
                    item[-1].source_row_number for item in members[:maximum]
                ),
                source_matches=tuple(
                    WeakLabelSourceMatch(
                        source_row_number=item[-1].source_row_number,
                        source_reaction_type=item[-1].reaction_type,
                        participant_roles=tuple(
                            str(query_participants[index].role or f"participant_{index + 1}")
                            for index in range(2)
                        ),
                        participant_display_labels=tuple(
                            item[-1].participants[source_index].display_label
                            for source_index in item[3]
                        ),
                        participant_signatures=tuple(
                            item[-1].participants[source_index].signature
                            for source_index in item[3]
                        ),
                        label_similarity=round(item[0], 6),
                        signature_similarity=round(item[1], 6),
                        qualifier_similarity=round(item[2], 6),
                    )
                    for item in members[:maximum]
                ),
                explanation=(
                    "Compatible graph-derived reaction-type hint",
                    (
                        "Exact unordered participant-signature pair"
                        if best_signature == 1.0
                        else "Compatible unordered participant-signature pair"
                    ),
                    f"Recipe compatibility score {compatibility.score:.2f}",
                    f"Supported by {len(members)} weak-label observation(s)",
                ),
                compatibility_evidence=compatibility.evidence,
                cautions=(
                    "Source reactions are not structure-verified",
                    "Reaction-type labels are weak evidence",
                ),
            )
        )
    ranked.sort(
        key=lambda item: (
            -item.score,
            -(item.expected_yield_pct or -1.0),
            item.recipe_id,
        )
    )
    return (
        tuple(replace(item, rank=index) for index, item in enumerate(ranked, start=1)),
        len(assessed),
        excluded,
    )


def recommend_weak_label_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = DEFAULT_WEAK_LABEL_RECORDS_PATH,
    top_k: int = 5,
    mode: str = "fallback",
    source_reaction_type_hint: str | None = None,
) -> WeakLabelRecommendationResult:
    """Recommend or diversify recipes from explicitly unverified label rows."""
    normalized_mode = str(mode).strip().casefold()
    result_mode = (
        "weak_label_screening"
        if normalized_mode == "screening"
        else "weak_label_fallback"
    )
    if normalized_mode not in {"fallback", "screening"}:
        return WeakLabelRecommendationResult(
            reaction_smiles,
            False,
            "weak_label_fallback",
            error="UNSUPPORTED_WEAK_LABEL_MODE",
        )
    if top_k < 1:
        return WeakLabelRecommendationResult(
            reaction_smiles,
            False,
            result_mode,
            error="TOP_K_MUST_BE_POSITIVE",
        )
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid or analysis.reaction_signature is None:
        return WeakLabelRecommendationResult(
            reaction_smiles,
            False,
            result_mode,
            warnings=("WEAK_LABEL_FALLBACK_REQUIRES_VERIFIED_GRAPH_EDITS",),
            error=analysis.error or "QUERY_HAS_NO_VERIFIED_REACTION_SIGNATURE",
        )
    interpretation = analysis.interpretation
    hints = tuple(interpretation.reaction_type_hints if interpretation else ())
    rules = load_weak_label_retrieval_rules()
    supported_scopes = set(rules["supported_reaction_scopes"])
    candidates = tuple(
        hint
        for hint in hints
        if hint.type_id in rules["reaction_type_sources"]
        and hint.reaction_scope in supported_scopes
        and len(hint.participants) == 2
        and not hint.warnings
    )
    if not candidates:
        return WeakLabelRecommendationResult(
            reaction_smiles,
            False,
            result_mode,
            warnings=("WEAK_LABEL_QUERY_TYPE_OR_PARTICIPANTS_UNSUPPORTED",),
            error="QUERY_NOT_SUPPORTED_BY_WEAK_LABEL_DATASET",
        )
    primary_id = (
        interpretation.primary_reaction_type_hint_id if interpretation else None
    )
    hint = next(
        (value for value in candidates if value.hint_id == primary_id),
        candidates[0],
    )
    source_types = tuple(rules["reaction_type_sources"][hint.type_id])
    explicit = str(source_reaction_type_hint or "").strip()
    if explicit:
        if explicit not in source_types:
            return WeakLabelRecommendationResult(
                reaction_smiles,
                False,
                result_mode,
                reaction_type_hint_id=hint.hint_id,
                reaction_type_id=hint.type_id,
                source_reaction_type_candidates=source_types,
                query_participants=tuple(asdict(value) for value in hint.participants),
                warnings=("SOURCE_REACTION_TYPE_HINT_CONTRADICTS_GRAPH",),
                error="INCOMPATIBLE_SOURCE_REACTION_TYPE_HINT",
            )
        source_types = (explicit,)
    index = load_weak_label_index(records_path)
    rows = index.select_types(source_types)
    ranked, compatible_count, excluded_count = _rank_recipes(
        hint.participants,
        rows,
        query_compatibility=_query_compatibility_payload(analysis),
    )
    recommendations = (
        _diverse_selection(ranked, size=top_k)
        if normalized_mode == "screening"
        else ranked[:top_k]
    )
    warnings = [
        "WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED",
        "REACTIVE_SITES_MATCHED_AS_UNORDERED_PARTICIPANTS",
        "SOURCE_REACTION_TYPES_ARE_COMPATIBILITY_HINTS",
        "WEAK_LABEL_RECOMMENDATIONS_REQUIRE_EXPERT_REVIEW",
    ]
    if len(source_types) > 1:
        warnings.append("AMBIGUOUS_SOURCE_REACTION_TYPES_RETAINED")
    if normalized_mode == "screening":
        warnings.append("SCREENING_ARRAY_DIVERSIFIED_ACROSS_CONDITION_COMPONENTS")
    return WeakLabelRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=bool(recommendations),
        recommendation_mode=result_mode,
        reaction_type_hint_id=hint.hint_id,
        reaction_type_id=hint.type_id,
        source_reaction_type_candidates=source_types,
        query_participants=tuple(asdict(value) for value in hint.participants),
        candidate_count=len(rows),
        compatible_candidate_count=compatible_count,
        excluded_candidate_count=excluded_count,
        recipe_count=len(ranked),
        recommendations=recommendations,
        warnings=tuple(warnings),
        error=None if recommendations else "NO_COMPATIBLE_WEAK_LABEL_RECIPE",
    )


def generate_weak_label_screening_array(
    reaction_smiles: str,
    *,
    records_path: str | Path = DEFAULT_WEAK_LABEL_RECORDS_PATH,
    array_size: int = 24,
    source_reaction_type_hint: str | None = None,
) -> WeakLabelRecommendationResult:
    """Generate a deterministic diversity-selected weak-label condition array."""
    return recommend_weak_label_conditions(
        reaction_smiles,
        records_path=records_path,
        top_k=array_size,
        mode="screening",
        source_reaction_type_hint=source_reaction_type_hint,
    )


__all__ = [
    "DEFAULT_WEAK_LABEL_RECORDS_PATH",
    "generate_weak_label_screening_array",
    "load_weak_label_retrieval_rules",
    "recommend_weak_label_conditions",
    "validate_weak_label_retrieval_rules",
]
