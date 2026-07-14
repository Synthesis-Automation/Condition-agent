"""Public type-agnostic condition recommendation API."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

from reactive_taxonomy import featurize_reaction

from .compatibility import CompatibilityAssessment
from .generic_indexing import GenericIndexedReaction, load_generic_index
from .generic_retrieval import (
    generic_signature_similarity,
    load_generic_retrieval_rules,
    retrieve_compatible_generic_pool,
)
from .models import GenericConditionRecommendation, GenericRecommendationResult


def _explanation(
    level: str,
    components: Dict[str, float],
    support: int,
    datasets: int,
    compatibility_score: float,
) -> Tuple[str, ...]:
    notes = [f"Retrieved at {level.replace('_', ' ')} level"]
    labels = {
        "edit_topology": "bond edits",
        "handles": "reactive handles",
        "contexts": "attachment contexts",
        "environment": "local environments",
        "spectators": "spectator groups",
        "transformation": "transformation class",
        "family": "named family",
    }
    matches = [labels[name] for name, score in components.items() if score >= 0.999]
    if matches:
        notes.append("Exact match: " + ", ".join(matches))
    notes.append(f"Recipe compatibility score: {compatibility_score:.2f}")
    notes.append(f"Supported by {support} precedent(s) from {datasets} dataset(s)")
    return tuple(notes)


def _aggregate(
    query: Dict[str, Any],
    assessed_pool: Tuple[
        Tuple[GenericIndexedReaction, CompatibilityAssessment], ...
    ],
    level: str,
    top_k: int,
) -> Tuple[GenericConditionRecommendation, ...]:
    rules = load_generic_retrieval_rules()
    maximum = int(rules["maximum_neighbors"])
    scored = []
    for row, assessment in assessed_pool:
        similarity, components = generic_signature_similarity(query, row.signature)
        scored.append((similarity, row, components, assessment))
    scored.sort(
        key=lambda item: (-item[0], -item[1].yield_pct, item[1].reaction_id)
    )
    groups: Dict[
        str,
        List[
            Tuple[
                float,
                GenericIndexedReaction,
                Dict[str, float],
                CompatibilityAssessment,
            ]
        ],
    ] = defaultdict(list)
    for item in scored[:maximum]:
        groups[item[1].recipe_id].append(item)
    pool_mean = sum(row.yield_pct for row, _ in assessed_pool) / len(assessed_pool)
    prior_strength = float(rules["yield_prior_strength"])
    ranking = rules["ranking_weights"]
    ranked = []
    for recipe_id, members in groups.items():
        similarity_weight = sum(max(item[0], 0.05) for item in members)
        weighted_yield = sum(
            max(similarity, 0.05) * row.yield_pct
            for similarity, row, _, _ in members
        )
        expected_yield = (
            weighted_yield + prior_strength * pool_mean
        ) / (similarity_weight + prior_strength)
        similarity_score = sum(item[0] for item in members) / len(members)
        compatibility_score = sum(item[3].score for item in members) / len(members)
        support_score = min(1.0, math.log1p(len(members)) / math.log1p(10))
        datasets = {item[1].source_dataset for item in members if item[1].source_dataset}
        diversity_score = min(1.0, len(datasets) / 3.0)
        score = (
            float(ranking["similarity"]) * similarity_score
            + float(ranking["yield"]) * expected_yield / 100.0
            + float(ranking["support"]) * support_score
            + float(ranking["dataset_diversity"]) * diversity_score
            + float(ranking["compatibility"]) * compatibility_score
        )
        ranked.append(
            (
                score,
                expected_yield,
                similarity_score,
                recipe_id,
                members,
                len(datasets),
                compatibility_score,
            )
        )
    ranked.sort(key=lambda item: (-item[0], -item[1], item[3]))
    recommendations = []
    for rank, item in enumerate(ranked[:top_k], start=1):
        (
            score,
            expected_yield,
            similarity_score,
            recipe_id,
            members,
            datasets,
            compatibility_score,
        ) = item
        best = members[0]
        uncertain = any(member[1].condition_uncertain for member in members)
        cautions = []
        if uncertain:
            cautions.append("Condition identity or contextual role is uncertain")
        if level.endswith("limited_support"):
            cautions.append("Retrieval pool is below the configured support threshold")
        compatibility_evidence = tuple(
            sorted(
                {
                    message
                    for member in members
                    for message in member[3].evidence
                }
            )
        )
        cautions.extend(compatibility_evidence)
        recommendations.append(
            GenericConditionRecommendation(
                rank=rank,
                recipe_id=recipe_id,
                resolved_recipe=best[1].resolved_recipe,
                score=round(score, 6),
                similarity_score=round(similarity_score, 6),
                compatibility_score=round(compatibility_score, 6),
                expected_yield_pct=round(expected_yield, 2),
                support=len(members),
                dataset_support=datasets,
                retrieval_level=level,
                precedent_reaction_ids=tuple(
                    member[1].reaction_id for member in members[:5]
                ),
                explanation=_explanation(
                    level,
                    best[2],
                    len(members),
                    datasets,
                    compatibility_score,
                ),
                compatibility_evidence=compatibility_evidence,
                cautions=tuple(cautions),
            )
        )
    return tuple(recommendations)


def recommend_generic_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = "results/generic_conversion/records.jsonl",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
) -> GenericRecommendationResult:
    """Recommend canonical recipes without requiring a named reaction family."""
    if top_k < 1:
        return GenericRecommendationResult(
            reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error=analysis.error or "INVALID_REACTION",
        )
    if analysis.reaction_signature is None:
        return GenericRecommendationResult(
            reaction_smiles,
            False,
            error="QUERY_HAS_NO_USABLE_REACTION_SIGNATURE",
        )
    signature = asdict(analysis.reaction_signature)
    index = load_generic_index(records_path)
    if not index.rows:
        return GenericRecommendationResult(
            reaction_smiles, False, error="EMPTY_GENERIC_INDEX"
        )
    level, compatible_pool, candidate_count, excluded_count = (
        retrieve_compatible_generic_pool(
            signature, index, minimum_pool_size=minimum_pool_size
        )
    )
    if not compatible_pool:
        compatibility_failure = level == "no_compatible_condition_precedent"
        return GenericRecommendationResult(
            query_reaction_smiles=reaction_smiles,
            valid=False,
            query_signature_id=str(signature.get("signature_id") or ""),
            named_family=signature.get("named_family"),
            transformation_class=signature.get("transformation_class"),
            retrieval_level=level,
            candidate_count=candidate_count,
            excluded_candidate_count=excluded_count,
            warnings=("ALL_RETRIEVED_RECIPES_FAILED_COMPATIBILITY",)
            if compatibility_failure
            else (),
            error="NO_COMPATIBLE_CONDITION_PRECEDENT"
            if compatibility_failure
            else "NO_CHEMICALLY_COMPATIBLE_PRECEDENT",
        )
    warnings = []
    if level.endswith("limited_support"):
        warnings.append("LIMITED_PRECEDENT_SUPPORT")
    if level not in {"exact_signature", "handle_signature", "named_family"}:
        warnings.append("TYPE_AGNOSTIC_FALLBACK_USED")
    if excluded_count:
        warnings.append(f"INCOMPATIBLE_PRECEDENTS_EXCLUDED:{excluded_count}")
    return GenericRecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_signature_id=str(signature.get("signature_id") or ""),
        named_family=signature.get("named_family"),
        transformation_class=signature.get("transformation_class"),
        retrieval_level=level,
        candidate_count=candidate_count,
        compatible_candidate_count=len(compatible_pool),
        excluded_candidate_count=excluded_count,
        recommendations=_aggregate(signature, compatible_pool, level, top_k),
        warnings=tuple(warnings),
    )


__all__ = ["recommend_generic_conditions"]
