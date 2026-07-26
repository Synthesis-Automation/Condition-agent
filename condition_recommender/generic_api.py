"""Public type-agnostic condition recommendation API."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

from reactive_taxonomy import featurize_reaction

from .compatibility import CompatibilityAssessment
from .generic_indexing import (
    GenericIndexedReaction,
    GenericReactionIndex,
    load_generic_index,
)
from .generic_retrieval import (
    generic_signature_similarity,
    load_generic_retrieval_rules,
    reaction_scope,
    retrieve_compatible_generic_pool,
)
from .models import GenericConditionRecommendation, GenericRecommendationResult


def _explanation(
    level: str,
    components: Dict[str, float],
    support: int,
    references: int,
    datasets: int,
    compatibility_score: float,
    query_scope: str,
    precedent_scope: str,
) -> Tuple[str, ...]:
    notes = [f"Retrieved at {level.replace('_', ' ')} level"]
    labels = {
        "edit_topology": "bond edits",
        "reaction_events": "reaction-event multiset",
        "reaction_topology": "reaction topology",
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
    if query_scope and precedent_scope and query_scope != precedent_scope:
        notes.append(
            f"Reaction-scope mismatch: query {query_scope}; precedent {precedent_scope}"
        )
    notes.append(f"Recipe compatibility score: {compatibility_score:.2f}")
    notes.append(
        f"Supported by {support} distinct reaction(s), "
        f"{references} independent reference(s), and {datasets} dataset(s)"
    )
    return tuple(notes)


def _evidence_unit(row: GenericIndexedReaction) -> str:
    """Return the publication-level unit used to limit correlated evidence."""
    if row.reference_id:
        return f"reference:{row.reference_id}"
    return "reaction:" + (
        row.canonical_reaction_id or row.observation_id or row.reaction_id
    )


def _best_by_evidence_unit(
    members: List[
        Tuple[
            float,
            GenericIndexedReaction,
            Dict[str, float],
            CompatibilityAssessment,
        ]
    ],
) -> list[
    Tuple[
        float,
        GenericIndexedReaction,
        Dict[str, float],
        CompatibilityAssessment,
    ]
]:
    selected = {}
    for member in sorted(
        members,
        key=lambda item: (
            -item[0],
            -(item[1].yield_pct if item[1].yield_pct is not None else -1.0),
            item[1].reaction_id,
        ),
    ):
        selected.setdefault(_evidence_unit(member[1]), member)
    return list(selected.values())


def _aggregate(
    query: Dict[str, Any],
    assessed_pool: Tuple[Tuple[GenericIndexedReaction, CompatibilityAssessment], ...],
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
        key=lambda item: (
            -item[0],
            -(item[1].yield_pct if item[1].yield_pct is not None else -1.0),
            item[1].reaction_id,
        )
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
    selected_units = set()
    for item in scored:
        unit = (item[1].recipe_core_id, _evidence_unit(item[1]))
        if unit not in selected_units and len(selected_units) >= maximum:
            continue
        selected_units.add(unit)
        groups[item[1].recipe_core_id].append(item)
    pool_yields = [
        member[1].yield_pct
        for member in _best_by_evidence_unit(
            [
                (
                    generic_signature_similarity(query, row.signature)[0],
                    row,
                    {},
                    assessment,
                )
                for row, assessment in assessed_pool
                if row.yield_pct is not None
            ]
        )
    ]
    pool_mean = (
        sum(value for value in pool_yields if value is not None) / len(pool_yields)
        if pool_yields
        else None
    )
    prior_strength = float(rules["yield_prior_strength"])
    ranking = rules["ranking_weights"]
    ranked = []
    for recipe_core_id, members in groups.items():
        independent_members = _best_by_evidence_unit(members)
        outcome_members = _best_by_evidence_unit(
            [member for member in members if member[1].yield_pct is not None]
        )
        if outcome_members and pool_mean is not None:
            similarity_weight = sum(max(item[0], 0.05) for item in outcome_members)
            weighted_yield = sum(
                max(similarity, 0.05) * float(row.yield_pct)
                for similarity, row, _, _ in outcome_members
            )
            expected_yield = (weighted_yield + prior_strength * pool_mean) / (
                similarity_weight + prior_strength
            )
        else:
            expected_yield = None
        similarity_score = sum(item[0] for item in independent_members) / len(
            independent_members
        )
        compatibility_score = sum(
            item[3].score for item in independent_members
        ) / len(independent_members)
        canonical_reactions = {
            item[1].canonical_reaction_id
            or item[1].observation_id
            or item[1].reaction_id
            for item in members
        }
        support = len(canonical_reactions)
        support_score = min(1.0, math.log1p(support) / math.log1p(10))
        datasets = {
            item[1].source_dataset for item in members if item[1].source_dataset
        }
        references = {
            item[1].reference_id for item in members if item[1].reference_id
        }
        condition_series = {
            item[1].reference_condition_series_id
            for item in members
            if item[1].reference_condition_series_id
        }
        recipe_variants = {item[1].recipe_id for item in members}
        diversity_score = min(1.0, len(datasets) / 3.0)
        score = (
            float(ranking["similarity"]) * similarity_score
            + float(ranking["yield"])
            * ((expected_yield or 0.0) / 100.0)
            + float(ranking["support"]) * support_score
            + float(ranking["dataset_diversity"]) * diversity_score
            + float(ranking["compatibility"]) * compatibility_score
        )
        ranked.append(
            (
                score,
                expected_yield,
                similarity_score,
                recipe_core_id,
                members,
                len(datasets),
                compatibility_score,
                support,
                len(references),
                len(condition_series),
                tuple(sorted(recipe_variants)),
            )
        )
    ranked.sort(
        key=lambda item: (
            -item[0],
            -(item[1] if item[1] is not None else -1.0),
            item[3],
        )
    )
    recommendations = []
    for rank, item in enumerate(ranked[:top_k], start=1):
        (
            score,
            expected_yield,
            similarity_score,
            recipe_core_id,
            members,
            datasets,
            compatibility_score,
            support,
            reference_support,
            condition_series_support,
            recipe_variants,
        ) = item
        best = members[0]
        uncertain = any(member[1].condition_uncertain for member in members)
        cautions = []
        if uncertain:
            cautions.append("Condition identity or contextual role is uncertain")
        if level.endswith("limited_support"):
            cautions.append("Retrieval pool is below the configured support threshold")
        if len(recipe_variants) > 1:
            cautions.append(
                f"Recipe core has {len(recipe_variants)} operating-condition variants"
            )
        if len(members) > len(_best_by_evidence_unit(members)):
            cautions.append(
                "Repeated observations from the same reference count as one "
                "independent evidence unit"
            )
        query_scope = reaction_scope(query)
        precedent_scopes = {
            reaction_scope(member[1].signature)
            for member in members
            if reaction_scope(member[1].signature)
        }
        mismatched_scopes = sorted(
            scope for scope in precedent_scopes if query_scope and scope != query_scope
        )
        if mismatched_scopes:
            cautions.append(
                "Reaction-scope mismatch: query "
                f"{query_scope}; precedent {', '.join(mismatched_scopes)}"
            )
        compatibility_evidence = tuple(
            sorted({message for member in members for message in member[3].evidence})
        )
        cautions.extend(compatibility_evidence)
        recommendations.append(
            GenericConditionRecommendation(
                rank=rank,
                recipe_id=best[1].recipe_id,
                recipe_core_id=recipe_core_id,
                recipe_variant_ids=recipe_variants,
                resolved_recipe=best[1].resolved_recipe,
                score=round(score, 6),
                similarity_score=round(similarity_score, 6),
                compatibility_score=round(compatibility_score, 6),
                expected_yield_pct=(
                    round(expected_yield, 2) if expected_yield is not None else None
                ),
                support=support,
                observation_support=len(members),
                reference_support=reference_support,
                condition_series_support=condition_series_support,
                dataset_support=datasets,
                retrieval_level=level,
                precedent_reaction_ids=tuple(
                    member[1].reaction_id for member in members[:5]
                ),
                explanation=_explanation(
                    level,
                    best[2],
                    support,
                    reference_support,
                    datasets,
                    compatibility_score,
                    query_scope,
                    reaction_scope(best[1].signature),
                ),
                compatibility_evidence=compatibility_evidence,
                cautions=tuple(cautions),
            )
        )
    return tuple(recommendations)


def recommend_indexed_signature(
    signature: Dict[str, Any],
    index: GenericReactionIndex,
    *,
    query_reaction_smiles: str = "",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
) -> GenericRecommendationResult:
    """Recommend from an existing signature and index without re-featurization."""
    if top_k < 1:
        return GenericRecommendationResult(
            query_reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE"
        )
    if not index.rows:
        return GenericRecommendationResult(
            query_reaction_smiles, False, error="EMPTY_GENERIC_INDEX"
        )
    if str(signature.get("schema_version") or "") != (
        index.reaction_signature_schema_version
    ):
        return GenericRecommendationResult(
            query_reaction_smiles,
            False,
            error="INCOMPATIBLE_REACTION_SIGNATURE_SCHEMA",
        )
    query_definitions = signature.get("definition_versions") or {}
    if (
        not isinstance(query_definitions, dict)
        or tuple(
            sorted((str(key), str(value)) for key, value in query_definitions.items())
        )
        != index.taxonomy_definition_versions
    ):
        return GenericRecommendationResult(
            query_reaction_smiles,
            False,
            error="INCOMPATIBLE_REACTION_TAXONOMY_DEFINITIONS",
        )
    level, compatible_pool, candidate_count, excluded_count = (
        retrieve_compatible_generic_pool(
            signature, index, minimum_pool_size=minimum_pool_size
        )
    )
    if not compatible_pool:
        compatibility_failure = level == "no_compatible_condition_precedent"
        return GenericRecommendationResult(
            query_reaction_smiles=query_reaction_smiles,
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
    recommendations = _aggregate(signature, compatible_pool, level, top_k)
    if any(
        caution.startswith("Reaction-scope mismatch:")
        for recommendation in recommendations
        for caution in recommendation.cautions
    ):
        warnings.append("REACTION_TOPOLOGY_FALLBACK_USED")
    return GenericRecommendationResult(
        query_reaction_smiles=query_reaction_smiles,
        valid=True,
        query_signature_id=str(signature.get("signature_id") or ""),
        named_family=signature.get("named_family"),
        transformation_class=signature.get("transformation_class"),
        retrieval_level=level,
        candidate_count=candidate_count,
        compatible_candidate_count=len(compatible_pool),
        excluded_candidate_count=excluded_count,
        recommendations=recommendations,
        warnings=tuple(warnings),
    )


def recommend_generic_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = "results/generic_conversion/records.jsonl",
    top_k: int = 5,
    minimum_pool_size: int | None = None,
) -> GenericRecommendationResult:
    """Featurize a reaction and recommend canonical resolved recipes."""
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
    return recommend_indexed_signature(
        signature,
        index,
        query_reaction_smiles=reaction_smiles,
        top_k=top_k,
        minimum_pool_size=minimum_pool_size,
    )


__all__ = ["recommend_generic_conditions", "recommend_indexed_signature"]
