"""Reference-aware aggregation and auditable ranking of condition recipes."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Tuple

from .compatibility import CompatibilityAssessment
from .generic_indexing import GenericIndexedReaction
from .models import (
    ChemistRankingPreferences,
    GenericConditionRecommendation,
    RecommendationScoreTrace,
)
from .preference_scoring import (
    aggregate_partner_category_evidence,
    assess_functional_group_tolerance,
)
from .protocols import protocol_draft_for_reaction
from .ranking_preferences import (
    RANKING_COMPONENTS,
    normalize_ranking_weights,
    resolve_ranking_preferences,
)
from .similarity import (
    SimilarityAssessment,
    assess_signature_similarity,
    reaction_scope,
)
from .support import (
    evidence_unit,
    load_evidence_support_rules,
    summarize_evidence_support,
)

_RULES_PATH = Path(__file__).with_name("definitions") / "generic_ranking.v1.json"
_DEFINITION_ID = "generic_ranking.v1"
_SCHEMA_VERSION = "1.1"
_RANKING_COMPONENTS = RANKING_COMPONENTS


@dataclass(frozen=True)
class ScoredPrecedent:
    similarity: SimilarityAssessment
    row: GenericIndexedReaction
    compatibility: CompatibilityAssessment


def _precedent_reaction_context(
    precedent: GenericIndexedReaction,
) -> Dict[str, Any]:
    """Return display context from one indexed structural precedent."""
    spectators = tuple(
        dict(group)
        for group in (precedent.signature.get("spectator_groups") or ())
        if isinstance(group, Mapping)
    )
    partners = tuple(
        dict(partner)
        for partner in (precedent.signature.get("partners") or ())
        if isinstance(partner, Mapping)
    )
    fallback = precedent.fallback_descriptor
    return {
        "reaction_id": precedent.reaction_id,
        "observation_id": precedent.observation_id,
        "reaction_smiles": precedent.reaction_smiles,
        "reference_id": precedent.reference_id,
        "reaction_label": precedent.reaction_label or None,
        "spectator_groups": spectators,
        "reaction_partners": partners,
        "fallback_evidence_mode": fallback.get("evidence_mode") or None,
        "observed_or_partial_edits": tuple(fallback.get("verified_edit_tokens") or ()),
        "reaction_center_core": tuple(
            fallback.get("reaction_center_core_tokens") or ()
        ),
        "partial_transformation_key": (
            fallback.get("partial_transformation_key") or None
        ),
        "fragment_source_support": tuple(precedent.fragment_source_support),
    }


@lru_cache(maxsize=1)
def load_generic_ranking_rules() -> Dict[str, Any]:
    """Load and validate the pilot ranking definition."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    rules["weights"] = _validated_ranking_weights(rules)
    return rules


def _validated_ranking_weights(
    rules: Mapping[str, Any],
) -> Dict[str, float]:
    if str(rules.get("schema_version") or "") != _SCHEMA_VERSION:
        raise ValueError("unsupported generic ranking definition schema")
    if str(rules.get("definition_id") or "") != _DEFINITION_ID:
        raise ValueError("unexpected generic ranking definition ID")
    if not str(rules.get("calibration_status") or "").strip():
        raise ValueError("generic ranking definition requires calibration status")
    normalized = _normalize_component_weights(
        rules.get("weights"),
        label="generic ranking",
    )
    baseline_weights = rules.get("evaluation_baseline_weights")
    if not isinstance(baseline_weights, Mapping) or set(baseline_weights) != {
        "transformation_prior",
        "legacy_pilot",
    }:
        raise ValueError("generic ranking baseline profiles are incomplete")
    for profile, weights in baseline_weights.items():
        _normalize_component_weights(
            weights,
            label=f"generic ranking baseline {profile}",
        )
    if rules.get("missing_yield_policy") != "renormalize_available_components":
        raise ValueError("unsupported generic ranking missing-yield policy")
    if int(rules["maximum_independent_neighbors"]) < 1:
        raise ValueError("maximum_independent_neighbors must be positive")
    if float(rules["yield_prior_strength"]) < 0.0:
        raise ValueError("yield_prior_strength must be non-negative")
    minimum_similarity_weight = float(rules["minimum_similarity_yield_weight"])
    if not 0.0 < minimum_similarity_weight <= 1.0:
        raise ValueError("minimum_similarity_yield_weight must be in (0, 1]")
    saturation = rules.get("support_saturation") or {}
    if any(
        int(saturation.get(name) or 0) < 1
        for name in ("independent_evidence", "canonical_reactions", "datasets")
    ):
        raise ValueError("generic ranking support saturation must be positive")
    return normalized


def _normalize_component_weights(
    weights: Any,
    *,
    label: str,
) -> Dict[str, float]:
    if not isinstance(weights, Mapping):
        raise ValueError(f"{label} weights must be an object")
    return normalize_ranking_weights(weights, label=label)


def validate_generic_ranking_rules(rules: Mapping[str, Any]) -> None:
    """Validate a generic recipe-ranking definition without loading state."""
    _validated_ranking_weights(rules)


def _best_by_evidence_unit(
    members: Iterable[ScoredPrecedent],
) -> list[ScoredPrecedent]:
    selected = {}
    for member in sorted(
        members,
        key=lambda item: (
            -item.similarity.score,
            -(item.row.yield_pct if item.row.yield_pct is not None else -1.0),
            item.row.reaction_id,
        ),
    ):
        selected.setdefault(evidence_unit(member.row), member)
    return list(selected.values())


def _saturated_log(count: int, saturation: int) -> float:
    return min(1.0, math.log1p(count) / math.log1p(saturation))


def _mean(values: Iterable[float]) -> float:
    items = tuple(values)
    return sum(items) / len(items)


def _mean_similarity_trace(
    members: Iterable[ScoredPrecedent],
) -> tuple[Dict[str, float], Dict[str, float]]:
    values = tuple(members)
    component_names = tuple(values[0].similarity.components)
    components = {
        name: round(
            _mean(item.similarity.components[name] for item in values),
            6,
        )
        for name in component_names
    }
    contributions = {
        name: round(
            _mean(item.similarity.contributions[name] for item in values),
            6,
        )
        for name in component_names
    }
    return components, contributions


def _pool_yield_prior(scored: Iterable[ScoredPrecedent]) -> Optional[float]:
    known = _best_by_evidence_unit(
        member for member in scored if member.row.yield_pct is not None
    )
    if not known:
        return None
    return _mean(float(member.row.yield_pct) for member in known)


def _expected_yield(
    members: Iterable[ScoredPrecedent],
    *,
    pool_prior: Optional[float],
    prior_strength: float,
    minimum_similarity_weight: float,
) -> tuple[Optional[float], int]:
    known = _best_by_evidence_unit(
        member for member in members if member.row.yield_pct is not None
    )
    if not known or pool_prior is None:
        return None, 0
    similarity_weight = sum(
        max(member.similarity.score, minimum_similarity_weight) for member in known
    )
    weighted_yield = sum(
        max(member.similarity.score, minimum_similarity_weight)
        * float(member.row.yield_pct)
        for member in known
    )
    expected = (weighted_yield + prior_strength * pool_prior) / (
        similarity_weight + prior_strength
    )
    return expected, len(known)


def _weighted_score(
    components: Mapping[str, Optional[float]],
    base_weights: Mapping[str, float],
) -> tuple[float, Dict[str, float], Dict[str, float]]:
    available = {name: value for name, value in components.items() if value is not None}
    denominator = sum(float(base_weights[name]) for name in available)
    if denominator <= 0.0:
        empty = {name: 0.0 for name in _RANKING_COMPONENTS}
        return 0.0, empty, empty
    applied = {
        name: (float(base_weights[name]) / denominator if name in available else 0.0)
        for name in _RANKING_COMPONENTS
    }
    contributions = {
        name: applied[name] * float(components[name])
        if components[name] is not None
        else 0.0
        for name in _RANKING_COMPONENTS
    }
    rounded_contributions = {
        name: round(contributions[name], 6) for name in _RANKING_COMPONENTS
    }
    return (
        round(sum(rounded_contributions.values()), 6),
        {name: round(applied[name], 6) for name in _RANKING_COMPONENTS},
        rounded_contributions,
    )


def _explanation(
    *,
    level: str,
    similarity_components: Mapping[str, float],
    reaction_support: int,
    reference_support: int,
    dataset_support: int,
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
        "environment": "R-group and attachment environments",
        "spectators": "spectator groups",
        "transformation": "transformation class",
        "family": "named family",
    }
    matches = [
        labels.get(name, name.replace("_", " "))
        for name, score in similarity_components.items()
        if score >= 0.999
    ]
    if matches:
        notes.append("Exact match: " + ", ".join(matches))
    local_scores = [
        (name, score)
        for name, score in similarity_components.items()
        if name.startswith("reaction_center_")
    ]
    if local_scores:
        notes.append(
            "Reaction-center similarity: "
            + ", ".join(
                f"{name.removeprefix('reaction_center_').replace('_', ' ')} {score:.2f}"
                for name, score in local_scores
            )
        )
    environment_score = similarity_components.get("environment")
    if environment_score is not None and environment_score < 0.5:
        notes.append(
            "R-group/attachment environment differs from the query "
            f"(similarity {environment_score:.2f})"
        )
    if query_scope and precedent_scope and query_scope != precedent_scope:
        notes.append(
            f"Reaction-scope mismatch: query {query_scope}; precedent {precedent_scope}"
        )
    notes.append(f"Recipe compatibility score: {compatibility_score:.2f}")
    notes.append(
        f"Supported by {reaction_support} distinct reaction(s), "
        f"{reference_support} independent reference(s), and "
        f"{dataset_support} dataset(s)"
    )
    return tuple(notes)


def rank_condition_recipes(
    query: Mapping[str, Any],
    assessed_pool: Tuple[
        Tuple[GenericIndexedReaction, CompatibilityAssessment],
        ...,
    ],
    *,
    retrieval_level: str,
    top_k: int,
    ranking_profile: str = "default",
    ranking_weights: Optional[Mapping[str, float]] = None,
    ranking_preferences: ChemistRankingPreferences | None = None,
    similarity_assessor: Optional[
        Callable[
            [Mapping[str, Any], GenericIndexedReaction],
            SimilarityAssessment,
        ]
    ] = None,
    query_reaction_core: Mapping[str, Any] | None = None,
    query_reaction_smiles: str = "",
) -> Tuple[GenericConditionRecommendation, ...]:
    """Aggregate recipe cores and rank them with a complete score trace."""
    rules = load_generic_ranking_rules()
    resolved_preferences = resolve_ranking_preferences(ranking_preferences)
    if ranking_weights is not None:
        ranking_weights = _normalize_component_weights(
            ranking_weights,
            label=f"generic ranking override {ranking_profile}",
        )
    elif ranking_preferences is not None or ranking_profile == "default":
        ranking_weights = resolved_preferences.weights
    else:
        profiles = rules.get("evaluation_baseline_weights") or {}
        if ranking_profile not in profiles:
            raise ValueError(f"Unsupported ranking profile: {ranking_profile}")
        ranking_weights = _normalize_component_weights(
            profiles[ranking_profile],
            label=f"generic ranking baseline {ranking_profile}",
        )
    assessor = similarity_assessor
    if assessor is None:

        def assessor(
            query_value: Mapping[str, Any],
            row: GenericIndexedReaction,
        ) -> SimilarityAssessment:
            return assess_signature_similarity(
                query_value,
                row.signature,
                query_reaction_core=query_reaction_core,
                precedent_reaction_core=row.reaction_core,
            )

    scored = [
        ScoredPrecedent(
            similarity=assessor(query, row),
            row=row,
            compatibility=compatibility,
        )
        for row, compatibility in assessed_pool
    ]
    scored.sort(
        key=lambda item: (
            -item.similarity.score,
            -(item.row.yield_pct if item.row.yield_pct is not None else -1.0),
            item.row.reaction_id,
        )
    )

    maximum = int(rules["maximum_independent_neighbors"])
    groups: Dict[str, list[ScoredPrecedent]] = defaultdict(list)
    selected_units = set()
    for item in scored:
        unit = (item.row.recipe_core_id, evidence_unit(item.row))
        if unit not in selected_units and len(selected_units) >= maximum:
            continue
        selected_units.add(unit)
        groups[item.row.recipe_core_id].append(item)

    pool_prior = _pool_yield_prior(scored)
    saturation = rules["support_saturation"]
    ranking_rows = []
    for recipe_core_id, members in groups.items():
        independent = _best_by_evidence_unit(members)
        expected_yield, outcome_count = _expected_yield(
            members,
            pool_prior=pool_prior,
            prior_strength=float(rules["yield_prior_strength"]),
            minimum_similarity_weight=float(rules["minimum_similarity_yield_weight"]),
        )
        similarity_score = _mean(member.similarity.score for member in independent)
        compatibility_score = _mean(
            member.compatibility.score for member in independent
        )
        condition_certainty = _mean(
            float(
                not member.row.condition_uncertain
                and member.row.condition_stage_status != "unassigned_multistage"
            )
            for member in independent
        )
        support = summarize_evidence_support(member.row for member in members)
        partner_category_score, partner_category_evidence = (
            aggregate_partner_category_evidence(
                query_reaction_core,
                (
                    (evidence_unit(member.row), member.row.reaction_core)
                    for member in independent
                ),
            )
        )
        tolerance_score, tolerance_evidence = assess_functional_group_tolerance(
            query,
            (
                (evidence_unit(member.row), member.row.signature)
                for member in independent
            ),
        )
        components: Dict[str, Optional[float]] = {
            "similarity": similarity_score,
            "partner_category": partner_category_score,
            "functional_group_tolerance": tolerance_score,
            "yield": expected_yield / 100.0 if expected_yield is not None else None,
            "independent_support": _saturated_log(
                len(independent),
                int(saturation["independent_evidence"]),
            ),
            "reaction_breadth": _saturated_log(
                support.reaction_count,
                int(saturation["canonical_reactions"]),
            ),
            "dataset_diversity": min(
                1.0,
                support.dataset_count / int(saturation["datasets"]),
            ),
            "compatibility": compatibility_score,
            "condition_certainty": condition_certainty,
        }
        score, applied_weights, ranking_contributions = _weighted_score(
            components,
            ranking_weights,
        )
        default_score, _, default_contributions = _weighted_score(
            components,
            rules["weights"],
        )
        similarity_components, similarity_contributions = _mean_similarity_trace(
            independent
        )
        recipe_variants = tuple(sorted({member.row.recipe_id for member in members}))
        ranking_rows.append(
            (
                score,
                default_score,
                expected_yield,
                recipe_core_id,
                members,
                independent,
                support,
                similarity_score,
                compatibility_score,
                condition_certainty,
                outcome_count,
                recipe_variants,
                components,
                applied_weights,
                ranking_contributions,
                similarity_components,
                similarity_contributions,
                default_contributions,
                {
                    "partner_category": partner_category_evidence,
                    "functional_group_tolerance": tolerance_evidence,
                },
            )
        )

    default_order = sorted(
        ranking_rows,
        key=lambda item: (
            -item[1],
            -(item[2] if item[2] is not None else -1.0),
            item[3],
        ),
    )
    default_ranks = {
        item[3]: rank for rank, item in enumerate(default_order, start=1)
    }
    ranking_rows.sort(
        key=lambda item: (
            -item[0],
            -(item[2] if item[2] is not None else -1.0),
            item[3],
        )
    )
    recommendations = []
    query_scope = reaction_scope(query)
    for rank, item in enumerate(ranking_rows[:top_k], start=1):
        (
            score,
            default_score,
            expected_yield,
            recipe_core_id,
            members,
            independent,
            support,
            similarity_score,
            compatibility_score,
            condition_certainty,
            outcome_count,
            recipe_variants,
            components,
            applied_weights,
            ranking_contributions,
            similarity_components,
            similarity_contributions,
            default_ranking_contributions,
            factor_evidence,
        ) = item
        best = members[0]
        cautions = []
        if any(
            member.row.precedent_tier.value == "review_core"
            for member in members
        ):
            cautions.append(
                "Recipe support includes core-qualified review precedents; "
                "expert review is required"
            )
        if condition_certainty < 1.0:
            cautions.append(
                "Condition identity, contextual role, or stage assignment is uncertain"
            )
        if any(
            member.row.condition_stage_status == "unassigned_multistage"
            for member in members
        ):
            cautions.append(
                "Ingredient set is observed, but assignment to ordered "
                "reaction stages is unavailable"
            )
        if expected_yield is None:
            cautions.append(
                "No usable yield evidence; ranking excludes the outcome component"
            )
        if (
            ranking_weights.get("partner_category", 0.0) > 0.0
            and components["partner_category"] is None
        ):
            cautions.append(
                "Reactant-category evidence is unavailable; its weight was "
                "renormalized over available factors"
            )
        if (
            ranking_weights.get("functional_group_tolerance", 0.0) > 0.0
            and components["functional_group_tolerance"] is None
        ):
            cautions.append(
                "No unchanged query functional groups were available for "
                "direct tolerance scoring"
            )
        unknown_tolerance = tuple(
            str(group.get("chemist_label") or group.get("group_id") or "")
            for group in (
                factor_evidence.get("functional_group_tolerance", {}).get("groups")
                or ()
            )
            if group.get("status") == "unknown_not_tolerant"
        )
        if unknown_tolerance:
            cautions.append(
                "No direct tolerance precedent for unchanged group(s): "
                + ", ".join(unknown_tolerance)
            )
        if retrieval_level.endswith("limited_support"):
            cautions.append("Retrieval pool is below the configured support threshold")
        if retrieval_level.startswith("edit_graph_neighbors"):
            cautions.append(
                "Retrieval uses a chemistry-gated approximate edit graph; "
                "the exact net bond-edit signature differs from the query"
            )
        if len(recipe_variants) > 1:
            cautions.append(
                f"Recipe core has {len(recipe_variants)} operating-condition variants"
            )
        if len(members) > len(independent):
            cautions.append(
                "Repeated observations from the same reference count as one "
                "independent evidence unit"
            )
        if support.mapping_deduplicated_independent_count:
            cautions.append(
                "Mapping-equivalent unreferenced observations count as one "
                "independent evidence unit"
            )
        inferred_correspondence = sorted(
            {
                str(member.row.signature.get("evidence_quality") or "")
                for member in members
                if str(member.row.signature.get("evidence_quality") or "")
                in {
                    "fragmented_scaffold_correspondence",
                    "global_atom_correspondence",
                    "unique_scaffold_correspondence",
                }
            }
        )
        if inferred_correspondence:
            cautions.append(
                "Precedent reaction edits use deterministic inferred atom "
                "correspondence and retain review-level chemistry confidence"
            )
        if any(
            not member.row.signature
            and member.row.fallback_descriptor.get("evidence_mode")
            == "partial_product_correspondence"
            for member in members
        ):
            cautions.append(
                "Precedent transformation uses partial product correspondence; "
                "the source of the otherwise missing product fragment is not "
                "verified"
            )
        supported_fragment_sources = tuple(
            sorted(
                {
                    raw_identifier
                    for member in members
                    for support in member.row.fragment_source_support
                    if str(support.get("status") or "") == "supported"
                    for raw_identifier in (
                        support.get("component_raw_identifiers") or ()
                    )
                }
            )
        )
        if supported_fragment_sources:
            cautions.append(
                "Installed product fragment is supported by reported condition "
                "component(s), but direct atom origin is not verified: "
                + ", ".join(supported_fragment_sources)
            )
        precedent_scopes = {
            reaction_scope(member.row.signature)
            for member in members
            if reaction_scope(member.row.signature)
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
            sorted(
                {
                    message
                    for member in members
                    for message in member.compatibility.evidence
                }
            )
        )
        cautions.extend(compatibility_evidence)
        score_trace = RecommendationScoreTrace(
            similarity_components=similarity_components,
            similarity_contributions=similarity_contributions,
            ranking_components={
                name: (round(value, 6) if value is not None else None)
                for name, value in components.items()
            },
            ranking_contributions=ranking_contributions,
            applied_ranking_weights=applied_weights,
            independent_evidence_count=len(independent),
            observed_outcome_count=outcome_count,
            pool_yield_prior_pct=(
                round(pool_prior, 6) if pool_prior is not None else None
            ),
            definition_versions={
                best.similarity.definition_id: best.similarity.definition_version,
                str(rules["definition_id"]): str(rules["schema_version"]),
                str(load_evidence_support_rules()["definition_id"]): str(
                    load_evidence_support_rules()["definition_version"]
                ),
                best.compatibility.definition_id: (
                    best.compatibility.definition_version
                ),
                resolved_preferences.definition_id: (
                    resolved_preferences.definition_version
                ),
            },
            ranking_profile=(
                resolved_preferences.profile_id
                if ranking_preferences is not None or ranking_profile == "default"
                else ranking_profile
            ),
            default_ranking_contributions=default_ranking_contributions,
        )
        recommendations.append(
            GenericConditionRecommendation(
                rank=rank,
                recipe_id=best.row.recipe_id,
                recipe_core_id=recipe_core_id,
                recipe_variant_ids=recipe_variants,
                resolved_recipe=best.row.resolved_recipe,
                synthesis_protocol=protocol_draft_for_reaction(
                    best.row.resolved_recipe,
                    query_reaction_smiles,
                ).to_dict(),
                score=round(score, 6),
                similarity_score=round(similarity_score, 6),
                compatibility_score=round(compatibility_score, 6),
                expected_yield_pct=(
                    round(expected_yield, 2) if expected_yield is not None else None
                ),
                support=support.reaction_count,
                observation_support=support.observation_count,
                reference_support=support.reference_count,
                condition_series_support=support.condition_series_count,
                dataset_support=support.dataset_count,
                retrieval_level=retrieval_level,
                precedent_reaction_ids=tuple(
                    member.row.reaction_id for member in members[:5]
                ),
                precedent_reaction_smiles=tuple(
                    member.row.reaction_smiles for member in members[:5]
                ),
                precedent_reaction_contexts=tuple(
                    _precedent_reaction_context(member.row) for member in members[:5]
                ),
                precedent_reference_ids=tuple(
                    sorted(
                        {
                            member.row.reference_id
                            for member in members
                            if member.row.reference_id
                        }
                    )[:5]
                ),
                explanation=_explanation(
                    level=retrieval_level,
                    similarity_components=similarity_components,
                    reaction_support=support.reaction_count,
                    reference_support=support.reference_count,
                    dataset_support=support.dataset_count,
                    compatibility_score=compatibility_score,
                    query_scope=query_scope,
                    precedent_scope=reaction_scope(best.row.signature),
                ),
                score_trace=score_trace,
                default_rank=default_ranks[recipe_core_id],
                default_score=round(default_score, 6),
                rank_change=default_ranks[recipe_core_id] - rank,
                factor_evidence=factor_evidence,
                compatibility_evidence=compatibility_evidence,
                cautions=tuple(cautions),
            )
        )
    return tuple(recommendations)


__all__ = [
    "load_generic_ranking_rules",
    "rank_condition_recipes",
    "validate_generic_ranking_rules",
]
