"""Shared candidate qualification inside the canonical condition engine."""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Mapping

from condition_registry import ConditionConstraintSet
from reactive_taxonomy.shared_reaction_core import (
    LEVELS,
    build_shared_reaction_core,
    compare_reaction_cores,
)

from .compatibility import filter_compatible_precedents
from .generic_indexing import GenericReactionIndex
from .models import (
    ChemistRankingPreferences,
    GenericRecommendationResult,
    RetrievalLevelTrace,
)
from .recipe_ranking import rank_condition_recipes
from .shared_core_index import SharedCoreIndex
from .support import summarize_evidence_support


def load_shared_retrieval_rules() -> dict[str, Any]:
    """Validate explicit candidate budgets and level display semantics."""
    rules = json.loads(
        (
            Path(__file__).with_name("definitions") / "shared_core_retrieval.v1.json"
        ).read_text(encoding="utf-8")
    )
    if (
        rules.get("definition_id") != "shared_core_retrieval.v1"
        or rules.get("schema_version") != "1.0"
        or rules.get("levels") != ["whole_reaction", *LEVELS]
        or len(rules.get("labels", [])) != 4
        or rules.get("status") != "experimental_pending_independent_review"
    ):
        raise ValueError("invalid shared core retrieval definition")
    for name in ("minimum_independent_support", "candidate_limit"):
        if type(rules.get(name)) is not int or rules[name] < 1:
            raise ValueError(f"invalid {name}")
    return rules


def recommend_from_shared_core(
    signature: Mapping[str, Any],
    index: GenericReactionIndex,
    shared_index: SharedCoreIndex,
    *,
    reaction_core: Mapping[str, Any],
    query_reaction_smiles: str,
    top_k: int,
    minimum_pool_size: int | None,
    search_scope: str,
    preferred_reaction_ids: tuple[str, ...] = (),
    ranking_preferences: ChemistRankingPreferences | None = None,
    ranking_weights: Mapping[str, float] | None = None,
    molecular_features: Mapping[str, Any] | None = None,
    condition_constraints: ConditionConstraintSet | None = None,
) -> GenericRecommendationResult:
    """Union direct and product-side seeds, then compare and filter once.

    Product-side retrieval reuses persisted observed product cores. Explicit
    retro operator precedent IDs enter through the same comparison, with no
    privileged ranking and no import of or recursion into the retro planner.
    """
    rules = load_shared_retrieval_rules()
    shared_index.validate_binding(index)
    minimum = (
        rules["minimum_independent_support"]
        if minimum_pool_size is None
        else minimum_pool_size
    )
    if minimum < 1:
        raise ValueError("minimum_pool_size must be positive")
    query = build_shared_reaction_core(query_reaction_smiles, signature, reaction_core)
    base = dict(
        query_reaction_smiles=query_reaction_smiles,
        query_signature_id=signature.get("signature_id"),
        query_reaction_core_id=query.observation_core_id,
        named_family=signature.get("named_family"),
        transformation_class=signature.get("transformation_class"),
        recommendation_mode="experimental_shared_core",
        search_scope=search_scope,
        retrieval_definition_version="shared_core_retrieval.v1@1.0;"
        + query.definition_hash,
        warnings=(
            "EXPERIMENTAL_SHARED_CORE_PENDING_INDEPENDENT_REVIEW",
            *query.warnings,
        ),
    )
    if not query.levels:
        return GenericRecommendationResult(
            **base,
            valid=False,
            error="SHARED_CORE_QUERY_UNAVAILABLE",
            shared_core_trace=(
                {"query_unavailable_reasons": query.unavailable_reasons},
            ),
        )
    limit = rules["candidate_limit"]
    channels: dict[int, set[str]] = defaultdict(set)
    audited: dict[int, dict[str, Any]] = {}
    eligible: dict[int, tuple[Any, Any]] = {}
    projections = {}
    comparisons = {}
    diagnostics: list[dict[str, Any]] = []
    traces = []
    truncated = False
    maximum_level = (
        "observed_local" if search_scope == "same_handle" else "retained_typed"
    )

    def add(kind: str, key: str, channel: str) -> None:
        nonlocal truncated
        positions, capped = shared_index.lookup(kind, key, limit)
        truncated |= capped
        for position in positions:
            if position not in channels and len(channels) >= limit:
                truncated = True
                continue
            channels[position].add(channel)

    def qualify() -> None:
        pending = sorted(set(channels) - set(audited))
        for position, row in zip(pending, index.select(pending)):
            projection = shared_index.projection(position, row)
            comparison = compare_reaction_cores(
                query, projection, maximum_level=maximum_level
            )
            projections[position], comparisons[position] = projection, comparison
            record = dict(
                position=position,
                reaction_id=row.reaction_id,
                observation_id=row.observation_id,
                comparison=asdict(comparison),
            )
            if comparison.eligible:
                accepted, excluded = filter_compatible_precedents(
                    signature,
                    (row,),
                    condition_constraints=condition_constraints,
                )
                assessment = (accepted or excluded)[0][1]
                record["condition_compatibility"] = asdict(assessment)
                if accepted:
                    eligible[position] = accepted[0]
            audited[position] = record

    def support() -> int:
        return summarize_evidence_support(
            tuple(row for row, _ in eligible.values())
        ).independent_count

    # Lookup each exact/local tier first. Abstraction is independent of recipe
    # preference and of how many rows the UI can display.
    direct = [("whole_reaction", query.reaction_identity)]
    direct += [(level.level, level.key) for level in query.levels]
    for kind, key in direct:
        if search_scope == "same_handle" and kind not in {
            "whole_reaction",
            "observed_local",
        }:
            break
        add(kind, key, "direct")
        qualify()
        count = support()
        traces.append(
            RetrievalLevelTrace(
                "shared_core:" + kind,
                len(channels),
                summarize_evidence_support(
                    index.select(sorted(channels))
                ).independent_count,
                len(eligible),
                count,
                len(channels) - len(eligible),
                minimum,
                "selected"
                if count >= minimum
                else "selected_limited_support"
                if eligible
                else "empty",
            )
        )
        if count >= minimum and search_scope != "broad":
            break

    if support() < minimum or search_scope == "broad":
        add("product_identity", query.product_identity, "product_side")
        for level in query.levels:
            add("product_core", level.product_side_key, "product_side")
    # Retro seeds share the budget and all original-query gates. Missing links
    # stay visible, including operators whose observations have no condition row.
    for reaction_id in dict.fromkeys(preferred_reaction_ids):
        positions = index.reaction_ids.get(reaction_id, ())
        if not positions:
            diagnostics.append(
                {
                    "reaction_id": reaction_id,
                    "reason": "RETRO_PRECEDENT_NOT_IN_CONDITION_INDEX",
                }
            )
        for position in positions:
            if position not in channels and len(channels) >= limit:
                truncated = True
                continue
            channels[position].add("retro_precedent")
    qualify()
    for position in sorted(audited):
        diagnostics.append(
            {**audited[position], "channels": sorted(channels[position])}
        )

    # Preserve input/source context when aggregating recipes. This conservative
    # first contract never makes a single recipe out of different required inputs.
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for position in eligible:
        groups[
            (comparisons[position].level, projections[position].input_identity)
        ].append(position)
    recommendations = []
    for (level, _), positions in sorted(groups.items()):
        ranked = rank_condition_recipes(
            signature,
            tuple(eligible[p] for p in positions),
            retrieval_level="shared_core:" + level,
            top_k=top_k,
            ranking_preferences=ranking_preferences,
            ranking_weights=ranking_weights,
            query_reaction_core=reaction_core,
            query_reaction_smiles=query_reaction_smiles,
            query_molecular_features=molecular_features,
        )
        ordinal = rules["levels"].index(level) + 1
        for item in ranked:
            members = [
                p
                for p in positions
                if eligible[p][0].reaction_id in item.precedent_reaction_ids
                and eligible[p][0].recipe_core_id == item.recipe_core_id
            ]
            differences = tuple(
                dict.fromkeys(d for p in members for d in comparisons[p].differences)
            )
            found_by = tuple(sorted({c for p in members for c in channels[p]}))
            requirements = tuple(
                {
                    "reaction_id": eligible[p][0].reaction_id,
                    "observation_id": eligible[p][0].observation_id,
                    "reaction_smiles": eligible[p][0].reaction_smiles,
                    "observed_ports": projections[p].realization_details,
                    "evidence_status": projections[p].evidence_status,
                }
                for p in members
            )
            recommendations.append(
                replace(
                    item,
                    match_level=ordinal,
                    match_label=rules["labels"][ordinal - 1],
                    match_namespace="shared_reaction_core.v1",
                    match_details=differences,
                    evidence_relation="same_setup"
                    if level == "whole_reaction"
                    else "analogue_evidence",
                    candidate_channels=found_by,
                    source_input_requirements=requirements,
                    explanation=(
                        *item.explanation,
                        "Shared graph comparison qualified this precedent against the original query",
                        "Retrieved through: " + ", ".join(found_by),
                        *differences,
                    ),
                    cautions=tuple(
                        dict.fromkeys(
                            (
                                *item.cautions,
                                *query.warnings,
                                *(w for p in members for w in projections[p].warnings),
                                *(
                                    (
                                        "Reported conditions are analogue evidence; transfer to the supplied inputs is unvalidated",
                                    )
                                    if level != "whole_reaction"
                                    else ()
                                ),
                            )
                        )
                    ),
                )
            )
    recommendations.sort(
        key=lambda item: (
            item.match_level,
            -item.score,
            item.recipe_id,
            item.precedent_reaction_ids,
        )
    )
    warnings = list(base.pop("warnings"))
    if support() < minimum:
        warnings.append("LIMITED_PRECEDENT_SUPPORT")
    if truncated:
        warnings.append("SHARED_CORE_CANDIDATE_BUDGET_REACHED")
    traces.append(
        RetrievalLevelTrace(
            "shared_core:qualified_union",
            len(channels),
            summarize_evidence_support(
                index.select(sorted(channels))
            ).independent_count,
            len(eligible),
            support(),
            len(channels) - len(eligible),
            minimum,
            "selected"
            if support() >= minimum
            else "selected_limited_support"
            if eligible
            else "empty",
        )
    )
    return GenericRecommendationResult(
        **base,
        valid=bool(recommendations),
        error=None if recommendations else "NO_QUALIFIED_SHARED_CORE_PRECEDENT",
        retrieval_level="shared_core:qualified_union",
        candidate_count=len(channels),
        independent_candidate_count=traces[-1].independent_candidate_count,
        compatible_candidate_count=len(eligible),
        independent_compatible_candidate_count=support(),
        excluded_candidate_count=len(channels) - len(eligible),
        retrieval_trace=tuple(traces),
        shared_core_trace=tuple(diagnostics),
        warnings=tuple(dict.fromkeys(warnings)),
        recommendations=tuple(
            replace(item, rank=i) for i, item in enumerate(recommendations[:top_k], 1)
        ),
    )
