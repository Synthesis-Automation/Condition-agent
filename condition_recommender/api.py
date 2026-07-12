"""Public chemistry-first condition recommendation API."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

from reactive_taxonomy import featurize_reaction

from .indexing import IndexedReaction, load_verified_index
from .models import ConditionRecommendation, RecommendationResult
from .retrieval import load_retrieval_rules, retrieve_pool, structured_similarity


def _query_features(analysis: Any) -> Dict[str, Any] | None:
    environment = analysis.family_environment
    if environment is None or environment.family_id != "suzuki_miyaura":
        return None
    partners = {partner.role: asdict(partner) for partner in environment.partners}
    if not {"electrophile", "transfer_partner"} <= set(partners):
        return None
    return {
        "electrophile": partners["electrophile"],
        "transfer_partner": partners["transfer_partner"],
        "spectator_group_ids": tuple(sorted({group.group_id for group in analysis.spectator_groups})),
        "family_flags": environment.flags,
    }


def _explanation(query: Dict[str, Any], row: IndexedReaction, components: Dict[str, float]) -> Tuple[str, ...]:
    notes = []
    for name, label in (
        ("electrophile_handle", "electrophile handle"),
        ("electrophile_context", "electrophile context"),
        ("transfer_handle", "transfer handle"),
        ("transfer_context", "transfer-partner context"),
    ):
        if components.get(name) == 1.0:
            notes.append(f"Exact {label} match")
    if components.get("spectators", 0.0) >= 0.5:
        notes.append("Similar spectator functional groups")
    return tuple(notes or ("Suzuki family-level precedent",))


def _aggregate(
    query: Dict[str, Any], pool: Tuple[IndexedReaction, ...], retrieval_level: str, top_k: int,
) -> Tuple[ConditionRecommendation, ...]:
    rules = load_retrieval_rules()
    maximum = int(rules.get("maximum_neighbors", 50))
    scored = []
    for row in pool:
        similarity, components = structured_similarity(query, row)
        scored.append((similarity, row, components))
    scored.sort(key=lambda item: (-item[0], -item[1].yield_pct, item[1].reaction_id))
    groups: Dict[str, List[Tuple[float, IndexedReaction, Dict[str, float]]]] = defaultdict(list)
    for item in scored[:maximum]:
        groups[item[1].recipe_id].append(item)
    ranking = rules["ranking_weights"]
    pool_yield_mean = sum(row.yield_pct for row in pool) / len(pool)
    prior_strength = float(rules.get("yield_prior_strength", 3.0))
    ranked = []
    for recipe_id, members in groups.items():
        weight_sum = sum(max(similarity, 0.05) for similarity, _, _ in members)
        weighted_yield_sum = sum(max(similarity, 0.05) * row.yield_pct for similarity, row, _ in members)
        expected_yield = (
            weighted_yield_sum + prior_strength * pool_yield_mean
        ) / (weight_sum + prior_strength)
        similarity_score = sum(similarity for similarity, _, _ in members) / len(members)
        support_score = min(1.0, math.log1p(len(members)) / math.log1p(10))
        score = (
            float(ranking["similarity"]) * similarity_score
            + float(ranking["yield"]) * expected_yield / 100.0
            + float(ranking["support"]) * support_score
        )
        best = members[0]
        ranked.append((score, expected_yield, similarity_score, recipe_id, members, best))
    ranked.sort(key=lambda item: (-item[0], -item[1], item[3]))
    cautions = tuple(sorted({
        flag.split(":", 1)[-1].replace("_", " ")
        for flag in query.get("family_flags", ())
        if any(token in flag for token in ("risk", "hindered", "challenging", "competing"))
    }))
    recommendations = []
    for rank, (score, expected_yield, similarity_score, recipe_id, members, best) in enumerate(ranked[:top_k], start=1):
        recommendations.append(ConditionRecommendation(
            rank=rank,
            conditions=best[1].conditions,
            recipe_id=recipe_id,
            score=round(score, 6),
            similarity_score=round(similarity_score, 6),
            expected_yield_pct=round(expected_yield, 2),
            support=len(members),
            retrieval_level=retrieval_level,
            precedent_reaction_ids=tuple(item[1].reaction_id for item in members[:5]),
            explanation=_explanation(query, best[1], best[2]),
            cautions=cautions,
        ))
    return tuple(recommendations)


def recommend_conditions(
    reaction_smiles: str,
    *,
    records_path: str | Path = "results/suzuki_recommendation_pilot/verified.csv",
    top_k: int = 5,
) -> RecommendationResult:
    """Recommend condition recipes for one product-specified Suzuki reaction."""
    if top_k < 1:
        return RecommendationResult(reaction_smiles, False, error="TOP_K_MUST_BE_POSITIVE")
    analysis = featurize_reaction(reaction_smiles)
    if not analysis.valid:
        return RecommendationResult(reaction_smiles, False, error=analysis.error or "INVALID_REACTION")
    if analysis.evidence_quality != "exact_product_reconstruction":
        return RecommendationResult(
            reaction_smiles, False, query_label=analysis.reaction_label,
            error="QUERY_REACTION_NOT_EXACTLY_VERIFIED",
        )
    query = _query_features(analysis)
    if query is None:
        return RecommendationResult(
            reaction_smiles, False, query_label=analysis.reaction_label,
            error="QUERY_NOT_VERIFIED_SUZUKI_MIYAURA",
        )
    rows = load_verified_index(records_path)
    if not rows:
        return RecommendationResult(reaction_smiles, False, error="EMPTY_VERIFIED_INDEX")
    retrieval_level, pool = retrieve_pool(query, rows)
    recommendations = _aggregate(query, pool, retrieval_level, top_k)
    return RecommendationResult(
        query_reaction_smiles=reaction_smiles,
        valid=True,
        query_label=analysis.reaction_label,
        product_connection_label=analysis.product_connection.concise_label if analysis.product_connection else None,
        retrieval_level=retrieval_level,
        candidate_count=len(pool),
        recommendations=recommendations,
        warnings=("TEMPERATURE_AND_TIME_NOT_MODELED",),
    )


__all__ = ["recommend_conditions"]
