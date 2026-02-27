"""High-level recommendation facade with plan-first execution and lazy loading."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Tuple

from .data_manager import RecommendationDataManager
from .models import (
    OutputView,
    QueryAnalysis,
    RecommendationRequest,
    RecommendationRunResult,
    RecommendationStrategy,
    SourceGroup,
)
from .planner import plan_sources
from .query_analysis import analyze_recommendation_query


def _normalize_source_group_label(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"", "any", "all"}:
        return "all"
    if text in {"literature", "dataset", "datasets", "lit"}:
        return "literature"
    if text in {"motif", "motifs", "experiments", "experiment", "experiements"}:
        return "motif"
    if text in {"protocols", "protocol"}:
        return "literature"
    if text == "rules":
        return "rules"
    if text == "precedent":
        return "precedent"
    return text


def _coerce_request(request: RecommendationRequest | str, **kwargs: Any) -> RecommendationRequest:
    if isinstance(request, RecommendationRequest):
        if kwargs:
            # Allow minor overrides without mutating caller-owned request.
            data = dict(request.__dict__)
            data.update(kwargs)
            return RecommendationRequest(**data)
        return request
    return RecommendationRequest(reaction_smiles=str(request or ""), **kwargs)


def _normalize_recommendations_by_source(source_map: Dict[str, Any]) -> Dict[str, List[Any]]:
    normalized_map: Dict[str, List[Any]] = {}
    for key, items in (source_map or {}).items():
        normalized = _normalize_source_group_label(key)
        normalized_map.setdefault(normalized, []).extend(list(items or []))
    return normalized_map


def _apply_precedent_only_view(result: Any) -> Any:
    source_map = _normalize_recommendations_by_source(getattr(result, "recommendations_by_source", {}) or {})
    precedent_recs = list(source_map.get("precedent") or [])
    result.recommendations_by_source = {"precedent": precedent_recs}
    result.recommendations = precedent_recs
    return result


def _interleave_recommendation_lists(
    by_source: Dict[str, List[Any]],
    *,
    source_order: Iterable[str],
    limit: int,
) -> List[Any]:
    order = [str(src) for src in source_order if str(src)]
    if not order:
        order = [str(k) for k in by_source.keys()]
    indices = {src: 0 for src in order}
    combined: List[Any] = []
    while limit <= 0 or len(combined) < limit:
        progressed = False
        for src in order:
            items = by_source.get(src) or []
            idx = indices.get(src, 0)
            if idx >= len(items):
                continue
            combined.append(items[idx])
            indices[src] = idx + 1
            progressed = True
            if limit > 0 and len(combined) >= limit:
                break
        if not progressed:
            break
    return combined


def _apply_source_subset_view(
    result: Any,
    *,
    allowed_sources: Iterable[str],
    rename_map: Optional[Dict[str, str]] = None,
    top_k: int = 10,
) -> Any:
    source_map = _normalize_recommendations_by_source(getattr(result, "recommendations_by_source", {}) or {})
    rename_map = dict(rename_map or {})
    filtered: Dict[str, List[Any]] = {}
    ordered_keys: List[str] = []
    for raw_key in allowed_sources:
        key = _normalize_source_group_label(raw_key)
        if key not in source_map:
            continue
        out_key = rename_map.get(key, key)
        filtered[out_key] = list(source_map.get(key) or [])
        ordered_keys.append(out_key)
    result.recommendations_by_source = filtered
    result.recommendations = _interleave_recommendation_lists(
        filtered,
        source_order=ordered_keys,
        limit=top_k,
    )
    return result


def _apply_strategy_view(result: Any, req: RecommendationRequest, plan_strategy: Optional[str]) -> Any:
    strategy = str(plan_strategy or "").strip().lower()
    if not strategy:
        return result
    if strategy == RecommendationStrategy.SIMILARITY.value:
        return _apply_source_subset_view(
            result,
            allowed_sources=("precedent",),
            rename_map={"precedent": "similarity"},
            top_k=req.top_k,
        )
    if strategy == RecommendationStrategy.LITERATURE.value:
        return _apply_source_subset_view(
            result,
            allowed_sources=(SourceGroup.LITERATURE.value,),
            top_k=req.top_k,
        )
    if strategy == RecommendationStrategy.RULES.value:
        return _apply_source_subset_view(
            result,
            allowed_sources=(SourceGroup.RULES.value,),
            top_k=req.top_k,
        )
    if strategy == RecommendationStrategy.MOTIF.value:
        return _apply_source_subset_view(
            result,
            allowed_sources=(SourceGroup.MOTIF.value,),
            top_k=req.top_k,
        )
    return result


def _run_single_pass(
    req: RecommendationRequest,
    dm: RecommendationDataManager,
    analysis: QueryAnalysis,
    *,
    source_group: Optional[str],
    force_precedent_search: bool = False,
) -> Tuple[Any, Dict[str, Any]]:
    recommender, info = dm.get_recommender(source_group=source_group)
    result = recommender.recommend(
        reactant_a_smiles=analysis.reactant_a_smiles,
        reactant_b_smiles=analysis.reactant_b_smiles,
        product_smiles=analysis.product_smiles,
        top_k=req.top_k,
        min_experiments=req.min_experiments,
        reaction_type_filter=req.reaction_type_filter,
        catalyst_filter=req.catalyst_filter,
        source_group=(None if source_group in {None, "", "all"} else source_group),
        reaction_key_only=req.reaction_key_only,
        use_aryl_steric_electronic_weighting=req.use_aryl_steric_electronic_weighting,
        use_spectator_groups=req.use_spectator_groups,
        prefer_mixfp_for_similarity=req.prefer_mixfp_for_similarity,
        similarity_mixfp_weight=req.similarity_mixfp_weight,
        force_precedent_search=force_precedent_search,
    )
    return result, {
        "hte_recommender": {
            "db_path": info.db_path,
            "source_group": info.source_group,
            "cache_hit": info.cache_hit,
        }
    }


def _run_per_source(
    req: RecommendationRequest,
    dm: RecommendationDataManager,
    analysis: QueryAnalysis,
    sources: Iterable[str],
) -> Tuple[Any, Dict[str, Any]]:
    baseline = None
    loaded: Dict[str, Any] = {"runs": []}
    merged_by_source: Dict[str, List[Any]] = {}
    source_order = [str(source) for source in sources]
    for source in source_order:
        recommender, info = dm.get_recommender(source_group=source)
        loaded["runs"].append(
            {
                "source_group": source,
                "db_path": info.db_path,
                "cache_hit": info.cache_hit,
            }
        )
        result = recommender.recommend(
            reactant_a_smiles=analysis.reactant_a_smiles,
            reactant_b_smiles=analysis.reactant_b_smiles,
            product_smiles=analysis.product_smiles,
            top_k=req.top_k,
            min_experiments=req.min_experiments,
            reaction_type_filter=req.reaction_type_filter,
            catalyst_filter=req.catalyst_filter,
            source_group=source,
            reaction_key_only=req.reaction_key_only,
            use_aryl_steric_electronic_weighting=req.use_aryl_steric_electronic_weighting,
            use_spectator_groups=req.use_spectator_groups,
            prefer_mixfp_for_similarity=req.prefer_mixfp_for_similarity,
            similarity_mixfp_weight=req.similarity_mixfp_weight,
        )
        if baseline is None:
            baseline = result
        source_map = _normalize_recommendations_by_source(getattr(result, "recommendations_by_source", {}) or {})
        group_recs = list(source_map.get(source) or [])
        if not group_recs:
            group_recs = list(getattr(result, "recommendations", []) or [])
        merged_by_source[source] = list(group_recs)
        if source == SourceGroup.LITERATURE.value:
            precedent = list(source_map.get("precedent") or [])
            if precedent:
                merged_by_source["precedent"] = precedent

    if baseline is None:
        return None, loaded
    baseline.recommendations_by_source = merged_by_source
    if req.normalized_output_view() is OutputView.BY_SOURCE:
        baseline.recommendations = []
    else:
        baseline.recommendations = _interleave_recommendation_lists(
            merged_by_source,
            source_order=source_order,
            limit=req.top_k,
        )
    return baseline, loaded


def _run_similarity_fast_pass(
    req: RecommendationRequest,
    analysis: QueryAnalysis,
) -> Tuple[Any, Dict[str, Any]]:
    """SIMILARITY fast path — skips loading HTE pkl files entirely.

    Calls the precedent KNN search directly (uses the lightweight disk-cached
    featurized CSV rows, ~5-10s) without touching the 966MB HTE literature pkl.
    """
    from .recommender import HTERecommendationResult, _run_precedent_knn

    # Similarity mode should be cross-family by default.
    # Only honor reaction type when the caller explicitly provides a filter.
    explicit_reaction_type_filter = (
        str(req.reaction_type_filter).strip() if req.reaction_type_filter is not None else ""
    ) or None

    recs = _run_precedent_knn(
        analysis.reactant_a_smiles,
        analysis.reactant_b_smiles,
        analysis.product_smiles,
        explicit_reaction_type_filter,
        req.top_k,
        source_group=None,
        prefer_mixfp_for_similarity=req.prefer_mixfp_for_similarity,
        similarity_mixfp_weight=req.similarity_mixfp_weight,
    )

    rec_obj = HTERecommendationResult(
        reactant_a_smiles=analysis.reactant_a_smiles,
        reactant_b_smiles=analysis.reactant_b_smiles,
        product_smiles=analysis.product_smiles,
        recommendations=list(recs),
        recommendations_by_source={"precedent": list(recs)},
    )
    loaded = {
        "similarity_fast_path": True,
        "precedent_count": len(recs),
        "similarity_family_filter": explicit_reaction_type_filter or "cross_family",
    }
    return rec_obj, loaded


def recommend(
    request: RecommendationRequest | str,
    *,
    data_manager: Optional[RecommendationDataManager] = None,
    **overrides: Any,
) -> RecommendationRunResult:
    """
    High-level recommendation entrypoint.

    This facade performs:
      1. Featurizer-based query analysis (no dataset load)
      2. Source planning (whether data loading is needed)
      3. Lazy resource loading via RecommendationDataManager
      4. Recommendation execution using the legacy HTERecommender backend
    """
    req = _coerce_request(request, **overrides)
    analysis = analyze_recommendation_query(req)
    plan = plan_sources(req)

    run_result = RecommendationRunResult(request=req, analysis=analysis, plan=plan)
    if not plan.needs_hte_data:
        run_result.loaded_resources = {"skipped": "no_dataset_load_required"}
        return run_result

    dm = data_manager or RecommendationDataManager(base_db_path=req.hte_db_path)

    if plan.recommendation_strategy == RecommendationStrategy.SIMILARITY.value:
        rec_obj, loaded = _run_similarity_fast_pass(req, analysis)
    elif plan.run_strategy.value == "per_source":
        rec_obj, loaded = _run_per_source(req, dm, analysis, plan.sources_to_run)
    else:
        source_group = plan.single_run_source_group
        rec_obj, loaded = _run_single_pass(
            req, dm, analysis, source_group=source_group, force_precedent_search=False
        )

    if rec_obj is not None and plan.output_view is OutputView.PRECEDENT_ONLY:
        rec_obj = _apply_precedent_only_view(rec_obj)
    if rec_obj is not None:
        rec_obj = _apply_strategy_view(rec_obj, req, plan.recommendation_strategy)

    run_result.recommendation = rec_obj
    run_result.loaded_resources = loaded
    return run_result


def warm_recommendation_cache(
    *,
    source_group: Optional[str | SourceGroup] = None,
    data_manager: Optional[RecommendationDataManager] = None,
    hte_db_path: str = "data/HTE_db",
    clear_memory_cache: bool = False,
) -> Dict[str, Any]:
    """Explicit preload hook that shares the same DataManager abstraction."""
    dm = data_manager or RecommendationDataManager(base_db_path=hte_db_path)
    return dm.warm(source_group=source_group, clear_memory_cache=clear_memory_cache)
