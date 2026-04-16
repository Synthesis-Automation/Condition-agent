"""
Condition recommendation tools for ChemCoworker.

Wraps the HTE-based condition recommender:
  - recommend_conditions : data-driven condition recommendations for a reaction
"""
from __future__ import annotations

import itertools
import pathlib
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence

from ._helpers import _clean_rxn_smiles, _error, _success, _to_jsonable, _validate_reaction_smiles
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Diversity helpers
# ---------------------------------------------------------------------------

_METALS = ["Pd", "Ni", "Cu", "Ir", "Rh", "Ru", "Fe", "Co", "Au", "Zn", "Pt"]
_CONDITION_SOURCE_ORDER = ("literature", "motif", "similarity", "rules")
_CONDITION_SOURCE_WEIGHTS: Dict[str, float] = {
    "literature": 1.00,
    "motif": 0.95,
    "similarity": 0.75,
    "rules": 0.45,
}


def _extract_metal(catalyst_str: str) -> str:
    """Return the primary metal symbol from a catalyst name, or 'Other'."""
    s = str(catalyst_str)
    for m in _METALS:
        if m in s:
            return m
    return "Other"


def _diversify(recs: List[Dict], top_k: int) -> List[Dict]:
    """Round-robin across catalyst-metal families to return diverse top_k results.

    Preserves best (lowest-rank / highest-Z) entry per metal family first,
    then fills remaining slots in rank order across families.
    """
    buckets: dict = defaultdict(list)
    for rec in recs:
        cond = rec.get("conditions") if isinstance(rec, dict) else {}
        catalyst = rec.get("catalyst", "") if isinstance(rec, dict) else ""
        if isinstance(cond, dict) and not catalyst:
            catalyst = cond.get("catalyst", "")
        metal = _extract_metal(catalyst)
        buckets[metal].append(rec)

    diverse: List[Dict] = []
    for rec in itertools.chain.from_iterable(
        itertools.zip_longest(*buckets.values())
    ):
        if rec is not None and len(diverse) < top_k:
            diverse.append(rec)
    return diverse


# Absolute path to HTE_db — resolved from this file's location so it works
# regardless of the current working directory when the CLI is invoked.
_HTE_DB_PATH = str(
    pathlib.Path(__file__).parent.parent.parent / "data" / "HTE_db"
)


def _conditions_cache_key(
    *,
    top_k: int,
    source_group: str = "",
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
    selection_mode: str = "best",
) -> str:
    return (
        f"top_k={int(top_k)}"
        f"|source_group={str(source_group or '').strip().lower()}"
        f"|reaction_key_only={1 if reaction_key_only else 0}"
        f"|use_spectator_groups={1 if use_spectator_groups else 0}"
        f"|prefer_mixfp_for_similarity={1 if prefer_mixfp_for_similarity else 0}"
        f"|similarity_mixfp_weight={float(similarity_mixfp_weight):.4f}"
        f"|selection_mode={str(selection_mode or 'best').strip().lower()}"
    )


def _condition_context_cache_key() -> str:
    return "condition_context:v1"


def _condition_evidence_cache_key(
    *,
    source: str,
    top_k: int,
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
    selection_mode: str = "best",
) -> str:
    suffix = _conditions_cache_key(
        top_k=int(top_k),
        source_group=source,
        reaction_key_only=reaction_key_only,
        use_spectator_groups=use_spectator_groups,
        prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
        similarity_mixfp_weight=similarity_mixfp_weight,
        selection_mode=selection_mode,
    )
    return (
        f"condition_evidence:v1|source={str(source or '').strip().lower()}"
        f"|{suffix}"
    )


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None:
            return default
        return float(value)
    except Exception:
        return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        if value is None:
            return default
        return int(value)
    except Exception:
        return default


def _condition_fingerprint(rec: Dict[str, Any]) -> str:
    """Stable fingerprint for deduplicating condition sets across sources."""
    fields = [
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "secondary_solvent",
        "additive",
        "coupling_reagent",
        "temperature",
        "atmosphere",
    ]
    parts = []
    for key in fields:
        val = str(rec.get(key) or "").strip().lower()
        parts.append(f"{key}={val}")
    return "|".join(parts)


def _condition_support_score(rec: Dict[str, Any]) -> float:
    """Deterministic support score from recommendation quality metrics."""
    confidence = max(0.0, min(1.0, _as_float(rec.get("confidence"), 0.0)))
    success_rate = max(0.0, min(1.0, _as_float(rec.get("success_rate"), 0.0)))
    avg_yield = _as_float(rec.get("avg_yield"), 0.0)
    if avg_yield > 1.0:
        avg_yield = avg_yield / 100.0
    avg_yield = max(0.0, min(1.0, avg_yield))
    median_yield = _as_float(rec.get("median_yield"), 0.0)
    if median_yield > 1.0:
        median_yield = median_yield / 100.0
    median_yield = max(0.0, min(1.0, median_yield))
    n_exp = _as_int(rec.get("num_experiments"), 0)
    n_exp_norm = max(0.0, min(1.0, n_exp / 20.0))
    return round(
        0.35 * confidence
        + 0.20 * success_rate
        + 0.20 * avg_yield
        + 0.10 * median_yield
        + 0.15 * n_exp_norm,
        4,
    )


def _dedupe_warnings(warnings: Sequence[Any]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for warning in warnings:
        text = str(warning or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return out


def _build_condition_reaction_context(reaction_smiles: str) -> Dict[str, Any]:
    """Build reusable taxonomy/feature context for condition recommendation."""
    try:
        reaction_smiles, rxn_err = _validate_reaction_smiles(reaction_smiles, require_product=True)
        if rxn_err:
            return _error(f"Invalid reaction_smiles: {rxn_err}")

        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None

        cache_key = _condition_context_cache_key()
        if _rtx is not None and hasattr(_rtx, "get_cached_condition_context"):
            cached = _rtx.get_cached_condition_context(reaction_smiles, cache_key=cache_key)
            if isinstance(cached, dict):
                return cached

        from chemtools.recommend import RecommendationRequest
        from chemtools.recommend.query_analysis import analyze_recommendation_query

        analysis = analyze_recommendation_query(
            RecommendationRequest(
                reaction_smiles=reaction_smiles,
                analysis_only=True,
            )
        )
        detected_reaction_type = (
            analysis.detected_reaction_type_id
            or analysis.detected_reaction_type
            or analysis.detected_reaction_type_name
            or ""
        )
        result = _success({
            "reaction_smiles": reaction_smiles,
            "normalized_reaction_smiles": analysis.reaction_smiles_normalized or reaction_smiles,
            "reactants": list(analysis.reactants),
            "agents": list(analysis.agents),
            "products": list(analysis.products),
            "reactant_a_smiles": analysis.reactant_a_smiles,
            "reactant_b_smiles": analysis.reactant_b_smiles or "",
            "product_smiles": analysis.product_smiles or "",
            "reaction_type": detected_reaction_type,
            "reaction_type_id": analysis.detected_reaction_type_id or detected_reaction_type,
            "reaction_type_name": analysis.detected_reaction_type_name or detected_reaction_type,
            "reaction_type_category": analysis.detected_reaction_type_category or "",
            "reaction_type_confidence": float(analysis.reaction_type_confidence or 0.0),
            "reaction_key": analysis.reaction_key or "",
            "reacted_motifs": list(analysis.reacted_motifs),
            "formed_motifs": list(analysis.formed_motifs),
            "spectator_motifs": list(analysis.spectator_motifs),
            "spectator_groups": list(analysis.spectator_groups),
            "requested_reaction_type_filter": analysis.requested_reaction_type_filter or "",
            "requested_reaction_type_filter_canonical": analysis.requested_reaction_type_filter_canonical or "",
            "feature_summary": dict(analysis.raw_feature_summary or {}),
            "warnings": list(analysis.warnings),
        })
        if analysis.warnings:
            result["_warnings"] = list(analysis.warnings)
        if _rtx is not None and hasattr(_rtx, "set_cached_condition_context"):
            _rtx.set_cached_condition_context(reaction_smiles, result, cache_key=cache_key)
        return result
    except Exception as exc:
        return _error(f"Condition reaction context build failed: {exc}")


def _project_condition_reaction_context(result: dict) -> Dict[str, Any]:
    """Project condition-context output into structured fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    out = {
        "reaction_type": result.get("reaction_type_id") or result.get("reaction_type"),
        "reaction_family": result.get("reaction_type_name") or result.get("reaction_type"),
    }
    return {k: v for k, v in out.items() if v is not None}


def _condition_slot_fragments(rec: Dict[str, Any]) -> Dict[str, Any]:
    keys = (
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "secondary_solvent",
        "additive",
        "coupling_reagent",
        "temperature",
        "atmosphere",
    )
    return {key: rec.get(key) for key in keys if str(rec.get(key) or "").strip()}


def _condition_evidence_hits_from_result(source: str, payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Normalize source-specific recommendations into a common evidence schema."""
    recommendations = list(payload.get("recommendations", []) or [])
    top_level_warnings = list(payload.get("_warnings", []) or [])
    payload_evidence = payload.get("evidence") or {}
    hits: List[Dict[str, Any]] = []
    for rec in recommendations:
        if not isinstance(rec, dict):
            continue
        rec_warnings = list(top_level_warnings)
        missing_required_fields = [
            str(value).strip()
            for value in (rec.get("missing_required_fields") or [])
            if str(value).strip()
        ]
        if missing_required_fields:
            rec_warnings.append(
                "Missing required condition fields: " + ", ".join(missing_required_fields)
            )
        num_experiments = _as_int(rec.get("num_experiments"), 0)
        avg_yield = rec.get("avg_yield")
        match_score = rec.get("match_score")
        summary_bits = []
        if num_experiments:
            summary_bits.append(f"{num_experiments} experiments")
        if avg_yield not in (None, ""):
            summary_bits.append(f"avg_yield={avg_yield}")
        if match_score not in (None, ""):
            summary_bits.append(f"match_score={match_score}")
        if not summary_bits:
            summary_bits.append("condition evidence")
        hits.append({
            "source": source,
            "rank": _as_int(rec.get("rank"), len(hits) + 1),
            "score": round(
                _condition_support_score(rec) * _CONDITION_SOURCE_WEIGHTS.get(source, 0.5),
                4,
            ),
            "confidence": _as_float(rec.get("confidence"), 0.0),
            "match_summary": "; ".join(summary_bits),
            "matched_features": {
                "reaction_type": payload.get("detected_reaction_type") or rec.get("reaction_type") or "",
                "reaction_category": rec.get("reaction_category") or "",
                "reactant_types": list(rec.get("reactant_types") or []),
                "matched_transformation": payload_evidence.get("matched_transformation") or "",
            },
            "proposed_fragments": _condition_slot_fragments(rec),
            "warnings": _dedupe_warnings(rec_warnings),
            "provenance": {
                "precedent_ids": rec.get("precedent_ids") or "",
                "num_experiments": num_experiments,
                "success_rate": rec.get("success_rate"),
                "avg_yield": avg_yield,
                "median_yield": rec.get("median_yield"),
                "database_coverage": payload_evidence.get("database_coverage"),
                "reaction_type_confidence": payload.get("reaction_type_confidence"),
            },
        })
    return hits


def _get_condition_evidence(
    reaction_smiles: str,
    *,
    source: str,
    top_k: int = 5,
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    """Internal source-specific evidence adapter for condition recommendation."""
    normalized_source = str(source or "").strip().lower()
    if normalized_source not in _CONDITION_SOURCE_ORDER:
        return _error(
            "source must be one of: literature, motif, similarity, rules"
        )

    reaction_smiles, rxn_err = _validate_reaction_smiles(reaction_smiles, require_product=True)
    if rxn_err:
        return _error(f"Invalid reaction_smiles: {rxn_err}")

    try:
        from chem_coworker.tool_runtime import get_current_tool_runtime_context
        _rtx = get_current_tool_runtime_context()
    except Exception:
        _rtx = None

    cache_key = _condition_evidence_cache_key(
        source=normalized_source,
        top_k=int(top_k),
        reaction_key_only=bool(reaction_key_only),
        use_spectator_groups=bool(use_spectator_groups),
        prefer_mixfp_for_similarity=bool(prefer_mixfp_for_similarity),
        similarity_mixfp_weight=float(similarity_mixfp_weight),
        selection_mode=selection_mode,
    )
    if _rtx is not None and hasattr(_rtx, "get_cached_condition_evidence"):
        cached = _rtx.get_cached_condition_evidence(reaction_smiles, cache_key=cache_key)
        if isinstance(cached, dict):
            return cached

    result = _recommend_conditions(
        reaction_smiles=reaction_smiles,
        top_k=top_k,
        source_group=normalized_source,
        reaction_key_only=reaction_key_only,
        use_spectator_groups=use_spectator_groups,
        prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
        similarity_mixfp_weight=similarity_mixfp_weight,
        selection_mode=selection_mode,
    )
    if not isinstance(result, dict) or not result.get("success"):
        return result

    out = _success({
        "reaction_smiles": reaction_smiles,
        "source": normalized_source,
        "recommendations": list(result.get("recommendations", []) or []),
        "condition_evidence": _condition_evidence_hits_from_result(normalized_source, result),
        "detected_reaction_type": result.get("detected_reaction_type") or "",
        "reaction_type_confidence": result.get("reaction_type_confidence"),
        "evidence": dict(result.get("evidence") or {}),
        "source_group": normalized_source,
        "selection_mode": result.get("selection_mode") or selection_mode,
        "total_available": _as_int(result.get("total_available"), 0),
        "hte_timing_ms": dict(result.get("hte_timing_ms") or {}),
        "hte_processing_time_ms": result.get("hte_processing_time_ms"),
        "hte_recommender_stage_timing_ms": dict(result.get("hte_recommender_stage_timing_ms") or {}),
    })
    warnings = _dedupe_warnings(result.get("_warnings", []))
    if warnings:
        out["_warnings"] = warnings
    if _rtx is not None and hasattr(_rtx, "set_cached_condition_evidence"):
        _rtx.set_cached_condition_evidence(reaction_smiles, out, cache_key=cache_key)
    return out


def _get_literature_condition_evidence(
    reaction_smiles: str,
    top_k: int = 5,
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    return _get_condition_evidence(
        reaction_smiles,
        source="literature",
        top_k=top_k,
        reaction_key_only=reaction_key_only,
        use_spectator_groups=use_spectator_groups,
        selection_mode=selection_mode,
    )


def _get_motif_condition_evidence(
    reaction_smiles: str,
    top_k: int = 5,
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    return _get_condition_evidence(
        reaction_smiles,
        source="motif",
        top_k=top_k,
        reaction_key_only=reaction_key_only,
        use_spectator_groups=use_spectator_groups,
        selection_mode=selection_mode,
    )


def _get_similarity_condition_evidence(
    reaction_smiles: str,
    top_k: int = 5,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    return _get_condition_evidence(
        reaction_smiles,
        source="similarity",
        top_k=top_k,
        reaction_key_only=False,
        use_spectator_groups=use_spectator_groups,
        prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
        similarity_mixfp_weight=similarity_mixfp_weight,
        selection_mode=selection_mode,
    )


def _get_rule_condition_evidence(
    reaction_smiles: str,
    top_k: int = 5,
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    return _get_condition_evidence(
        reaction_smiles,
        source="rules",
        top_k=top_k,
        reaction_key_only=reaction_key_only,
        use_spectator_groups=use_spectator_groups,
        selection_mode=selection_mode,
    )


def _compose_condition_candidates(
    reaction_smiles: str,
    top_k: int = 5,
    sources: Optional[List[str]] = None,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    """Deterministically merge cached/source-specific condition evidence into ranked candidates."""
    reaction_smiles, rxn_err = _validate_reaction_smiles(reaction_smiles, require_product=True)
    if rxn_err:
        return _error(f"Invalid reaction_smiles: {rxn_err}")

    requested_sources = [
        str(source).strip().lower()
        for source in (sources or list(_CONDITION_SOURCE_ORDER))
        if str(source).strip()
    ]
    if not requested_sources:
        requested_sources = list(_CONDITION_SOURCE_ORDER)
    normalized_sources: List[str] = []
    for source in requested_sources:
        if source not in _CONDITION_SOURCE_ORDER:
            return _error("sources must be chosen from: literature, motif, similarity, rules")
        if source not in normalized_sources:
            normalized_sources.append(source)

    per_source: Dict[str, Dict[str, Any]] = {}
    merged: Dict[str, Dict[str, Any]] = {}
    warnings: List[str] = []

    for source in normalized_sources:
        payload = _get_condition_evidence(
            reaction_smiles,
            source=source,
            top_k=top_k,
            selection_mode=selection_mode,
        )
        per_source[source] = payload if isinstance(payload, dict) else _error("Unexpected evidence result")
        if not isinstance(payload, dict):
            continue
        warnings.extend(payload.get("_warnings", []) or [])
        if not payload.get("success"):
            continue
        for rec in list(payload.get("recommendations", []) or []):
            if not isinstance(rec, dict):
                continue
            fp = _condition_fingerprint(rec)
            bucket = merged.get(fp)
            if bucket is None:
                bucket = {
                    "fingerprint": fp,
                    "condition": dict(rec),
                    "sources": [],
                    "source_scores": {},
                    "support_score": 0.0,
                    "ensemble_score": 0.0,
                }
                merged[fp] = bucket

            if source in bucket["sources"]:
                continue

            source_weight = _CONDITION_SOURCE_WEIGHTS.get(source, 0.5)
            support_score = _condition_support_score(rec)
            bucket["sources"].append(source)
            bucket["source_scores"][source] = {
                "source_weight": source_weight,
                "support_score": support_score,
            }
            bucket["support_score"] = round(max(float(bucket["support_score"]), support_score), 4)

            current_best = _condition_support_score(bucket["condition"])
            if support_score > current_best:
                bucket["condition"] = dict(rec)

    ranked_rows: List[Dict[str, Any]] = []
    for row in merged.values():
        weighted_sum = 0.0
        for source in row["sources"]:
            meta = row["source_scores"].get(source) or {}
            weighted_sum += _as_float(meta.get("support_score"), 0.0) * _as_float(meta.get("source_weight"), 0.0)
        consensus_bonus = 0.10 * max(0, len(row["sources"]) - 1)
        row["ensemble_score"] = round(weighted_sum + consensus_bonus, 4)
        ranked_rows.append(row)

    ranked_rows.sort(
        key=lambda row: (
            row.get("ensemble_score", 0.0),
            len(row.get("sources", [])),
            row.get("support_score", 0.0),
        ),
        reverse=True,
    )

    recommendations: List[Dict[str, Any]] = []
    for rank, row in enumerate(ranked_rows[: max(1, min(int(top_k or 5), 10))], 1):
        rec = dict(row["condition"])
        rec["rank"] = rank
        rec["source_consensus"] = list(row["sources"])
        rec["consensus_count"] = len(row["sources"])
        rec["ensemble_score"] = row["ensemble_score"]
        rec["full_support_score"] = row["support_score"]
        rec["balanced_support_score"] = row["support_score"]
        recommendations.append(rec)

    source_counts = {
        source: len((per_source.get(source) or {}).get("recommendations", []) or [])
        for source in normalized_sources
    }
    failed_sources = [
        source for source in normalized_sources
        if not (per_source.get(source) or {}).get("success", False)
    ]

    result = _success({
        "reaction_smiles": reaction_smiles,
        "recommendations": recommendations,
        "sources_run": normalized_sources,
        "source_counts": source_counts,
        "per_source_results": per_source,
        "consensus_summary": {
            "unique_condition_sets": len(ranked_rows),
            "failed_sources": failed_sources,
            "sources_run": normalized_sources,
        },
    })
    deduped_warnings = _dedupe_warnings(warnings)
    if deduped_warnings:
        result["_warnings"] = deduped_warnings
    return result


def _project_compose_condition_candidates(result: dict) -> Dict[str, Any]:
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "conditions": result.get("recommendations", []),
    }


# ---------------------------------------------------------------------------
# Tool: recommend_conditions
# ---------------------------------------------------------------------------

def _recommend_conditions(
    reaction_smiles: str,
    top_k: int = 5,
    source_group: str = "",
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
    selection_mode: str = "best",
) -> Dict[str, Any]:
    """Recommend reaction conditions based on HTE experimental data.

    Uses a high-throughput experimentation (HTE) database to find conditions
    that have worked for reactions with similar motifs and bond changes.
    Ranks by Z-score (avg_yield, success_rate, num_experiments) with direct
    reaction-key matches prioritized before fallback matching.

    Call this after detect_reaction_type to ensure the reaction is properly
    classified before condition matching.

    Args:
        reaction_smiles: Reaction SMILES in 'reactants>>products' format.
        top_k: Number of top recommendations to return (default 5).
        source_group: Filter to a specific data source. One of:
            "literature" — published literature conditions (default pool),
            "motif"      — experimental HTE motif screen results,
            "rules"      — rule-based fallback conditions,
            "similarity" — precedent KNN similarity search (fastest, no HTE pkl load),
            ""           — all sources combined (default).
        reaction_key_only: If True, only return exact reaction-key matches
            (no fallback to motif/similarity). Use for high-confidence reaction types.
        use_spectator_groups: Weight recommendations by spectator group match
            (e.g. steric/electronic groups on arene). Default True.
        prefer_mixfp_for_similarity: Use mixed-fingerprint (MixFP) for precedent
            similarity ranking instead of standard FP. Default False.
        similarity_mixfp_weight: Weight of MixFP component in similarity score
            (0.0–1.0, default 0.3). Only used when prefer_mixfp_for_similarity=True.
        selection_mode: "best" keeps evidence-ranked recommendations as-is.
            "diverse" round-robins by catalyst metal for screening-style output.

    Returns:
        dict with recommendations list, each containing:
          rank, catalyst, ligand, base, solvent, secondary_solvent, additive,
          coupling_reagent, temperature, atmosphere, confidence,
          success_rate, avg_yield, median_yield, match_score, num_experiments,
          reaction_type, source, precedent_ids.
    """
    try:
        try:
            from chem_coworker.tool_runtime import get_current_tool_runtime_context
            _rtx = get_current_tool_runtime_context()
        except Exception:
            _rtx = None
        cache_key = _conditions_cache_key(
            top_k=int(top_k),
            source_group=source_group,
            reaction_key_only=bool(reaction_key_only),
            use_spectator_groups=bool(use_spectator_groups),
            prefer_mixfp_for_similarity=bool(prefer_mixfp_for_similarity),
            similarity_mixfp_weight=float(similarity_mixfp_weight),
            selection_mode=selection_mode,
        )
        if _rtx is not None and hasattr(_rtx, "get_cached_conditions"):
            try:
                cached = _rtx.get_cached_conditions(reaction_smiles, int(top_k), cache_key=cache_key)
            except TypeError:
                cached = _rtx.get_cached_conditions(reaction_smiles, int(top_k))
            if isinstance(cached, dict):
                return cached

        from chemtools.recommend.hte_adapter import recommend_from_reaction

        reaction_smiles = _clean_rxn_smiles(reaction_smiles)
        normalized_source_group = str(source_group or "").strip().lower()
        normalized_selection_mode = str(selection_mode or "best").strip().lower()
        if normalized_selection_mode not in {"best", "diverse"}:
            return _error("selection_mode must be one of: best, diverse")
        raw: Any = None
        if normalized_source_group == "similarity":
            # Align similarity behavior with the recommendation facade used by GUI full-mode.
            try:
                from chemtools.recommend import RecommendationRequest
                from chemtools.recommend.api import recommend as recommend_facade

                req = RecommendationRequest(
                    reaction_smiles=reaction_smiles,
                    strategy="similarity",
                    source_group="any",
                    top_k=max(top_k * 2, 25),
                    min_experiments=1,
                    hte_db_path=_HTE_DB_PATH,
                    use_spectator_groups=use_spectator_groups,
                    prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
                    similarity_mixfp_weight=float(similarity_mixfp_weight),
                )
                run_result = recommend_facade(req)
                rec_obj = getattr(run_result, "recommendation", None)
                facade_recs = [_to_jsonable(item) for item in list(getattr(rec_obj, "recommendations", []) or [])]
                if facade_recs:
                    raw = facade_recs
            except Exception:
                raw = None

        if raw is None:
            raw = recommend_from_reaction(
                reaction_smiles,
                k=max(top_k * 2, 25),
                hte_db_path=_HTE_DB_PATH,
                source_group=source_group or None,
                reaction_key_only=reaction_key_only,
                use_spectator_groups=use_spectator_groups,
                prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
                similarity_mixfp_weight=float(similarity_mixfp_weight),
            )
        if not raw:
            return _success({
                "reaction_smiles": reaction_smiles,
                "recommendations": [],
                "note": "No HTE precedents found for this reaction type.",
            })

        # Normalize raw recommendations into a clean list
        # "recommended_conditions" is the primary key; "recommendations" is a mirror
        recs = raw if isinstance(raw, list) else (
            raw.get("recommended_conditions") or raw.get("recommendations", [])
        )
        hte_timing_ms = None
        hte_processing_time_ms = None
        hte_recommender_stage_timing_ms = None
        if isinstance(raw, dict):
            meta = raw.get("meta") or {}
            extras = raw.get("extras") or {}
            hte_extra = extras.get("hte") or {}
            hte_timing_ms = meta.get("timing_ms") or hte_extra.get("timing_ms")
            hte_processing_time_ms = meta.get("processing_time_ms")
            hte_recommender_stage_timing_ms = hte_extra.get("recommender_stage_timing_ms")
        selected_recs = (
            _diversify(recs, top_k)
            if normalized_selection_mode == "diverse"
            else list(recs[:top_k])
        )
        cleaned = []
        for i, rec in enumerate(selected_recs, 1):
            if not isinstance(rec, dict):
                continue
            # Conditions are nested in rec["conditions"] sub-dict
            cond = rec.get("conditions") or {}
            scores = rec.get("component_scores") or {}
            confidence_val = rec.get("confidence")
            if confidence_val is None:
                confidence_val = rec.get("confidence_score", 0.0)
            try:
                confidence_num = float(confidence_val or 0.0)
            except Exception:
                confidence_num = 0.0
            if confidence_num > 1.0:
                confidence_num = confidence_num / 100.0
            entry = {
                "rank":              rec.get("rank", i),
                "catalyst":          cond.get("catalyst")          or rec.get("catalyst")          or "",
                "ligand":            cond.get("ligand")            or rec.get("ligand")            or "",
                "base":              cond.get("base")              or rec.get("base")              or "",
                "solvent":           cond.get("solvent")           or rec.get("solvent")           or "",
                "secondary_solvent": cond.get("secondary_solvent") or rec.get("secondary_solvent") or "",
                "additive":          cond.get("additive")          or rec.get("additive")          or "",
                "coupling_reagent":  cond.get("coupling_reagent")  or rec.get("coupling_reagent")  or "",
                "temperature":       cond.get("temperature")       or rec.get("temperature")       or "",
                "atmosphere":        cond.get("atmosphere")        or rec.get("atmosphere")        or "",
                "confidence":        round(confidence_num, 2),
                "success_rate":      scores.get("success_rate", rec.get("success_rate")),
                "avg_yield":         scores.get("avg_yield", rec.get("avg_yield")),
                "median_yield":      scores.get("median_yield", rec.get("median_yield")),
                "match_score":       scores.get("match_score", rec.get("match_score")),
                "condition_quality_score": scores.get("condition_quality_score", rec.get("condition_quality_score", 1.0)),
                "num_experiments":   int(scores.get("num_experiments", rec.get("num_experiments", 0)) or 0),
                "reaction_type":     rec.get("reaction") or rec.get("reaction_type") or "",
                "reaction_category": rec.get("reaction_category") or "",
                "reactant_types":    rec.get("reactant_types") or [],
                "missing_required_fields": rec.get("missing_required_fields") or [],
                "source":            rec.get("source") or rec.get("source_group") or (source_group or ""),
                "precedent_ids":     rec.get("reaction_id") or rec.get("precedent_ids") or "",
            }
            cleaned.append(entry)

        detection = raw.get("detection") or {} if isinstance(raw, dict) else {}
        extras = raw.get("extras") or {} if isinstance(raw, dict) else {}
        hte_extra = extras.get("hte") or {} if isinstance(extras, dict) else {}
        matched_motifs = hte_extra.get("matched_motifs") or []
        matched_transformation = ""
        if isinstance(matched_motifs, list) and matched_motifs:
            matched_transformation = str(matched_motifs[0] or "").strip()
        reaction_type_confidence = detection.get("confidence")
        if reaction_type_confidence in (None, ""):
            reaction_type_confidence = hte_extra.get("reaction_type_confidence")
        try:
            reaction_type_confidence = (
                None if reaction_type_confidence in (None, "") else float(reaction_type_confidence)
            )
        except Exception:
            reaction_type_confidence = None
        total_matching_experiments = int(hte_extra.get("total_matching_experiments") or 0)
        database_coverage = hte_extra.get("database_coverage")
        try:
            database_coverage = None if database_coverage is None else float(database_coverage)
        except Exception:
            database_coverage = None
        is_fallback_match = bool(hte_extra.get("is_fallback_match", False))
        catalyst_requirement_enforced = bool(hte_extra.get("catalyst_requirement_enforced", False))
        required_catalyst_family = str(hte_extra.get("required_catalyst_family") or "").strip() or None
        required_catalyst_classes = [
            str(value).strip()
            for value in (hte_extra.get("required_catalyst_classes") or [])
            if str(value).strip()
        ]
        filtered_missing_catalyst_rows = int(hte_extra.get("filtered_missing_catalyst_rows") or 0)
        retained_missing_catalyst_rows = int(hte_extra.get("retained_missing_catalyst_rows") or 0)
        condition_quality_family = str(hte_extra.get("condition_quality_family") or "").strip() or None
        penalized_incomplete_condition_rows = int(hte_extra.get("penalized_incomplete_condition_rows") or 0)
        missing_required_condition_fields = hte_extra.get("missing_required_condition_fields") or {}
        if not isinstance(missing_required_condition_fields, dict):
            missing_required_condition_fields = {}
        plausibility_family = str(hte_extra.get("plausibility_family") or "").strip() or None
        filtered_implausible_catalyst_rows = int(hte_extra.get("filtered_implausible_catalyst_rows") or 0)
        penalized_implausible_mechanism_rows = int(hte_extra.get("penalized_implausible_mechanism_rows") or 0)
        plausibility_issue_counts = hte_extra.get("plausibility_issue_counts") or {}
        if not isinstance(plausibility_issue_counts, dict):
            plausibility_issue_counts = {}

        warnings: List[str] = []
        if reaction_type_confidence is not None and reaction_type_confidence < 0.5:
            warnings.append(
                f"Low reaction-type confidence for condition retrieval ({reaction_type_confidence:.2f}); treat recommendations as tentative."
            )
        if is_fallback_match:
            warnings.append(
                "Condition retrieval used fallback transformation matching rather than a stronger exact structured match."
            )
        if cleaned and int(cleaned[0].get("num_experiments", 0) or 0) < 2:
            warnings.append(
                "Top condition recommendation has limited experimental support (<2 experiments)."
            )
        if retained_missing_catalyst_rows > 0 and catalyst_requirement_enforced:
            catalyst_desc = ", ".join(required_catalyst_classes) if required_catalyst_classes else "a catalyst"
            family_desc = required_catalyst_family or "the detected reaction family"
            warnings.append(
                f"Matched source rows for {family_desc} are missing required catalyst annotations ({catalyst_desc}); source data is incomplete."
            )
        elif filtered_missing_catalyst_rows > 0 and catalyst_requirement_enforced:
            warnings.append(
                f"Filtered {filtered_missing_catalyst_rows} matched rows with missing catalyst annotations for a catalyst-required reaction family."
            )
        if penalized_incomplete_condition_rows > 0:
            field_bits = ", ".join(
                f"{field}={count}" for field, count in sorted(missing_required_condition_fields.items())
            )
            family_desc = condition_quality_family or required_catalyst_family or "the detected reaction family"
            suffix = f" ({field_bits})" if field_bits else ""
            warnings.append(
                f"Penalized {penalized_incomplete_condition_rows} matched rows with incomplete reaction-critical conditions for {family_desc}{suffix}."
            )
        if filtered_implausible_catalyst_rows > 0:
            family_desc = plausibility_family or condition_quality_family or required_catalyst_family or "the detected reaction family"
            warnings.append(
                f"Filtered {filtered_implausible_catalyst_rows} matched rows with catalyst annotations that were chemically implausible for {family_desc}."
            )
        if penalized_implausible_mechanism_rows > 0:
            family_desc = plausibility_family or condition_quality_family or required_catalyst_family or "the detected reaction family"
            warnings.append(
                f"Penalized {penalized_implausible_mechanism_rows} matched rows whose SNAr-style mechanism labels were not plausible for the matched leaving-group pattern in {family_desc}."
            )
        if cleaned and cleaned[0].get("missing_required_fields"):
            family_desc = condition_quality_family or required_catalyst_family or "this reaction family"
            missing_desc = ", ".join(str(value) for value in cleaned[0]["missing_required_fields"] if str(value).strip())
            warnings.append(
                f"Top condition recommendation is still missing required fields for {family_desc}: {missing_desc}."
            )

        result = _success({
            "reaction_smiles": reaction_smiles,
            "recommendations": cleaned,
            "total_available": len(recs),
            "source_group": source_group or "all",
            "selection_mode": normalized_selection_mode,
            "catalyst_families": sorted({_extract_metal(r["catalyst"]) for r in cleaned if r.get("catalyst")}),
            "detected_reaction_type": detection.get("detected_reaction_type") if isinstance(detection, dict) else None,
            "reaction_type_confidence": reaction_type_confidence,
            "evidence": {
                "reaction_type_confidence": reaction_type_confidence,
                "is_fallback_match": is_fallback_match,
                "matched_transformation": matched_transformation,
                "total_matching_experiments": total_matching_experiments,
                "database_coverage": database_coverage,
                "catalyst_requirement_enforced": catalyst_requirement_enforced,
                "required_catalyst_family": required_catalyst_family,
                "required_catalyst_classes": required_catalyst_classes,
                "filtered_missing_catalyst_rows": filtered_missing_catalyst_rows,
                "retained_missing_catalyst_rows": retained_missing_catalyst_rows,
                "condition_quality_family": condition_quality_family,
                "penalized_incomplete_condition_rows": penalized_incomplete_condition_rows,
                "missing_required_condition_fields": missing_required_condition_fields,
                "plausibility_family": plausibility_family,
                "filtered_implausible_catalyst_rows": filtered_implausible_catalyst_rows,
                "penalized_implausible_mechanism_rows": penalized_implausible_mechanism_rows,
                "plausibility_issue_counts": plausibility_issue_counts,
            },
            "hte_timing_ms": hte_timing_ms or {},
            "hte_processing_time_ms": hte_processing_time_ms,
            "hte_recommender_stage_timing_ms": hte_recommender_stage_timing_ms or {},
        })
        if warnings:
            result["_warnings"] = warnings
        if _rtx is not None and hasattr(_rtx, "set_cached_conditions"):
            try:
                _rtx.set_cached_conditions(reaction_smiles, int(top_k), result, cache_key=cache_key)
            except TypeError:
                _rtx.set_cached_conditions(reaction_smiles, int(top_k), result)
        return result
    except Exception as exc:
        return _error(f"Condition recommendation failed: {exc}")


# ---------------------------------------------------------------------------
# Validator for recommend_conditions (Phase 1)
# Moves _check_hypothesis logic onto the tool that owns it.
# ---------------------------------------------------------------------------

def _validate_recommend_conditions(result: dict) -> object:
    """Post-execution validator: surfaces no-HTE-data warnings at tool time."""
    if not result.get("success"):
        return None
    recs = result.get("recommendations", [])
    if not recs:
        return (
            "recommend_conditions returned NO HTE precedents. "
            "State clearly that no experimental data was found. "
            "Do NOT invent conditions."
        )
    top = recs[0]
    n_exp = int(top.get("num_experiments", 0))
    conf = float(top.get("confidence", 1.0))
    if n_exp == 0 and conf < 0.3:
        return (
            "Top HTE recommendation has no experimental support "
            f"(0 experiments, confidence={conf:.2f}). "
            "Treat suggested conditions as tentative and state this explicitly."
        )
    return None


def _project_recommend_conditions(result: dict) -> Dict[str, Any]:
    """Project condition recommendations into structured output fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "conditions": result.get("recommendations", []),
    }


recommend_conditions_tool = ToolPlugin(
    name="recommend_conditions",
    category="conditions",
    description=(
        "Recommend reaction conditions (catalyst, ligand, base, solvent, temperature) from HTE experimental data. "
        "Ranks results by Z-score using success_rate, avg_yield, and num_experiments, with direct reaction-key "
        "matches prioritized before fallback matching. "
        "Optional source_group: 'literature' (published conditions), 'motif' (HTE motif screen results), "
        "'rules' (rule-based fallback), 'similarity' (fast KNN precedent search, no HTE pkl load), "
        "or '' for all sources combined. "
        "Use reaction_key_only=True for exact matches only. "
        "Output includes secondary_solvent, coupling_reagent, match_score, and median_yield."
    ),
    prerequisites=["detect_reaction_type"],
    fn=_recommend_conditions,
    provides=["conditions", "recommendations"],
    requires=["reaction_type"],
    validators=[_validate_recommend_conditions],
)

recommend_conditions_tool.structured_projection = _project_recommend_conditions


build_reaction_context_tool = ToolPlugin(
    name="build_reaction_context",
    category="conditions",
    description=(
        "Build reusable reaction context for condition recommendation: normalized reaction SMILES, "
        "electrophile/nucleophile assignment, taxonomy reaction type guess, reaction key, reacted/formed motifs, "
        "and spectator groups. Use this before source-specific condition evidence tools."
    ),
    prerequisites=[],
    fn=_build_condition_reaction_context,
    provides=["reaction_type", "reaction_type_id", "reaction_key", "reacted_motifs", "formed_motifs"],
)
build_reaction_context_tool.structured_projection = _project_condition_reaction_context


get_literature_condition_evidence_tool = ToolPlugin(
    name="get_literature_condition_evidence",
    category="conditions",
    description=(
        "Get literature-backed condition evidence only. Returns source-specific condition recommendations plus a "
        "normalized condition_evidence list with provenance, matched features, and warnings."
    ),
    prerequisites=["build_reaction_context"],
    fn=_get_literature_condition_evidence,
    provides=["condition_evidence", "recommendations"],
    requires=["reaction_type"],
)


get_motif_condition_evidence_tool = ToolPlugin(
    name="get_motif_condition_evidence",
    category="conditions",
    description=(
        "Get motif-screen/HTE condition evidence only. Returns source-specific recommendations plus normalized "
        "condition_evidence for substrate-class and motif-level matches."
    ),
    prerequisites=["build_reaction_context"],
    fn=_get_motif_condition_evidence,
    provides=["condition_evidence", "recommendations"],
    requires=["reaction_type"],
)


get_similarity_condition_evidence_tool = ToolPlugin(
    name="get_similarity_condition_evidence",
    category="conditions",
    description=(
        "Get similarity/predecessent-based condition evidence only. Uses structurally similar precedents and returns "
        "normalized condition_evidence with provenance and match summaries."
    ),
    prerequisites=["build_reaction_context"],
    fn=_get_similarity_condition_evidence,
    provides=["condition_evidence", "recommendations"],
    requires=["reaction_type"],
)


get_rule_condition_evidence_tool = ToolPlugin(
    name="get_rule_condition_evidence",
    category="conditions",
    description=(
        "Get rule-based condition evidence only. Use as a fast deterministic fallback when motif or literature data is sparse."
    ),
    prerequisites=["build_reaction_context"],
    fn=_get_rule_condition_evidence,
    provides=["condition_evidence", "recommendations"],
    requires=["reaction_type"],
)


compose_condition_candidates_tool = ToolPlugin(
    name="compose_condition_candidates",
    category="conditions",
    description=(
        "Merge cached/source-specific condition_evidence from literature, motif, similarity, and rules into a "
        "consensus-ranked recommendation list. Run after one or more source-specific condition evidence tools."
    ),
    prerequisites=[
        "get_literature_condition_evidence",
        "get_motif_condition_evidence",
        "get_similarity_condition_evidence",
        "get_rule_condition_evidence",
    ],
    fn=_compose_condition_candidates,
    provides=["conditions", "recommendations"],
    requires=["condition_evidence"],
)
compose_condition_candidates_tool.structured_projection = _project_compose_condition_candidates


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

CONDITIONS_TOOLS = [
    build_reaction_context_tool,
    get_literature_condition_evidence_tool,
    get_motif_condition_evidence_tool,
    get_similarity_condition_evidence_tool,
    get_rule_condition_evidence_tool,
    compose_condition_candidates_tool,
    recommend_conditions_tool,
]
