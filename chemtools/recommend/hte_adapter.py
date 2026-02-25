"""
HTE-backed recommendation adapter.

Provides recommend_from_reaction / recommend_conditions_structured
using the chemtools.recommend recommender as the primary backend.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
import threading
import time

from .utils import pick_electrophile_nucleophile

_DEFAULT_DB_PATH: Optional[str] = None
_DEFAULT_RECOMMENDER = None
_DEFAULT_RECOMMENDER_LOCK = threading.Lock()


def _get_default_recommender(hte_db_path: Optional[str] = None):
    global _DEFAULT_RECOMMENDER, _DEFAULT_DB_PATH
    path = hte_db_path or "data/HTE_db"
    if _DEFAULT_RECOMMENDER is None or _DEFAULT_DB_PATH != path:
        with _DEFAULT_RECOMMENDER_LOCK:
            if _DEFAULT_RECOMMENDER is None or _DEFAULT_DB_PATH != path:
                from .recommender import HTERecommender
                _DEFAULT_RECOMMENDER = HTERecommender(path)
                _DEFAULT_DB_PATH = path
    return _DEFAULT_RECOMMENDER


def _ms(start: float, end: float) -> float:
    """Convert perf_counter interval to milliseconds with stable rounding."""
    return round((end - start) * 1000, 2)


def _normalize_smiles_list(items: List[Dict[str, Any]]) -> List[str]:
    out: List[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        smi = item.get("smiles_norm") or item.get("largest_smiles") or item.get("input")
        if smi:
            out.append(str(smi))
    return out


def _extract_reaction_parts(reaction: str) -> Tuple[List[str], List[str], str]:
    if not reaction:
        return [], [], ""
    try:
        from chemtools.smiles import normalize_reaction
        norm = normalize_reaction(reaction)
    except Exception:
        return [], [], reaction

    reactants = _normalize_smiles_list(norm.get("reactants") or [])
    products = _normalize_smiles_list(norm.get("products") or [])
    normalized = norm.get("normalized") or reaction
    return reactants, products, normalized


def _select_reactants(
    reactants: List[str],
    reactant_a_smiles: Optional[str],
    reactant_b_smiles: Optional[str],
) -> Tuple[str, Optional[str]]:
    reactant_a = reactant_a_smiles or ""
    reactant_b = reactant_b_smiles or None

    if not reactant_a and reactants:
        if len(reactants) >= 2:
            reactant_a, reactant_b = pick_electrophile_nucleophile(reactants)
        else:
            reactant_a = reactants[0]

    if reactant_a and reactant_b is None and reactants:
        for candidate in reactants:
            if candidate != reactant_a:
                reactant_b = candidate
                break

    return reactant_a, reactant_b


def _rec_key(rec: Any) -> Tuple[str, str, str, str, str, str, str]:
    return (
        rec.catalyst or "",
        rec.ligand or "",
        rec.base or "",
        rec.solvent or "",
        rec.secondary_solvent or "",
        rec.additive or "",
        rec.coupling_reagent or "",
    )


def _serialize_recommendation(rec: Any) -> Dict[str, Any]:
    return {
        "catalyst": rec.catalyst,
        "ligand": rec.ligand,
        "base": rec.base,
        "solvent": rec.solvent,
        "secondary_solvent": rec.secondary_solvent,
        "additive": rec.additive,
        "coupling_reagent": rec.coupling_reagent,
        "spectator_groups": rec.spectator_groups,
        "spectator_score": rec.spectator_score,
        "success_rate": rec.success_rate,
        "avg_yield": rec.avg_yield,
        "median_yield": rec.median_yield,
        "num_experiments": rec.num_experiments,
        "avg_z_score": rec.avg_z_score,
        "confidence_score": rec.confidence_score,
        "match_score": rec.match_score,
        "reaction_type": rec.reaction_type,
        "reaction_category": rec.reaction_category,
        "reaction_id": rec.reaction_id,
        "reaction_key": rec.reaction_key,
        "reaction_events": rec.reaction_events,
        "reactant_types": list(rec.reactant_types) if rec.reactant_types else [],
        "z_score_range": list(rec.z_score_range) if rec.z_score_range else [],
    }


def _build_source_lookup(result: Any) -> Dict[Tuple[str, str, str, str, str, str, str], str]:
    source_lookup: Dict[Tuple[str, str, str, str, str, str, str], str] = {}
    by_source = getattr(result, "recommendations_by_source", {}) or {}
    if not isinstance(by_source, dict):
        return source_lookup
    for source, recs in by_source.items():
        if not recs:
            continue
        for rec in recs:
            source_lookup.setdefault(_rec_key(rec), str(source))
    return source_lookup


def _build_recommendation_entries(result: Any, requested_type: Optional[str]) -> List[Dict[str, Any]]:
    source_lookup = _build_source_lookup(result)
    output: List[Dict[str, Any]] = []

    for rank, rec in enumerate(getattr(result, "recommendations", []) or [], start=1):
        conditions = {
            "catalyst": rec.catalyst,
            "ligand": rec.ligand,
            "base": rec.base,
            "solvent": rec.solvent,
            "secondary_solvent": rec.secondary_solvent,
            "additive": rec.additive,
            "coupling_reagent": rec.coupling_reagent,
        }
        if rec.spectator_groups:
            conditions["spectator_groups"] = rec.spectator_groups

        conditions = {
            key: value for key, value in conditions.items()
            if value not in (None, "", [])
        }

        component_scores = {
            "avg_z_score": rec.avg_z_score,
            "match_score": rec.match_score,
            "success_rate": rec.success_rate,
            "avg_yield": rec.avg_yield,
            "median_yield": rec.median_yield,
            "num_experiments": rec.num_experiments,
            "spectator_score": rec.spectator_score,
        }

        entry: Dict[str, Any] = {
            "rank": rank,
            "conditions": conditions,
            "confidence": rec.confidence_score,
            "component_scores": component_scores,
            "source": source_lookup.get(_rec_key(rec)) or "hte",
            "reaction": rec.reaction_type or requested_type,
            "reaction_category": rec.reaction_category,
            "reaction_id": rec.reaction_id,
            "reactant_types": list(rec.reactant_types) if rec.reactant_types else [],
            "z_score_range": list(rec.z_score_range) if rec.z_score_range else [],
        }

        output.append(entry)

    return output


def recommend_from_reaction(
    reaction: str,
    k: int = 5,
    *_args: Any,
    reaction_type: Optional[str] = None,
    reactant_a_smiles: Optional[str] = None,
    reactant_b_smiles: Optional[str] = None,
    product_smiles: Optional[str] = None,
    hte_db_path: Optional[str] = None,
    min_experiments: int = 2,
    catalyst_filter: Optional[str] = None,
    source_group: Optional[str] = None,
    reaction_key_only: bool = False,
    use_aryl_steric_electronic_weighting: bool = False,
    use_spectator_groups: bool = True,
    **_: Any,
) -> Dict[str, Any]:
    t0 = time.perf_counter()

    t_parse_0 = time.perf_counter()
    reactants, products, normalized = _extract_reaction_parts(reaction)
    reactant_a, reactant_b = _select_reactants(
        reactants, reactant_a_smiles, reactant_b_smiles
    )
    t_parse_1 = time.perf_counter()
    if product_smiles is None and products:
        product_smiles = products[0]

    if not reactant_a:
        t_end = time.perf_counter()
        timing_ms = {
            "input_parse_ms": _ms(t_parse_0, t_parse_1),
            "recommender_get_ms": 0.0,
            "recommend_compute_ms": 0.0,
            "postprocess_ms": 0.0,
            "total_ms": _ms(t0, t_end),
        }
        return {
            "meta": {
                "model": "hte_recommender",
                "status": "error",
                "processing_time_ms": timing_ms["total_ms"],
                "timing_ms": timing_ms,
            },
            "input": {
                "reaction_smiles": normalized or reaction,
                "requested_reaction_type": reaction_type,
                "reactant_a_smiles": reactant_a_smiles,
                "reactant_b_smiles": reactant_b_smiles,
                "product_smiles": product_smiles,
            },
            "detection": {
                "detected_reaction_type": reaction_type,
                "confidence": None,
                "method": "hte",
            },
            "recommended_conditions": [],
            "recommendations": [],
            "error": "No reactants found in reaction SMILES.",
        }

    t_get_0 = time.perf_counter()
    recommender = _get_default_recommender(hte_db_path)
    t_get_1 = time.perf_counter()

    t_rec_0 = time.perf_counter()
    result = recommender.recommend(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        product_smiles=product_smiles,
        top_k=max(1, int(k or 1)),
        min_experiments=min_experiments,
        reaction_type_filter=reaction_type,
        catalyst_filter=catalyst_filter,
        source_group=source_group,
        reaction_key_only=reaction_key_only,
        use_aryl_steric_electronic_weighting=use_aryl_steric_electronic_weighting,
        use_spectator_groups=use_spectator_groups,
    )
    t_rec_1 = time.perf_counter()

    t_post_0 = time.perf_counter()
    recommendations = _build_recommendation_entries(result, reaction_type)
    t_post_1 = time.perf_counter()
    t_end = time.perf_counter()
    timing_ms = {
        "input_parse_ms": _ms(t_parse_0, t_parse_1),
        "recommender_get_ms": _ms(t_get_0, t_get_1),
        "recommend_compute_ms": _ms(t_rec_0, t_rec_1),
        "postprocess_ms": _ms(t_post_0, t_post_1),
        "total_ms": _ms(t0, t_end),
    }
    processing_time_ms = timing_ms["total_ms"]

    detected_type = getattr(result, "predicted_reaction_type", None) or reaction_type
    confidence = getattr(result, "reaction_type_confidence", 0.0) or None
    recommender_stage_timing_ms = getattr(result, "timing_ms", None) or {}

    extras = {
        "hte": {
            "reactant_a_type": getattr(result, "reactant_a_type", None),
            "reactant_b_type": getattr(result, "reactant_b_type", None),
            "reactant_a_category": getattr(result, "reactant_a_category", None),
            "reactant_b_category": getattr(result, "reactant_b_category", None),
            "predicted_reaction_type": getattr(result, "predicted_reaction_type", None),
            "reaction_type_confidence": getattr(result, "reaction_type_confidence", 0.0),
            "total_matching_experiments": getattr(result, "total_matching_experiments", 0),
            "database_coverage": getattr(result, "database_coverage", 0.0),
            "is_fallback_match": getattr(result, "is_fallback_match", False),
            "matched_motifs": list(getattr(result, "matched_motifs", ()) or ()),
            "reacted_motifs": list(getattr(result, "reacted_motifs", ()) or ()),
            "formed_motifs": list(getattr(result, "formed_motifs", ()) or ()),
            "spectator_motifs": list(getattr(result, "spectator_motifs", ()) or ()),
            "recommendations_by_source": {
                str(source): [_serialize_recommendation(rec) for rec in recs]
                for source, recs in (getattr(result, "recommendations_by_source", {}) or {}).items()
            },
            "timing_ms": timing_ms,
            "recommender_stage_timing_ms": recommender_stage_timing_ms,
        }
    }

    return {
        "meta": {
            "model": "hte_recommender",
            "status": "success",
            "processing_time_ms": processing_time_ms,
            "timing_ms": timing_ms,
        },
        "input": {
            "reaction_smiles": normalized or reaction,
            "requested_reaction_type": reaction_type,
            "reactant_a_smiles": reactant_a,
            "reactant_b_smiles": reactant_b,
            "product_smiles": product_smiles,
        },
        "detection": {
            "detected_reaction_type": detected_type,
            "confidence": confidence,
            "method": "hte",
        },
        "recommended_conditions": recommendations,
        "recommendations": recommendations,
        "extras": extras,
    }


def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 5,
    limit: int = 5,
    **kwargs: Any,
) -> Dict[str, Any]:
    results = recommend_from_reaction(
        reaction,
        k=k,
        reaction_type=reaction_type,
        **kwargs,
    )
    recommendations = results.get("recommended_conditions") or []
    results["recommended_conditions"] = recommendations[: max(1, int(limit or 1))]
    results["recommendations"] = results["recommended_conditions"]
    return results
