"""
Condition recommendation tools for ChemCoworker.

Wraps the HTE-based condition recommender:
  - recommend_conditions : data-driven condition recommendations for a reaction
"""
from __future__ import annotations

import itertools
import pathlib
from collections import defaultdict
from typing import Any, Dict, List

from ._helpers import _clean_rxn_smiles, _error, _success, _to_jsonable
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Diversity helpers
# ---------------------------------------------------------------------------

_METALS = ["Pd", "Ni", "Cu", "Ir", "Rh", "Ru", "Fe", "Co", "Au", "Zn", "Pt"]


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
        metal = _extract_metal(rec.get("catalyst", ""))
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
) -> str:
    return (
        f"top_k={int(top_k)}"
        f"|source_group={str(source_group or '').strip().lower()}"
        f"|reaction_key_only={1 if reaction_key_only else 0}"
        f"|use_spectator_groups={1 if use_spectator_groups else 0}"
        f"|prefer_mixfp_for_similarity={1 if prefer_mixfp_for_similarity else 0}"
        f"|similarity_mixfp_weight={float(similarity_mixfp_weight):.4f}"
    )


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
        diverse_recs = _diversify(recs, top_k)
        cleaned = []
        for i, rec in enumerate(diverse_recs, 1):
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
                "num_experiments":   int(scores.get("num_experiments", rec.get("num_experiments", 0)) or 0),
                "reaction_type":     rec.get("reaction") or rec.get("reaction_type") or "",
                "reaction_category": rec.get("reaction_category") or "",
                "reactant_types":    rec.get("reactant_types") or [],
                "source":            rec.get("source") or rec.get("source_group") or (source_group or ""),
                "precedent_ids":     rec.get("reaction_id") or rec.get("precedent_ids") or "",
            }
            cleaned.append(entry)

        result = _success({
            "reaction_smiles": reaction_smiles,
            "recommendations": cleaned,
            "total_available": len(recs),
            "source_group": source_group or "all",
            "catalyst_families": sorted({_extract_metal(r["catalyst"]) for r in cleaned if r.get("catalyst")}),
            "detected_reaction_type": (raw.get("detection") or {}).get("detected_reaction_type") if isinstance(raw, dict) else None,
            "hte_timing_ms": hte_timing_ms or {},
            "hte_processing_time_ms": hte_processing_time_ms,
            "hte_recommender_stage_timing_ms": hte_recommender_stage_timing_ms or {},
        })
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


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

CONDITIONS_TOOLS = [
    recommend_conditions_tool,
]
