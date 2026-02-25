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


# ---------------------------------------------------------------------------
# Tool: recommend_conditions
# ---------------------------------------------------------------------------

def _recommend_conditions(reaction_smiles: str, top_k: int = 5) -> Dict[str, Any]:
    """Recommend reaction conditions based on HTE experimental data.

    Uses a high-throughput experimentation (HTE) database to find conditions
    that have worked for reactions with similar motifs and bond changes.
    Ranks by success rate (yield > 50%) and motif similarity.

    Call this after detect_reaction_type to ensure the reaction is properly
    classified before condition matching.

    Args:
        reaction_smiles: Reaction SMILES in 'reactants>>products' format.
        top_k: Number of top recommendations to return (default 5).

    Returns:
        dict with recommendations list, each containing:
          rank, catalyst, ligand, base, solvent, temperature, confidence,
          precedent_count.
    """
    try:
        from chemtools.recommend.hte_adapter import recommend_from_reaction

        reaction_smiles = _clean_rxn_smiles(reaction_smiles)
        raw = recommend_from_reaction(
            reaction_smiles,
            k=max(top_k * 2, 25),
            hte_db_path=_HTE_DB_PATH,
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
        diverse_recs = _diversify(recs, top_k)
        cleaned = []
        for i, rec in enumerate(diverse_recs, 1):
            if not isinstance(rec, dict):
                continue
            # Conditions are nested in rec["conditions"] sub-dict
            cond = rec.get("conditions") or {}
            scores = rec.get("component_scores") or {}
            entry = {
                "rank": rec.get("rank", i),
                "catalyst":    cond.get("catalyst")    or rec.get("catalyst")    or "",
                "ligand":      cond.get("ligand")      or rec.get("ligand")      or "",
                "base":        cond.get("base")        or rec.get("base")        or "",
                "solvent":     cond.get("solvent")     or rec.get("solvent")     or "",
                "additive":    cond.get("additive")    or rec.get("additive")    or "",
                "temperature": cond.get("temperature") or rec.get("temperature") or "",
                "atmosphere":  cond.get("atmosphere")  or rec.get("atmosphere")  or "",
                "confidence":  round(float(rec.get("confidence", 0.0)), 1),
                "success_rate": scores.get("success_rate"),
                "avg_yield":    scores.get("avg_yield"),
                "num_experiments": int(scores.get("num_experiments", 0)),
                "reaction_type": rec.get("reaction") or "",
                "reactant_types": rec.get("reactant_types") or [],
                "source": rec.get("source") or "",
                "precedent_ids": rec.get("reaction_id") or "",
            }
            cleaned.append(entry)

        return _success({
            "reaction_smiles": reaction_smiles,
            "recommendations": cleaned,
            "total_available": len(recs),
            "catalyst_families": sorted({_extract_metal(r["catalyst"]) for r in cleaned if r.get("catalyst")}),
        })
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
    description="Recommend reaction conditions (catalyst, ligand, base, solvent, temperature) from HTE experimental data.",
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
