"""
Condition recommendation tools for ChemCoworker.

Wraps the HTE-based condition recommender:
  - recommend_conditions : data-driven condition recommendations for a reaction
"""
from __future__ import annotations

from typing import Any, Dict

from ._helpers import _error, _success, _to_jsonable
from ._base import ToolPlugin


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

        raw = recommend_from_reaction(reaction_smiles, k=max(top_k * 2, 25))
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
        cleaned = []
        for i, rec in enumerate(recs[:top_k], 1):
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
        })
    except Exception as exc:
        return _error(f"Condition recommendation failed: {exc}")


recommend_conditions_tool = ToolPlugin(
    name="recommend_conditions",
    category="conditions",
    description="Recommend reaction conditions (catalyst, ligand, base, solvent, temperature) from HTE experimental data.",
    prerequisites=["detect_reaction_type"],
    fn=_recommend_conditions,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

CONDITIONS_TOOLS = [
    recommend_conditions_tool,
]
