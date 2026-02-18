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
        recs = raw if isinstance(raw, list) else raw.get("recommendations", [])
        cleaned = []
        for i, rec in enumerate(recs[:top_k], 1):
            if isinstance(rec, dict):
                entry = {
                    "rank": i,
                    "catalyst": rec.get("catalyst") or rec.get("metal_catalyst") or "",
                    "ligand": rec.get("ligand") or "",
                    "base": rec.get("base") or "",
                    "solvent": rec.get("solvent") or "",
                    "temperature": rec.get("temperature") or rec.get("temp") or "",
                    "atmosphere": rec.get("atmosphere") or "",
                    "confidence": round(float(rec.get("confidence", 0.0)), 3),
                    "precedent_count": int(rec.get("precedent_count", 0)),
                    "reaction_key": rec.get("reaction_key") or "",
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
