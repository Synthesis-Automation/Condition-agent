"""
Taxonomy search tools for ChemCoworker.

Wraps the reaction taxonomy and motif databases:
  - search_reaction_types : keyword search across 60+ reaction type definitions
  - search_motifs         : search compound motif labels (Ar-Br, Alkyl-NH2, etc.)

Always run after analyze_bond_changes — bond type constrains which reaction
families are relevant and produces better search results.
"""
from __future__ import annotations

from typing import Any, Dict

from ._helpers import _error, _success
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Tool: search_reaction_types
# ---------------------------------------------------------------------------

def _search_reaction_types(query: str) -> Dict[str, Any]:
    """Search the reaction taxonomy for types matching a keyword query.

    Best practice: include the KEY BOND TYPE and NUCLEOPHILE CLASS in the query.
      Good: "C-N coupling amine aryl halide"
      Bad:  "Pd coupling"  (electrophile only → Suzuki bias)

    Always run analyze_bond_changes before this tool to know the bond type.

    Args:
        query: Search keywords, e.g. 'C-C coupling boronic acid Pd',
               'amide formation', 'SNAr amine electron-poor arene'.

    Returns:
        dict with matches list (id, name, description, score, catalysts).
    """
    try:
        from chem_coworker.taxonomy_prompts import TaxonomyContext
        ctx = TaxonomyContext()
        results = ctx.search_reaction_types(query)
        if not results:
            return _success({
                "query": query,
                "matches": [],
                "note": "No taxonomy matches found. Consider a broader search term.",
            })
        return _success({
            "query": query,
            "matches": results,
            "total_found": len(results),
        })
    except Exception as exc:
        return _error(f"Reaction type search failed: {exc}")


def _project_search_reaction_types(result: dict) -> Dict[str, Any]:
    """Project taxonomy search matches into structured output fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "taxonomy_matches": result.get("matches", []),
    }


search_reaction_types_tool = ToolPlugin(
    name="search_reaction_types",
    category="taxonomy",
    description="Search the reaction taxonomy by keyword (e.g. 'C-N coupling amine aryl halide'). Run AFTER analyze_bond_changes.",
    prerequisites=["analyze_bond_changes"],
    fn=_search_reaction_types,
    llm_exposed=False,
)

search_reaction_types_tool.structured_projection = _project_search_reaction_types


# ---------------------------------------------------------------------------
# Tool: search_motifs
# ---------------------------------------------------------------------------

def _search_motifs(scaffold: str = "", substituent: str = "") -> Dict[str, Any]:
    """Search compound motif labels used for reaction classification.

    Motif labels follow the format 'Scaffold-Substituent', e.g.:
      'Ar-Br', 'Ar-B(OH)2', 'Alkyl-NH2', 'Alkenyl-Cl', 'HetAr-OTf'

    Args:
        scaffold: Scaffold part, e.g. 'Ar', 'Alkyl', 'Alkenyl', 'HetAr'. Empty = all.
        substituent: Substituent part, e.g. 'Br', 'NH2', 'B(OH)2'. Empty = all.

    Returns:
        dict with matches list (motif_id, scaffold, substituent).
    """
    try:
        from chem_coworker.taxonomy_prompts import TaxonomyContext
        ctx = TaxonomyContext()
        results = ctx.search_motifs(scaffold=scaffold, substituent=substituent)
        return _success({
            "scaffold_query": scaffold,
            "substituent_query": substituent,
            "matches": results[:30],
            "total_found": len(results),
        })
    except Exception as exc:
        return _error(f"Motif search failed: {exc}")


search_motifs_tool = ToolPlugin(
    name="search_motifs",
    category="taxonomy",
    description="Search compound motif labels (e.g. scaffold='Ar', substituent='Br' → 'Ar-Br'). Useful for taxonomy and classification.",
    prerequisites=[],
    fn=_search_motifs,
    llm_exposed=False,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

TAXONOMY_TOOLS = [
    search_reaction_types_tool,
    search_motifs_tool,
]
