"""
Reagent database tools for ChemCoworker.

Wraps the chemtools reagent registry:
  - lookup_reagent        : find a reagent by name, abbreviation, or CAS
  - list_reagents_by_role : list all reagents for a given role
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from ._helpers import _error, _success, _to_jsonable
from ._base import ToolPlugin


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_reagent(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Return a compact, readable version of a reagent entry."""
    return {
        "name": entry.get("name", ""),
        "abbreviation": entry.get("abbreviation", []),
        "cas": entry.get("cas", ""),
        "smiles": entry.get("smiles", ""),
        "role": entry.get("role", ""),
        "family_id": entry.get("family_id", ""),
        "tag": entry.get("tag", ""),
    }


# ---------------------------------------------------------------------------
# Tool: lookup_reagent
# ---------------------------------------------------------------------------

def _lookup_reagent(query: str) -> Dict[str, Any]:
    """Search the reagent database by name, abbreviation, or CAS number.

    The reagent database contains catalysts, ligands, bases, solvents,
    additives, and other common synthetic chemistry reagents.

    Args:
        query: Reagent name, common abbreviation (e.g. 'dppf', 'DMAP', 'DMF'),
               or CAS number (e.g. '7440-05-3' for Pd).

    Returns:
        dict with matched reagents (name, abbreviation, CAS, SMILES, role, family).
    """
    try:
        from chemtools.reagent.lookup import _load_all_reagents

        all_reagents = _load_all_reagents()
        query_lower = query.lower().strip()

        matches: List[Dict[str, Any]] = []
        for entry in all_reagents:
            name = (entry.get("name") or "").lower()
            cas = (entry.get("cas") or "").lower()
            abbrevs = [a.lower() for a in (entry.get("abbreviation") or [])]

            if (
                query_lower in name
                or query_lower == cas
                or query_lower in abbrevs
                or any(query_lower in a for a in abbrevs)
            ):
                matches.append(_clean_reagent(entry))

        if not matches:
            return _success({
                "query": query,
                "matches": [],
                "note": f"No reagent found for '{query}'. Try a broader name or abbreviation.",
            })

        return _success({
            "query": query,
            "matches": matches[:20],
            "total_found": len(matches),
        })
    except Exception as exc:
        return _error(f"Reagent lookup failed: {exc}")


lookup_reagent_tool = ToolPlugin(
    name="lookup_reagent",
    category="reagent",
    description="Look up a reagent by name, abbreviation, or CAS number. Returns role, family, SMILES, and other details.",
    prerequisites=[],
    fn=_lookup_reagent,
)


# ---------------------------------------------------------------------------
# Tool: list_reagents_by_role
# ---------------------------------------------------------------------------

def _list_reagents_by_role(role: str, family: Optional[str] = None) -> Dict[str, Any]:
    """List available reagents for a given role, optionally filtered by family.

    Valid roles include: metal_catalyst, ligand, base, solvent, additive,
    oxidant, reductant, acid, fluoride_source, and more.

    Common role aliases are also accepted:
      'catalyst' → 'metal_catalyst', 'pd' → 'metal_catalyst', etc.

    Args:
        role: Reagent role (e.g. 'base', 'ligand', 'solvent', 'metal_catalyst').
        family: Optional family filter (e.g. 'phosphine', 'carbonate', 'amine').

    Returns:
        dict with reagent list for that role/family.
    """
    try:
        from chemtools.reagent.lookup import load_reagent_database

        entries = load_reagent_database(role)
        if not entries:
            return _success({
                "role": role,
                "family": family,
                "reagents": [],
                "note": f"No reagents found for role='{role}'. Check available roles.",
            })

        # Filter by family if requested
        if family:
            family_lower = family.lower()
            entries = [
                e for e in entries
                if family_lower in (e.get("family_id") or "").lower()
                or family_lower in (e.get("tag") or "").lower()
            ]

        cleaned = [_clean_reagent(e) for e in entries[:50]]
        return _success({
            "role": role,
            "family_filter": family or "",
            "reagents": cleaned,
            "total_found": len(entries),
        })
    except Exception as exc:
        return _error(f"Reagent listing failed: {exc}")


list_reagents_by_role_tool = ToolPlugin(
    name="list_reagents_by_role",
    category="reagent",
    description="List all reagents for a role (base, ligand, solvent, metal_catalyst, additive). Optional family filter.",
    prerequisites=[],
    fn=_list_reagents_by_role,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

REAGENT_TOOLS = [
    lookup_reagent_tool,
    list_reagents_by_role_tool,
]
