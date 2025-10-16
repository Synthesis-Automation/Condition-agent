"""Catalyst class detection and matching logic.

Provides functionality to classify reactions by catalyst type (metal, enzyme, etc.)
and match precedents based on catalyst class.
"""
from typing import Dict, Any, Optional


_METAL_NAME_TO_SYMBOL = {
    # common transition metals and aliases
    "palladium": "Pd", "pd": "Pd",
    "nickel": "Ni", "ni": "Ni",
    "cobalt": "Co", "co": "Co",
    "copper": "Cu", "cu": "Cu",
    "iron": "Fe", "fe": "Fe",
    "ruthenium": "Ru", "ru": "Ru",
    "rhodium": "Rh", "rh": "Rh",
    "iridium": "Ir", "ir": "Ir",
    "gold": "Au", "au": "Au",
    "silver": "Ag", "ag": "Ag",
    "zinc": "Zn", "zn": "Zn",
    "magnesium": "Mg", "mg": "Mg",
    "manganese": "Mn", "mn": "Mn",
    "chromium": "Cr", "cr": "Cr",
    "vanadium": "V", "v": "V",
    "titanium": "Ti", "ti": "Ti",
    "zirconium": "Zr", "zr": "Zr",
    "molybdenum": "Mo", "mo": "Mo",
    "tungsten": "W", "w": "W",
    "rhenium": "Re", "re": "Re",
    "osmium": "Os", "os": "Os",
    "platinum": "Pt", "pt": "Pt",
}

_METAL_SYMBOLS = set(_METAL_NAME_TO_SYMBOL.values())
_ENZYME_KEYWORDS = {"enzyme", "protein", "lipase", "oxidase", "dehydrogenase", "transferase"}


def _normalize_symbol(token: str) -> Optional[str]:
    """Normalize a metal name/symbol to standard symbol (e.g., 'palladium' -> 'Pd')."""
    t = (token or "").strip()
    if not t:
        return None
    lo = t.lower()
    if lo in _METAL_NAME_TO_SYMBOL:
        return _METAL_NAME_TO_SYMBOL[lo]
    # Try exact symbol match case-insensitively
    up = t[0].upper() + t[1:].lower()
    if up in _METAL_SYMBOLS:
        return up
    return None


def _row_catalyst_class(row: Dict[str, Any]) -> str:
    """Heuristically classify a precedent row into a catalyst class.

    Returns one of metal symbols (e.g., 'Pd', 'Ni', ...), 'enzyme', 'organo_catalyst', or 'other'.
    
    Note: Bases and acids (e.g., "Base: K2CO3", "Acid: HCl") are classified as 'other'
    because they are stoichiometric reagents, not catalysts.
    """
    # 1) Use condition_core like 'Pd/XPhos'
    core = (row.get("condition_core") or "").strip()
    if core:
        # Check if it's a base/acid reagent (not a true catalyst)
        if core.startswith("Base:") or core.startswith("Acid:"):
            return "other"  # Bases/acids are reagents, not catalysts
        
        head = core.split("/", 1)[0].strip()
        sym = _normalize_symbol(head)
        if sym:
            return sym
    
    # 2) Scan full_system names for explicit metal names
    fs = row.get("full_system")
    if isinstance(fs, list):
        names = [str((it or {}).get("name") or "") for it in fs]
        text = " ".join(names).lower()
        # enzyme detection first
        if any(k in text for k in _ENZYME_KEYWORDS):
            return "enzyme"
        # metal names
        for key, sym in _METAL_NAME_TO_SYMBOL.items():
            if len(key) <= 2:  # skip symbol aliases here; handled below
                continue
            if key in text:
                return sym
        # symbol occurrence as word-ish
        for sym in _METAL_SYMBOLS:
            if f" {sym.lower()}" in text or f"({sym.lower()}" in text or f"[{sym.lower()}" in text:
                return sym
    
    # 3) Fallback enzyme detection in catalyst dict
    cat = row.get("catalyst") or {}
    if isinstance(cat, dict):
        nm = str(cat.get("name") or "").lower()
        if any(k in nm for k in _ENZYME_KEYWORDS):
            return "enzyme"
    
    # 4) If no metal detected, check for true organocatalysts
    # True organocatalysts: proline, DMAP, DBU, thioureas, squaramides, etc.
    # But NOT simple bases/acids which are used stoichiometrically
    if core and not core.startswith("Base:") and not core.startswith("Acid:"):
        # Only classify as organo_catalyst if it's not a common base/acid
        if (isinstance(fs, list) and fs):
            return "organo_catalyst"
    
    return "other"


def _match_catalyst_class(selected: str, row_cls: str) -> bool:
    """Check if a row's catalyst class matches the selected filter.
    
    Args:
        selected: Desired catalyst class (metal symbol, 'enzyme', 'organo_catalyst', etc.)
        row_cls: Catalyst class of the precedent row
        
    Returns:
        True if the row matches the filter, False otherwise
    """
    sel = (selected or "").strip().lower()
    if not sel:
        return True
    if sel in {"organo_catalyst", "enzyme", "other"}:
        return row_cls == sel
    sym = _normalize_symbol(sel)
    if sym:
        return row_cls == sym
    return True  # unknown filter -> do not exclude
