"""Core utility functions for precedent search.

Includes family name normalization, bin parsing, and other helper utilities.
"""
from typing import Dict, Tuple, Optional


def _family_text(family: str) -> str:
    """Map API family tokens to canonical family names (new naming convention)."""
    f = (family or "").strip()
    fl = f.lower()

    # Taxonomy v2 identifiers
    if fl == "c_n_cross_coupling":
        return "C_N_Coupling"
    if fl == "suzuki_miyaura":
        return "Suzuki"
    if fl == "sonogashira":
        return "Sonogashira_coupling"
    if fl == "heck":
        return "HeckMizoroki_coupling"
    if fl == "snar_cn":
        return "SNAr-CN"
    if fl == "amide_coupling":
        return "Amide_formation"
    
    # New systematic naming
    if fl in {"c_n_coupling_cu", "c_n_coupling_cu_ullmann"}:
        return "C_N_Coupling_Cu"
    if fl in {"c_n_coupling_pd", "c_n_coupling_pd_buchwald"}:
        return "C_N_Coupling_Pd"
    if fl == "c_n_coupling_ni":
        return "C_N_Coupling_Ni"
    
    # Legacy naming → new naming
    if fl in {"ullmann_cn", "ullmann c鈥搉", "ullmann c-n", "ullmann"}:
        return "C_N_Coupling_Cu"
    if fl in {"buchwald_cn", "buchwald c鈥搉", "buchwald c-n", "buchwald"}:
        return "C_N_Coupling_Pd"
    
    # Amide formation aliases
    if fl in {"amide_coupling", "amidation", "amide"}:
        return "Amide_formation"
    
    return f


def _proto_family_id(family_txt: str) -> str:
    """Normalize family text for use in prototype_id."""
    txt = str(family_txt).replace(' ', '_')
    for old in ('–', '—', '−'):
        txt = txt.replace(old, '-')
    return txt.replace('/', '_')


def _parse_bin(bin_str: str) -> Dict[str, str]:
    """Parse a bin string like 'LG:Br|NUC:aniline' into a dict."""
    out: Dict[str, str] = {}
    if not bin_str:
        return out
    for part in str(bin_str).split('|'):
        if ':' in part:
            k, v = part.split(':', 1)
            out[k.strip()] = v.strip()
    return out


def _parse_core_tokens(core_text: str) -> Tuple[str, str]:
    """Parse a condition_core string like 'Pd/XPhos' into (metal, ligand) tokens (lowercased).

    Returns (metal, ligand) where any missing part is an empty string.
    """
    t = (core_text or "").strip()
    if not t:
        return "", ""
    if "/" in t:
        a, b = t.split("/", 1)
        return (a or "").strip().lower(), (b or "").strip().lower()
    return (t.strip().lower(), "")


def _norm_family(fam: Optional[str]) -> Optional[str]:
    """Normalize family string using _family_text, or return None."""
    if fam is None:
        return None
    f = (fam or "").strip()
    if not f:
        return None
    return _family_text(f)
