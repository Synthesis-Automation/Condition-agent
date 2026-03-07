"""
Chemical name → SMILES resolution utility.

Resolution chain (tried in priority order):
  0. Local reagent registry (CSV) — instant, offline, no rate-limit
  1. PubChem REST API   — most comprehensive; supports IUPAC, trivial, CAS, trade names
  2. OPSIN (Cambridge)  — excellent for systematic IUPAC names, free parser
  3. NCI CACTUS         — broad fallback for names not in PubChem/OPSIN

All returned SMILES are validated and canonicalized with RDKit before caching.
Results are held in a module-level in-memory cache for the session lifetime.

Usage:
    from chemtools.util.name_to_smiles import resolve_name, resolve_names

    smiles, source = resolve_name("4-methylbiphenyl")
    # → ("Cc1ccc(-c2ccccc2)cc1", "pubchem")

    results = resolve_names(["phenol", "toluene"])
    # → {"phenol": ("Oc1ccccc1", "pubchem"), "toluene": ("Cc1ccccc1", "pubchem")}
"""
from __future__ import annotations

import json
import random
import re
import time
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote_plus

# ---------------------------------------------------------------------------
# In-memory cache: {normalized_name_lower: (smiles_or_None, source_str)}
# ---------------------------------------------------------------------------
_NAME_CACHE: Dict[str, Tuple[Optional[str], str]] = {}

_HTTP_TIMEOUT = 10   # seconds per request
_RATE_DELAY   = 0.3  # minimum pause between API calls (polite crawling)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _clean_name(name: str) -> str:
    """Normalize a chemical name for API lookup."""
    name = re.sub(r"\s+", " ", str(name).strip())
    # Take only the first component when multiple are semicolon-separated
    name = re.split(r"\s*;\s*", name)[0].strip()
    # Strip metadata prefixes that would confuse APIs
    for prefix in ("drug:", "compound:", "chemical:", "cas:"):
        if name.lower().startswith(prefix):
            name = name[len(prefix):].strip()
    return name


def _get(url: str) -> Optional[str]:
    """HTTP GET with timeout. Returns response body text or None."""
    try:
        import requests
        resp = requests.get(
            url,
            timeout=_HTTP_TIMEOUT,
            headers={"User-Agent": "ChemCoworker-NameResolver/1.0"},
        )
        if resp.status_code == 200:
            return resp.text
        if resp.status_code == 404:
            return None
        # 429 rate-limit: brief pause then give up (single retry attempt)
        if resp.status_code == 429:
            time.sleep(5)
            resp2 = requests.get(url, timeout=_HTTP_TIMEOUT,
                                 headers={"User-Agent": "ChemCoworker-NameResolver/1.0"})
            return resp2.text if resp2.status_code == 200 else None
    except Exception:
        pass
    return None


def _validate_smiles(smiles: str) -> bool:
    """Return True if RDKit can parse the SMILES."""
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return bool(smiles and smiles.strip())


def _canonical(smiles: str) -> str:
    """Return RDKit canonical SMILES, or the original string if RDKit fails."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smiles


# ---------------------------------------------------------------------------
# Individual resolution backends
# ---------------------------------------------------------------------------

def _local(name: str) -> Optional[str]:
    """Search the local CSV reagent registry via a pre-built name index.

    Uses a plain case-insensitive comparison (NOT normalize_name) to avoid
    false positives from entries like '(Piperidinyl)aniline' whose
    normalize_name result collapses to 'aniline'.
    """
    try:
        from ..reagent.lookup import get_name_index
        reagent = get_name_index().get(name.strip().lower())
        if reagent:
            smiles = (reagent.get("smiles") or "").strip()
            return smiles if smiles and _validate_smiles(smiles) else None
    except Exception:
        pass
    return None


def _pubchem(name: str) -> Optional[str]:
    """PubChem compound name → IsomericSMILES."""
    encoded = quote_plus(name)
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/IsomericSMILES/JSON"
    )
    text = _get(url)
    if not text:
        return None
    try:
        data = json.loads(text)
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("IsomericSMILES")
    except Exception:
        pass
    return None


def _opsin(name: str) -> Optional[str]:
    """OPSIN IUPAC name parser (University of Cambridge)."""
    encoded = quote_plus(name)
    url = f"https://opsin.ch.cam.ac.uk/opsin/{encoded}.json"
    text = _get(url)
    if not text:
        return None
    try:
        data = json.loads(text)
        smiles = data.get("smiles", "")
        return smiles or None
    except Exception:
        pass
    return None


def _cactus(name: str) -> Optional[str]:
    """NCI CACTUS chemical structure resolver (broad fallback)."""
    encoded = quote_plus(name)
    url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/smiles"
    text = _get(url)
    if text and text.strip() and not text.strip().lower().startswith("error"):
        return text.strip().split()[0]  # CACTUS may return multiple; take first
    return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def resolve_name(name: str) -> Tuple[Optional[str], str]:
    """
    Resolve a chemical name to SMILES using a three-backend cascade.

    Tries PubChem → OPSIN → CACTUS, stopping at the first valid hit.
    Result is validated with RDKit and canonicalized before caching.

    Args:
        name: Chemical name (IUPAC, trivial, trade name, or CAS number).

    Returns:
        (smiles, source) tuple where:
          - smiles  : canonical SMILES string, or None if not found
          - source  : "pubchem", "opsin", "cactus", or "not_found"
    """
    if not name or not isinstance(name, str):
        return None, "not_found"

    cleaned = _clean_name(name)
    if not cleaned:
        return None, "not_found"

    cache_key = cleaned.lower()
    if cache_key in _NAME_CACHE:
        return _NAME_CACHE[cache_key]

    # Check local registry first — instant, no network, no rate-limit
    local_smiles = _local(cleaned)
    if local_smiles and _validate_smiles(local_smiles):
        result: Tuple[Optional[str], str] = (_canonical(local_smiles), "local")
        _NAME_CACHE[cache_key] = result
        return result

    backends = [
        (_pubchem, "pubchem"),
        (_opsin,   "opsin"),
        (_cactus,  "cactus"),
    ]

    for backend_fn, source_label in backends:
        try:
            time.sleep(_RATE_DELAY + random.uniform(0, 0.1))
            smiles = backend_fn(cleaned)
            if smiles and _validate_smiles(smiles):
                result: Tuple[Optional[str], str] = (_canonical(smiles), source_label)
                _NAME_CACHE[cache_key] = result
                return result
        except Exception:
            continue

    _NAME_CACHE[cache_key] = (None, "not_found")
    return None, "not_found"


def resolve_names(names: List[str]) -> Dict[str, Tuple[Optional[str], str]]:
    """
    Resolve multiple chemical names to SMILES (sequential, rate-limited).

    Args:
        names: List of chemical names to resolve.

    Returns:
        {name: (smiles, source)} mapping for each input name.
    """
    return {name: resolve_name(name) for name in names}


def clear_cache() -> None:
    """Clear the in-memory name resolution cache."""
    _NAME_CACHE.clear()
