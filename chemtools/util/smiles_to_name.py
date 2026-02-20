"""
SMILES → chemical name resolution utility.

Resolution chain (tried in priority order):
  1. PubChem REST API  — returns IUPAC name, preferred name, and synonyms
  2. NCI CACTUS        — fallback IUPAC name and common name list

The input SMILES is canonicalized with RDKit before lookup to maximize
API match rates. Results are cached in-memory for the session.

Usage:
    from chemtools.util.smiles_to_name import resolve_smiles

    result = resolve_smiles("CC(=O)Oc1ccccc1C(=O)O")
    # → {
    #     "iupac_name": "2-acetyloxybenzoic acid",
    #     "preferred_name": "Aspirin",
    #     "synonyms": ["aspirin", "acetylsalicylic acid", ...],
    #     "molecular_formula": "C9H8O4",
    #     "pubchem_cid": 2244,
    #     "source": "pubchem",
    #   }
"""
from __future__ import annotations

import json
import random
import re
import time
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote_plus

# In-memory cache: {canonical_smiles: result_dict}
_SMILES_CACHE: Dict[str, Dict] = {}

_HTTP_TIMEOUT = 10
_RATE_DELAY   = 0.3


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _canonical(smiles: str) -> str:
    """Return RDKit canonical SMILES, or original if RDKit fails."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smiles


def _get(url: str) -> Optional[str]:
    """HTTP GET with timeout. Returns response text or None."""
    try:
        import requests
        resp = requests.get(
            url,
            timeout=_HTTP_TIMEOUT,
            headers={"User-Agent": "ChemCoworker-SMILESResolver/1.0"},
        )
        if resp.status_code == 200:
            return resp.text
        if resp.status_code == 429:
            time.sleep(5)
            resp2 = requests.get(url, timeout=_HTTP_TIMEOUT,
                                 headers={"User-Agent": "ChemCoworker-SMILESResolver/1.0"})
            return resp2.text if resp2.status_code == 200 else None
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Individual backends
# ---------------------------------------------------------------------------

def _pubchem(smiles: str) -> Optional[Dict]:
    """PubChem SMILES → CID → Title (common name) + IUPACName + synonyms.

    Two-step: resolve SMILES to CID first, then fetch properties by CID.
    This avoids PubChem's SMILES-as-property-input limitations.
    """
    encoded = quote_plus(smiles)

    # Step 1: resolve SMILES → CID
    cid_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{encoded}/cids/JSON"
    )
    cid_text = _get(cid_url)
    if not cid_text:
        return None
    try:
        cid_data = json.loads(cid_text)
        cid_list = cid_data.get("IdentifierList", {}).get("CID", [])
        if not cid_list:
            return None
        cid = cid_list[0]
    except Exception:
        return None

    # Step 2: fetch Title + IUPACName + MolecularFormula by CID
    time.sleep(_RATE_DELAY)
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Title,IUPACName,MolecularFormula/JSON"
    )
    prop_text = _get(prop_url)
    if not prop_text:
        return None
    try:
        data = json.loads(prop_text)
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        p = props[0]
    except Exception:
        return None

    result = {
        "iupac_name": p.get("IUPACName", ""),
        "preferred_name": p.get("Title", ""),        # "Aspirin", "Caffeine", etc.
        "molecular_formula": p.get("MolecularFormula", ""),
        "pubchem_cid": cid,
        "synonyms": [],
        "source": "pubchem",
    }

    # Step 3: fetch top synonyms by CID (non-critical, best effort)
    try:
        time.sleep(_RATE_DELAY)
        syn_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/synonyms/JSON"
        )
        syn_text = _get(syn_url)
        if syn_text:
            syn_data = json.loads(syn_text)
            info = syn_data.get("InformationList", {}).get("Information", [])
            if info:
                all_syns = info[0].get("Synonym", [])
                # Keep first 8, skip long registry/catalog entries
                result["synonyms"] = [
                    s for s in all_syns
                    if len(s) < 60 and not re.match(r"^\d{2,7}-\d{2}-\d$", s)
                ][:8]
    except Exception:
        pass

    return result


def _cactus(smiles: str) -> Optional[Dict]:
    """NCI CACTUS SMILES → IUPAC name + name list."""
    encoded = quote_plus(smiles)

    iupac_url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/iupac_name"
    iupac_text = _get(iupac_url)
    iupac = iupac_text.strip() if iupac_text and not iupac_text.lower().startswith("error") else ""

    # Also try common names
    names_url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/names"
    names_text = _get(names_url)
    synonyms: List[str] = []
    preferred = ""
    if names_text and not names_text.lower().startswith("error"):
        all_names = [n.strip() for n in names_text.strip().splitlines() if n.strip()]
        # Skip catalog/registry identifiers (contain _, all-caps codes, CAS-like patterns)
        _catalog_re = re.compile(
            r"(_\d|^[A-Z]{2,}-\d|^\d{2,7}-\d{2}-\d$|[A-Z]\d{4,}|NCGC|SDCC|AKOS|MFCD)"
        )
        clean_names = [n for n in all_names if not _catalog_re.search(n)]
        for n in clean_names[:20]:
            # Prefer short names without special chars as the common name
            if len(n) <= 40 and not re.match(r"^\d", n) and n == n.rstrip(")"):
                preferred = preferred or n
            synonyms.append(n)
        synonyms = synonyms[:8]

    if not iupac and not preferred:
        return None

    return {
        "iupac_name": iupac,
        "preferred_name": preferred,
        "molecular_formula": "",
        "pubchem_cid": None,
        "synonyms": synonyms,
        "source": "cactus",
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def resolve_smiles(smiles: str) -> Dict:
    """
    Resolve a SMILES string to chemical names using PubChem → CACTUS cascade.

    The SMILES is canonicalized with RDKit before lookup. Results include
    the IUPAC name, preferred (common) name, synonyms, molecular formula,
    and PubChem CID (when available).

    Args:
        smiles: SMILES string (molecule or reaction component).
                Reaction SMILES (containing >>) are not supported.

    Returns:
        dict with keys:
          - iupac_name       : IUPAC systematic name (may be empty)
          - preferred_name   : Common/trade name (may be empty)
          - synonyms         : List of up to 8 alternative names
          - molecular_formula: e.g. "C9H8O4"
          - pubchem_cid      : int or None
          - source           : "pubchem", "cactus", or "not_found"
          - canonical_smiles : canonicalized input SMILES
    """
    if not smiles or not isinstance(smiles, str):
        return {"source": "not_found", "error": "Empty SMILES"}

    # Strip reaction arrow components if accidentally passed
    mol_smiles = smiles.split(">>")[0].split(">")[0].strip()
    canonical = _canonical(mol_smiles)

    if canonical in _SMILES_CACHE:
        return _SMILES_CACHE[canonical]

    backends = [
        (_pubchem, "pubchem"),
        (_cactus,  "cactus"),
    ]

    for backend_fn, _ in backends:
        try:
            time.sleep(_RATE_DELAY + random.uniform(0, 0.1))
            result = backend_fn(canonical)
            if result and (result.get("iupac_name") or result.get("preferred_name")):
                result["canonical_smiles"] = canonical
                _SMILES_CACHE[canonical] = result
                return result
        except Exception:
            continue

    fallback = {
        "iupac_name": "",
        "preferred_name": "",
        "synonyms": [],
        "molecular_formula": "",
        "pubchem_cid": None,
        "source": "not_found",
        "canonical_smiles": canonical,
    }
    _SMILES_CACHE[canonical] = fallback
    return fallback


def resolve_smiles_batch(smiles_list: List[str]) -> Dict[str, Dict]:
    """
    Resolve multiple SMILES strings to names (sequential, rate-limited).

    Returns:
        {smiles: result_dict} mapping.
    """
    return {s: resolve_smiles(s) for s in smiles_list}


def clear_cache() -> None:
    """Clear the in-memory SMILES→name cache."""
    _SMILES_CACHE.clear()
