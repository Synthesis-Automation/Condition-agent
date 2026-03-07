"""
CAS number resolver — converts CAS Registry Numbers to compound identifiers.

Provides two complementary functions:

  resolve_cas(cas)   → full compound info (name, SMILES, formula, IUPAC, CID, synonyms)
  smiles_to_cas(smi) → look up the CAS number(s) for a known SMILES

Resolution cascade for resolve_cas:
  0. Local reagent registry (CSV) — instant, offline, no rate-limit
  1. PubChem   — CAS treated as synonym: CAS → CID → Title + IUPAC + SMILES + formula
  2. NCI CACTUS — CAS directly to SMILES + names

CAS number validation is performed before any API call using the standard
checksum algorithm (CASRN format XXXXXXX-YY-Z, check digit = sum(i*d) mod 10).

Usage:
    from chemtools.util.cas_resolver import resolve_cas, smiles_to_cas, is_cas

    # CAS → full info
    info = resolve_cas("50-78-2")
    # → {"name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    #    "iupac_name": "2-acetyloxybenzoic acid", "molecular_formula": "C9H8O4",
    #    "pubchem_cid": 2244, "synonyms": [...], "source": "pubchem"}

    # SMILES → CAS
    cas_list = smiles_to_cas("CC(=O)Oc1ccccc1C(=O)O")
    # → ["50-78-2"]

    # Validate format
    is_cas("50-78-2")   # True
    is_cas("50-78-3")   # False  (wrong check digit)
"""
from __future__ import annotations

import json
import random
import re
import time
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote_plus

# In-memory cache: {cas_normalized: result_dict}
_CAS_CACHE: Dict[str, Dict] = {}

_HTTP_TIMEOUT = 10
_RATE_DELAY   = 0.3

_CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")


# ---------------------------------------------------------------------------
# CAS validation
# ---------------------------------------------------------------------------

def is_cas(value: str) -> bool:
    """
    Return True if *value* is a syntactically valid CAS Registry Number.

    Validates both the format (XXXXXXX-YY-Z) and the checksum digit.
    Does NOT verify that the number exists in any database.
    """
    value = value.strip()
    if not _CAS_PATTERN.match(value):
        return False

    digits = value.replace("-", "")
    check_digit = int(digits[-1])
    body = digits[:-1]

    total = sum((i + 1) * int(d) for i, d in enumerate(reversed(body)))
    return total % 10 == check_digit


def normalize_cas(value: str) -> Optional[str]:
    """
    Normalize a CAS string to standard format (drop leading zeros, validate).

    Returns the normalized CAS string if valid, or None if invalid.
    """
    value = value.strip()
    # Handle common variants: no dashes, extra spaces
    digits_only = re.sub(r"[^0-9]", "", value)
    if len(digits_only) < 5:
        return None

    # Reconstruct CAS from pure digits: XXXXX-YY-Z
    reconstructed = f"{digits_only[:-3]}-{digits_only[-3:-1]}-{digits_only[-1]}"
    return reconstructed if is_cas(reconstructed) else None


# ---------------------------------------------------------------------------
# HTTP helper
# ---------------------------------------------------------------------------

def _get(url: str) -> Optional[str]:
    """HTTP GET with timeout. Returns response text or None."""
    try:
        import requests
        resp = requests.get(
            url,
            timeout=_HTTP_TIMEOUT,
            headers={"User-Agent": "ChemCoworker-CASResolver/1.0"},
        )
        if resp.status_code == 200:
            return resp.text
        if resp.status_code == 429:
            time.sleep(5)
            resp2 = requests.get(url, timeout=_HTTP_TIMEOUT,
                                 headers={"User-Agent": "ChemCoworker-CASResolver/1.0"})
            return resp2.text if resp2.status_code == 200 else None
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Helpers: SMILES canonicalization, validation
# ---------------------------------------------------------------------------

def _canonical(smiles: str) -> str:
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return smiles


def _validate_smiles(smiles: str) -> bool:
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return bool(smiles and smiles.strip())


# ---------------------------------------------------------------------------
# Backend 0: Local reagent registry (CSV, offline)
# ---------------------------------------------------------------------------

def _local_cas(cas: str) -> Optional[Dict]:
    """Search the local CSV reagent registry for a matching CAS number.

    Returns a result dict (same shape as the online backends) or None.
    """
    try:
        from ..reagent.lookup import load_reagent_database
        for reagent in load_reagent_database("*"):
            r_cas = normalize_cas((reagent.get("cas") or "").strip())
            if r_cas and r_cas == cas:
                smiles = (reagent.get("smiles") or "").strip()
                if smiles and _validate_smiles(smiles):
                    smiles = _canonical(smiles)
                abbrs: List[str] = [
                    a for a in (reagent.get("abbreviation") or []) if a
                ]
                return {
                    "name": (reagent.get("name") or "").strip(),
                    "iupac_name": "",
                    "molecular_formula": "",
                    "smiles": smiles,
                    "pubchem_cid": None,
                    "synonyms": abbrs,
                    "source": "local",
                }
    except Exception:
        pass
    return None


def _local_smiles_to_cas(smiles: str) -> List[str]:
    """Search the local registry for all CAS numbers matching a SMILES string."""
    try:
        from ..reagent.lookup import get_canonical_smiles_index
        reagent = get_canonical_smiles_index().get(smiles)
        if reagent:
            r_cas = (reagent.get("cas") or "").strip()
            if r_cas and is_cas(r_cas):
                return [r_cas]
    except Exception:
        pass
    return []


# ---------------------------------------------------------------------------
# Backend 1: PubChem (CAS as synonym → CID → properties)
# ---------------------------------------------------------------------------

def _pubchem_cas(cas: str) -> Optional[Dict]:
    """
    PubChem: treat CAS number as a compound synonym to get the CID,
    then fetch Title, IUPACName, MolecularFormula, and SMILES by CID.
    Also fetches synonyms (filtered to remove registry/catalog noise).
    """
    encoded = quote_plus(cas)

    # Step 1: CAS → CID
    cid_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/cids/JSON"
    )
    cid_text = _get(cid_url)
    if not cid_text:
        return None
    try:
        cid_list = json.loads(cid_text).get("IdentifierList", {}).get("CID", [])
        if not cid_list:
            return None
        cid = cid_list[0]
    except Exception:
        return None

    # Step 2: CID → Title + IUPACName + Formula + SMILES
    time.sleep(_RATE_DELAY)
    prop_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Title,IUPACName,MolecularFormula/JSON"
    )
    prop_text = _get(prop_url)
    if not prop_text:
        return None
    try:
        props = json.loads(prop_text).get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        p = props[0]
    except Exception:
        return None

    # Step 3: SMILES via separate endpoint (more reliable than property list)
    time.sleep(_RATE_DELAY)
    smiles_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/IsomericSMILES,CanonicalSMILES/JSON"
    )
    smiles = ""
    smiles_text = _get(smiles_url)
    if smiles_text:
        try:
            smiles_props = (
                json.loads(smiles_text)
                .get("PropertyTable", {})
                .get("Properties", [{}])
            )
            sp = smiles_props[0] if smiles_props else {}
            # PubChem may return the field as SMILES, IsomericSMILES, or CanonicalSMILES
            smiles = (
                sp.get("IsomericSMILES")
                or sp.get("CanonicalSMILES")
                or sp.get("SMILES")
                or sp.get("ConnectivitySMILES")
                or ""
            )
        except Exception:
            pass

    # Canonicalize with RDKit if we got a SMILES
    if smiles and _validate_smiles(smiles):
        smiles = _canonical(smiles)

    result: Dict = {
        "name": p.get("Title", ""),
        "iupac_name": p.get("IUPACName", ""),
        "molecular_formula": p.get("MolecularFormula", ""),
        "smiles": smiles,
        "pubchem_cid": cid,
        "synonyms": [],
        "source": "pubchem",
    }

    # Step 4: synonyms by CID (non-critical)
    try:
        time.sleep(_RATE_DELAY)
        syn_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/synonyms/JSON"
        )
        syn_text = _get(syn_url)
        if syn_text:
            info = (
                json.loads(syn_text)
                .get("InformationList", {})
                .get("Information", [{}])
            )
            all_syns = info[0].get("Synonym", []) if info else []
            result["synonyms"] = [
                s for s in all_syns
                if len(s) < 60 and not re.match(r"^[A-Z]{2,}\d", s)
            ][:10]
    except Exception:
        pass

    return result


# ---------------------------------------------------------------------------
# Backend 2: NCI CACTUS (CAS → SMILES + names directly)
# ---------------------------------------------------------------------------

def _cactus_cas(cas: str) -> Optional[Dict]:
    """CACTUS: resolve CAS directly to SMILES and names."""
    encoded = quote_plus(cas)

    # SMILES
    smiles_text = _get(
        f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/smiles"
    )
    smiles = ""
    if smiles_text and not smiles_text.lower().startswith("error"):
        raw = smiles_text.strip().split()[0]
        if _validate_smiles(raw):
            smiles = _canonical(raw)

    # IUPAC name
    iupac_text = _get(
        f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/iupac_name"
    )
    iupac = (
        iupac_text.strip()
        if iupac_text and not iupac_text.lower().startswith("error")
        else ""
    )

    # Names (synonyms)
    _catalog_re = re.compile(
        r"(_\d|^[A-Z]{2,}-\d|[A-Z]\d{4,}|NCGC|SDCC|AKOS|MFCD)"
    )
    synonyms: List[str] = []
    preferred_name = ""
    names_text = _get(
        f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/names"
    )
    if names_text and not names_text.lower().startswith("error"):
        all_names = [n.strip() for n in names_text.strip().splitlines() if n.strip()]
        clean = [n for n in all_names if not _catalog_re.search(n)]
        for n in clean[:20]:
            if len(n) <= 40 and not re.match(r"^\d", n):
                preferred_name = preferred_name or n
            synonyms.append(n)
        synonyms = synonyms[:8]

    if not smiles and not iupac and not preferred_name:
        return None

    return {
        "name": preferred_name,
        "iupac_name": iupac,
        "molecular_formula": "",  # CACTUS doesn't return formula cleanly
        "smiles": smiles,
        "pubchem_cid": None,
        "synonyms": synonyms,
        "source": "cactus",
    }


# ---------------------------------------------------------------------------
# Public API — CAS → compound info
# ---------------------------------------------------------------------------

def resolve_cas(cas: str) -> Dict:
    """
    Resolve a CAS Registry Number to full compound information.

    Validates the CAS checksum, then queries PubChem → CACTUS cascade.
    Returns name, SMILES, IUPAC name, molecular formula, PubChem CID,
    and a list of synonyms.

    Args:
        cas: CAS Registry Number string. Format: XXXXXXX-YY-Z
             Examples: "50-78-2" (aspirin), "58-08-2" (caffeine)
             Also accepts whitespace padding or variants without dashes.

    Returns:
        dict with keys: cas, name, iupac_name, smiles, molecular_formula,
        pubchem_cid, synonyms, source, valid_cas.
        source is "pubchem", "cactus", or "not_found".
    """
    cas = cas.strip()
    normalized = normalize_cas(cas) or cas  # keep original if normalization fails
    valid = is_cas(normalized)

    if normalized in _CAS_CACHE:
        return _CAS_CACHE[normalized]

    base = {
        "cas": normalized,
        "valid_cas": valid,
        "name": "",
        "iupac_name": "",
        "smiles": "",
        "molecular_formula": "",
        "pubchem_cid": None,
        "synonyms": [],
    }

    if not valid:
        result = {**base, "source": "not_found",
                  "error": f"'{cas}' is not a valid CAS number (format or checksum error)"}
        _CAS_CACHE[normalized] = result
        return result

    backends = [(_pubchem_cas, "pubchem"), (_cactus_cas, "cactus")]

    # Check local registry first — instant, no network, no rate-limit
    try:
        local_info = _local_cas(normalized)
        if local_info and (local_info.get("name") or local_info.get("smiles")):
            result = {**base, **local_info, "cas": normalized}
            _CAS_CACHE[normalized] = result
            return result
    except Exception:
        pass

    for backend_fn, _ in backends:
        try:
            time.sleep(_RATE_DELAY + random.uniform(0, 0.1))
            info = backend_fn(normalized)
            if info and (info.get("name") or info.get("smiles") or info.get("iupac_name")):
                result = {**base, **info, "cas": normalized}
                _CAS_CACHE[normalized] = result
                return result
        except Exception:
            continue

    result = {**base, "source": "not_found",
              "error": f"CAS '{normalized}' not found in PubChem or CACTUS"}
    _CAS_CACHE[normalized] = result
    return result


# ---------------------------------------------------------------------------
# Public API — SMILES → CAS number(s)
# ---------------------------------------------------------------------------

def smiles_to_cas(smiles: str) -> List[str]:
    """
    Look up CAS Registry Number(s) for a given SMILES via PubChem synonyms.

    Canonicalizes the SMILES, resolves to a PubChem CID, then fetches the
    raw synonym list and filters for strings matching the CAS number pattern.
    Returns all matching CAS numbers (usually 1, sometimes multiple for salts).

    Args:
        smiles: Molecule SMILES string.

    Returns:
        List of CAS number strings (may be empty if not found in PubChem).
    """
    canonical = _canonical(smiles)
    if not canonical:
        return []

    # Check local registry first (instant, no network)
    local_hits = _local_smiles_to_cas(canonical)
    if local_hits:
        return local_hits

    try:
        encoded = quote_plus(canonical)

        # Step 1: SMILES → CID
        cid_text = _get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            f"{encoded}/cids/JSON"
        )
        if not cid_text:
            # Retry with stereo-stripped SMILES
            stripped = re.sub(r"[/@\\]", "", canonical)
            cid_text = _get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
                f"{quote_plus(stripped)}/cids/JSON"
            )
        if not cid_text:
            return []

        cid_list = json.loads(cid_text).get("IdentifierList", {}).get("CID", [])
        if not cid_list:
            return []
        cid = cid_list[0]

        # Step 2: CID → raw synonyms (no filtering — we need the CAS entries)
        time.sleep(_RATE_DELAY)
        syn_text = _get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/synonyms/JSON"
        )
        if not syn_text:
            return []

        info = (
            json.loads(syn_text)
            .get("InformationList", {})
            .get("Information", [{}])
        )
        all_syns = info[0].get("Synonym", []) if info else []

        # Filter: valid CAS numbers only (format + checksum)
        return [s for s in all_syns if is_cas(s.strip())]

    except Exception:
        return []


def clear_cache() -> None:
    """Clear the in-memory CAS resolution cache."""
    _CAS_CACHE.clear()
