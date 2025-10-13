"""
Utility functions for reagent taxonomy management.

Provides:
- CAS number validation and normalization
- Text tokenization and sanitization
- Identity resolution from CAS (PubChem, Cactus)
- Entry building and synonym deduplication
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set
from urllib.parse import quote

try:
    import requests
except Exception:  # pragma: no cover
    requests = None  # type: ignore[assignment]

DEFAULT_RESOLVER_TIMEOUT = 6.0

# Regex patterns for CAS numbers and tokens
TOKEN_PATTERN = re.compile(r"[a-z0-9]{3,}")
CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
CAS_INLINE_PATTERN = re.compile(r"\d{2,7}-\d{2}-\d")


def _sanitize_text(text: str) -> str:
    """Remove all non-alphanumeric characters and lowercase."""
    return re.sub(r"[^a-z0-9]+", "", (text or "").lower())


def _tokenize_text(text: str) -> Set[str]:
    """Extract tokens from text for matching."""
    tokens: Set[str] = set()
    if not text:
        return tokens
    lower = (text or "").lower()
    compact = re.sub(r"[^a-z0-9]+", "", lower)
    if len(compact) >= 3:
        tokens.add(compact)
    tokens.update(tok for tok in TOKEN_PATTERN.findall(lower) if len(tok) >= 3)
    return tokens


def tokenize_all(texts: Iterable[str]) -> Set[str]:
    """Extract all tokens from multiple texts."""
    tokens: Set[str] = set()
    for text in texts:
        tokens.update(_tokenize_text(text))
    return tokens


def _valid_cas_checksum(cas: str) -> bool:
    """Validate CAS number checksum."""
    digits = cas.replace("-", "")
    total = 0
    multiplier = 1
    for digit in reversed(digits[:-1]):
        total += int(digit) * multiplier
        multiplier += 1
    return total % 10 == int(digits[-1])


def normalize_cas(cas: str) -> str:
    """
    Normalize and validate CAS number.
    
    Args:
        cas: CAS number (with or without hyphens)
        
    Returns:
        Normalized CAS number (e.g., "14221-01-3")
        
    Raises:
        ValueError: If CAS is invalid
        
    Example:
        >>> normalize_cas("14221013")
        '14221-01-3'
        >>> normalize_cas("14221-01-3")
        '14221-01-3'
    """
    digits = re.sub(r"[^0-9]", "", (cas or ""))
    if len(digits) < 5:
        raise ValueError(f"CAS '{cas}' does not have enough digits")
    prefix = str(int(digits[:-3]))
    mid = digits[-3:-1]
    check = digits[-1]
    normalized = f"{prefix}-{mid}-{check}"
    if not CAS_PATTERN.fullmatch(normalized):
        raise ValueError(f"CAS '{cas}' is not in a valid format")
    if not _valid_cas_checksum(normalized):
        raise ValueError(f"CAS '{cas}' failed checksum validation")
    return normalized


def _http_get_json(session: Any, url: str, timeout: float) -> Optional[Dict[str, Any]]:
    """HTTP GET returning JSON."""
    try:
        response = session.get(url, timeout=timeout)
    except Exception:
        return None
    if response.status_code != 200:
        return None
    try:
        return response.json()
    except Exception:
        return None


def _http_get_text(session: Any, url: str, timeout: float) -> Optional[str]:
    """HTTP GET returning text."""
    try:
        response = session.get(url, timeout=timeout)
    except Exception:
        return None
    if response.status_code != 200:
        return None
    return response.text


def _normalized_cas_tokens(synonyms: Sequence[str]) -> Set[str]:
    """Extract normalized CAS numbers from synonyms."""
    tokens: Set[str] = set()
    for synonym in synonyms:
        if not synonym:
            continue
        inline_matches = CAS_INLINE_PATTERN.findall(synonym)
        if inline_matches:
            for match in inline_matches:
                try:
                    tokens.add(normalize_cas(match))
                except ValueError:
                    continue
            continue
        try:
            normalized = normalize_cas(synonym)
        except ValueError:
            continue
        tokens.add(normalized)
    return tokens


def _pubchem_cid_rn_map(session: Any, cids: Sequence[Any], timeout: float) -> Dict[Any, Set[str]]:
    """Get CAS registry numbers for PubChem CIDs."""
    if not cids:
        return {}
    cid_tokens: List[str] = []
    for cid in cids:
        try:
            cid_int = int(str(cid))
        except (TypeError, ValueError):
            continue
        if cid_int < 0:
            continue
        cid_tokens.append(str(cid_int))
    if not cid_tokens:
        return {}
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        + ",".join(cid_tokens)
        + "/xrefs/RN/JSON"
    )
    data = _http_get_json(session, url, timeout)
    mapping: Dict[Any, Set[str]] = {}
    info_list = data.get("InformationList", {}).get("Information", []) if data else []
    for block in info_list:
        cid_val = block.get("CID")
        rn_list = block.get("RN") or []
        normalized: Set[str] = set()
        for rn in rn_list:
            try:
                normalized.add(normalize_cas(rn))
            except ValueError:
                continue
        if normalized:
            mapping[cid_val] = normalized
    return mapping


def _resolve_via_pubchem(session: Any, cas: str, timeout: float) -> Optional[Dict[str, Any]]:
    """Resolve reagent identity from CAS using PubChem API."""
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{quote(cas)}"
    props_url = base + "/property/Title,IUPACName,IsomericSMILES,CanonicalSMILES/JSON"
    props = _http_get_json(session, props_url, timeout)
    entries = props.get("PropertyTable", {}).get("Properties", []) if props else []
    if not entries:
        return None

    syn_url = base + "/synonyms/JSON"
    syn_data = _http_get_json(session, syn_url, timeout)
    synonyms_by_cid: Dict[Any, List[str]] = {}
    if syn_data:
        info = syn_data.get("InformationList", {}).get("Information", [])
        for block in info:
            cid = block.get("CID")
            synonyms = block.get("Synonym", []) or []
            deduped_synonyms = dedupe_synonyms([syn for syn in synonyms if syn])
            if cid is not None:
                synonyms_by_cid[cid] = deduped_synonyms

    try:
        normalized_query = normalize_cas(cas)
    except ValueError:
        normalized_query = None

    selected_record = entries[0]
    selected_synonyms = synonyms_by_cid.get(selected_record.get("CID"), [])
    match_found = False
    if normalized_query:
        for record in entries:
            cid = record.get("CID")
            record_synonyms = synonyms_by_cid.get(cid, [])
            if not record_synonyms:
                continue
            if normalized_query in _normalized_cas_tokens(record_synonyms):
                selected_record = record
                selected_synonyms = record_synonyms
                match_found = True
                break
        if len(entries) > 1 and not match_found:
            rn_map = _pubchem_cid_rn_map(
                session, [record.get("CID") for record in entries], timeout
            )
            for record in entries:
                cid = record.get("CID")
                if cid is None:
                    continue
                rn_set = rn_map.get(cid, set())
                if normalized_query in rn_set:
                    selected_record = record
                    selected_synonyms = synonyms_by_cid.get(cid, [])
                    match_found = True
                    break
        if len(entries) > 1 and not match_found:
            return None

    smiles = selected_record.get("IsomericSMILES") or selected_record.get("CanonicalSMILES")
    primary_name = selected_record.get("Title") or selected_record.get("IUPACName")

    raw: List[str] = []
    if primary_name:
        raw.append(primary_name)
    raw.extend(selected_synonyms)
    deduped = dedupe_synonyms(raw)
    if not deduped and primary_name:
        deduped = [primary_name]
    if not deduped:
        return None
    if len(deduped) > 16:
        deduped = deduped[:16]

    name = primary_name or deduped[0]
    return {"name": name, "synonyms": deduped, "smiles": smiles}


def _resolve_via_cactus(session: Any, cas: str, timeout: float) -> Optional[Dict[str, Any]]:
    """Resolve reagent identity from CAS using Cactus NCI API."""
    base = f"https://cactus.nci.nih.gov/chemical/structure/{quote(cas)}"
    names_text = _http_get_text(session, base + "/names", timeout)
    synonyms = []
    if names_text:
        synonyms = [line.strip() for line in names_text.splitlines() if line.strip()]
    deduped = dedupe_synonyms(synonyms)
    if not deduped:
        return None
    smiles_text = _http_get_text(session, base + "/smiles", timeout)
    smiles = (smiles_text.strip() if smiles_text else None) or None
    if len(deduped) > 16:
        deduped = deduped[:16]
    return {"name": deduped[0], "synonyms": deduped, "smiles": smiles}


def resolve_identity_from_cas(
    cas: str,
    *,
    timeout: float = DEFAULT_RESOLVER_TIMEOUT,
    session: Optional[Any] = None,
) -> Optional[Dict[str, Any]]:
    """
    Resolve reagent name and SMILES from CAS number.
    
    Tries PubChem first, then Cactus NCI as fallback.
    
    Args:
        cas: CAS registry number
        timeout: HTTP request timeout in seconds
        session: Optional requests session (creates one if None)
        
    Returns:
        Dict with 'name', 'synonyms', 'smiles', 'source' or None
        
    Example:
        >>> info = resolve_identity_from_cas("14221-01-3")
        >>> print(info['name'])
        'Tetrakis(triphenylphosphine)palladium(0)'
    """
    if requests is None:
        return None
    if not cas:
        return None
    own_session = session is None
    sess = session or requests.Session()
    if own_session:
        sess.headers.setdefault("User-Agent", "ConditionAgent/TaxonomyGenerator")
    try:
        identity = _resolve_via_pubchem(sess, cas, timeout)
        if identity:
            identity["source"] = "pubchem"
            return identity
        identity = _resolve_via_cactus(sess, cas, timeout)
        if identity:
            identity["source"] = "cactus"
            return identity
        return None
    finally:
        if own_session:
            sess.close()


def build_entry(
    role: str,
    family: Dict[str, Any],
    cas: str,
    name: str,
    abbr: str,
    synonyms: List[str],
    smiles: Optional[str],
    numeric_features: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Build a taxonomy entry dictionary.
    
    Args:
        role: Reagent role (ligand, base, etc.)
        family: Family data dictionary
        cas: CAS registry number
        name: Primary name
        abbr: Abbreviation
        synonyms: List of synonyms
        smiles: Optional SMILES string
        numeric_features: Optional numeric features dict
        
    Returns:
        Complete entry dictionary with embedding text
    """
    from .taxonomy_store import build_embedding_text
    
    entry: Dict[str, Any] = {
        "name": name,
        "abbr": abbr,
        "cas": cas,
        "synonyms": synonyms,
    }
    if numeric_features:
        entry["numeric_features"] = numeric_features
    if smiles:
        entry["smiles"] = smiles
    alias_candidates = [syn for syn in synonyms if syn.lower() not in {name.lower(), abbr.lower()}]
    if alias_candidates:
        entry["aliases"] = alias_candidates
    entry["embedding_text"] = build_embedding_text(role, family, entry, synonyms)
    return entry


def dedupe_synonyms(values: Iterable[str]) -> List[str]:
    """
    Deduplicate synonyms (case-insensitive).
    
    Args:
        values: Iterable of synonym strings
        
    Returns:
        Deduplicated list preserving original case and order
        
    Example:
        >>> dedupe_synonyms(["PPh3", "pph3", "Triphenylphosphine"])
        ['PPh3', 'Triphenylphosphine']
    """
    result: List[str] = []
    seen: Set[str] = set()
    for value in values:
        if not value:
            continue
        clean = value.strip()
        if not clean:
            continue
        key = clean.lower()
        if key in seen:
            continue
        seen.add(key)
        result.append(clean)
    return result
