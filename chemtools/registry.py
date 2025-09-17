from __future__ import annotations
"""
chemtools.registry

Lightweight in-memory registry built from the compound taxonomy (data/compound_taxonomy).

- Default source: taxonomy JSON files (ligand, base, solvent, coupling_reagent, catalysts_precursor)
- Override path: set CHEMTOOLS_TAXONOMY_DIR to an alternate taxonomy directory
- Optional custom registry: set CHEMTOOLS_REGISTRY_PATH to a JSON/JSONL to load instead

No dependency on a merged CAS registry JSONL by default.
"""

import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple
import json as _json

# Local optional properties enrichment
try:
    from . import properties as _properties
except Exception:  # pragma: no cover
    _properties = None  # type: ignore


# --- Roles mapping ---
ROLE_MAP = {
    "ligand": "LIGAND",
    "base": "BASE",
    "solvent": "SOLVENT",
    "catalyst_core": "CATALYST",
    "preformed metal–ligand catalyst": "CATALYST",
    "preformed metal-ligand catalyst": "CATALYST",
    "metal": "CATALYST",
    "catalyst": "CATALYST",
    # Map various additives to ADDITIVE
    "oxidant": "ADDITIVE",
    "activator": "ADDITIVE",
    "acid": "ADDITIVE",
    "salt": "ADDITIVE",
    "additive": "ADDITIVE",
    "reductant": "ADDITIVE",
}


def _normalize_alias(s: str) -> str:
    """Loose canonical token for alias matching.

    - Lowercase
    - Normalize dashes and strip punctuation/whitespace
    - Collapse to alphanumerics plus limited symbols (letters+digits only)
    """
    s0 = (s or "").strip().lower()
    if not s0:
        return ""
    s0 = (
        s0.replace("–", "-")
        .replace("—", "-")
        .replace("−", "-")
    )
    # Remove all non-alphanumeric characters
    return re.sub(r"[^a-z0-9]+", "", s0)


def _looks_like_cas(s: str) -> bool:
    s = (s or "").strip()
    return bool(re.fullmatch(r"\d{2,7}-\d{2}-\d", s))


def _maybe_canonicalize_cas(s: str) -> Optional[str]:
    """Attempt to canonicalize CAS numbers.

    Accepts inputs like '108-88-3' or raw digits '108883' (best-effort).
    Returns a hyphenated CAS if pattern is plausible; otherwise None.
    """
    s = (s or "").strip()
    if _looks_like_cas(s):
        return s
    digits = re.sub(r"[^0-9]", "", s)
    if len(digits) < 5:
        return None
    # Heuristic split: last is check digit, prior two are middle, rest prefix
    prefix = digits[:-3]
    mid = digits[-3:-1]
    check = digits[-1]
    if not prefix or len(mid) != 2:
        return None
    return f"{int(prefix)}-{mid}-{check}"


class _RegistryIndex:
    def __init__(self) -> None:
        self.by_uid: Dict[str, Dict[str, Any]] = {}
        self.alias_to_uid: Dict[str, str] = {}
        self.uid_to_aliases: Dict[str, List[str]] = {}

    def add_record(self, rec: Dict[str, Any]) -> None:
        uid = rec.get("uid")
        if not uid:
            return
        self.by_uid.setdefault(uid, rec)
        aliases = rec.get("aliases", []) or []
        # Maintain return-time alias list (dedup, preserve order best-effort)
        seen: set[str] = set()
        out_aliases: List[str] = []
        for a in [rec.get("name"), rec.get("abbreviation"), rec.get("token"), rec.get("generic_core"), uid, *aliases]:
            if not a:
                continue
            aa = str(a)
            if aa in seen:
                continue
            seen.add(aa)
            out_aliases.append(aa)
            na = _normalize_alias(aa)
            if na:
                # First writer wins to avoid noisy overwrites
                self.alias_to_uid.setdefault(na, uid)
        self.uid_to_aliases[uid] = out_aliases


_INDEX: Optional[_RegistryIndex] = None


def _role_from_record(r: Dict[str, Any]) -> str:
    # Prefer explicit hint if maps cleanly
    hint = (r.get("category_hint") or "").strip().upper()
    if hint in {"CATALYST", "LIGAND", "BASE", "SOLVENT", "ADDITIVE"}:
        return hint
    # Fall back to compound_type mapping
    ctype = (r.get("compound_type") or "").strip().lower()
    return ROLE_MAP.get(ctype, "ADDITIVE" if ctype else "ADDITIVE")


def _taxonomy_dir() -> Optional[str]:
    env = os.environ.get("CHEMTOOLS_TAXONOMY_DIR")
    if env and os.path.isdir(env):
        return env
    try:
        base = os.path.dirname(os.path.dirname(__file__))
        d = os.path.join(base, "data", "compound_taxonomy")
        return d if os.path.isdir(d) else None
    except Exception:
        return None


def _add_alias(idx: _RegistryIndex, uid: str, text: Optional[str]) -> None:
    if not text:
        return
    # Route via add_record’s alias path by appending to record later; for indexing, update alias_to_uid directly here
    na = _normalize_alias(str(text))
    if na:
        idx.alias_to_uid.setdefault(na, uid)


def _load_registry_from_taxonomy() -> _RegistryIndex:
    idx = _RegistryIndex()

    def load(rel: str) -> Dict[str, Any] | None:
        d = _taxonomy_dir()
        if not d:
            return None
        p = os.path.join(d, rel)
        if not os.path.exists(p):
            return None
        try:
            with open(p, "r", encoding="utf-8") as f:
                return _json.load(f)
        except Exception:
            return None

    # Ligands
    lig = load("taxonomy_ligand.json")
    if lig:
        for fam in lig.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "LIGAND",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": None,
                    "compound_type": "ligand",
                    "aliases": (em.get("synonyms") or []),
                }
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)

    # Bases
    bas = load("taxonomy_base.json")
    if bas:
        for fam in bas.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "BASE",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": None,
                    "compound_type": "base",
                    "aliases": (em.get("synonyms") or []),
                }
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)

    # Solvents
    solv = load("taxonomy_solvent.json")
    if solv:
        for fam in solv.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "SOLVENT",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": None,
                    "compound_type": "solvent",
                    "aliases": (em.get("synonyms") or []),
                }
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)

    # Coupling reagents -> ADDITIVE role, explicit type
    cou = load("taxonomy_coupling_reagent.json")
    if cou:
        for fam in cou.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "ADDITIVE",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": None,
                    "compound_type": "coupling_reagent",
                    "aliases": (em.get("synonyms") or []),
                }
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)

    # Reductants -> ADDITIVE role with reductant compound_type
    red = load("taxonomy_reductant.json")
    if red:
        for fam in red.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            fam_class = fam.get("class")
            mechanism = fam.get("mechanism")
            strength_index = fam.get("strength_index")
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "ADDITIVE",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": None,
                    "compound_type": "reductant",
                    "aliases": (em.get("synonyms") or []),
                }
                if fam_class:
                    rec["reductant_class"] = fam_class
                if mechanism:
                    rec["mechanism"] = mechanism
                if strength_index is not None:
                    rec["strength_index"] = strength_index
                if em.get("recommended_pairs"):
                    rec["recommended_pairs"] = em["recommended_pairs"]
                if em.get("properties"):
                    rec["properties"] = em["properties"]
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)
                _add_alias(idx, uid, fam_class)

    # Catalysts precursors -> CATALYST, generic_core from family metal
    cat = load("taxonomy_catalysts_precursor.json")
    if cat:
        for fam in cat.get("families", []) or []:
            fam_label = fam.get("label")
            fam_id = fam.get("family_id")
            metal = (fam.get("metal") or "").strip()
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                if not uid:
                    continue
                name = (em.get("name") or uid).strip()
                abbr = (em.get("abbr") or None)
                token = (abbr or name or uid)
                rec = {
                    "uid": uid,
                    "cas": uid,
                    "role": "CATALYST",
                    "name": name,
                    "abbreviation": abbr,
                    "token": token,
                    "generic_core": metal or None,
                    "compound_type": "catalyst_core",
                    "aliases": (em.get("synonyms") or []),
                }
                idx.add_record(rec)
                _add_alias(idx, uid, fam_label)
                _add_alias(idx, uid, fam_id)

    return idx


def _load_registry() -> _RegistryIndex:
    global _INDEX
    if _INDEX is not None:
        return _INDEX
    idx = _RegistryIndex()
    # Allow override via environment variable for testing/alternate datasets
    env_path = os.environ.get("CHEMTOOLS_REGISTRY_PATH")
    if env_path and os.path.exists(env_path):
        # JSON/JSONL fallback path if explicitly provided
        try:
            if env_path.lower().endswith(".json"):
                with open(env_path, "r", encoding="utf-8") as f:
                    data = _json.load(f)
                if isinstance(data, list):
                    items = data
                elif isinstance(data, dict):
                    # assume {uid: rec}
                    items = [dict(v, uid=k) for k, v in data.items()]
                else:
                    items = []
            else:
                items = []
                with open(env_path, "r", encoding="utf-8") as f:
                    for line in f:
                        s = line.strip()
                        if not s:
                            continue
                        try:
                            items.append(_json.loads(s))
                        except Exception:
                            continue
            for raw in items:
                cas = (raw.get("cas") or "").strip()
                uid = cas or (raw.get("uid") or "").strip()
                if not uid:
                    continue
                role_raw = (raw.get("role") or "").strip().upper()
                role = role_raw if role_raw in {"CATALYST", "LIGAND", "BASE", "SOLVENT", "ADDITIVE"} else _role_from_record(raw)
                name = raw.get("name") or uid
                rec: Dict[str, Any] = {
                    "uid": uid,
                    "cas": cas or None,
                    "role": role,
                    "name": name,
                    "abbreviation": raw.get("abbreviation") or None,
                    "token": raw.get("token") or None,
                    "generic_core": raw.get("generic_core") or None,
                    "compound_type": raw.get("compound_type") or None,
                    "aliases": [],
                }
                aliases: List[str] = []
                for k in ("name", "abbreviation", "token", "generic_core", "uid", "cas"):
                    v = raw.get(k)
                    if v:
                        aliases.append(str(v))
                rec["aliases"] = aliases
                idx.add_record(rec)
            _INDEX = idx
            return idx
        except Exception:
            # Fall through to taxonomy
            pass
    # Default: taxonomy-based registry
    _INDEX = _load_registry_from_taxonomy()
    return _INDEX


def categories() -> Dict[str, List[str]]:
    """List distinct roles and compound types present in the registry.

    Returns {roles: [...], compound_types: [...]} with sorted unique values.
    """
    idx = _load_registry()
    roles: set[str] = set()
    ctypes: set[str] = set()
    for _, rec in idx.by_uid.items():
        r = (rec.get("role") or "").strip().upper()
        if r:
            roles.add(r)
        ct = (rec.get("compound_type") or "").strip()
        if ct:
            ctypes.add(ct)
    return {
        "roles": sorted(roles),
        "compound_types": sorted(ctypes),
    }


def _enrich_props(uid: str, name: Optional[str] = None) -> Dict[str, Any]:
    props: Dict[str, Any] = {}
    if _properties is None:
        return props
    # Try lookup by uid (CAS), then by provided name token
    for q in [uid, name or ""]:
        q = (q or "").strip()
        if not q:
            continue
        try:
            # Avoid recursive calls back into registry.resolve while enriching
            res = _properties.lookup(q, allow_registry=False)
        except Exception:
            continue
        if res and res.get("found") and res.get("record"):
            # Copy all fields except uid and role (role is handled from registry)
            rec = dict(res["record"])  # type: ignore
            rec.pop("uid", None)
            rec.pop("role", None)
            props.update(rec)
            break
    return props


def resolve(query: str) -> Dict[str, Any]:
    """Resolve a registry item by CAS, name, alias, or token.

    Returns a record with fields: {uid, role, name, aliases, props}
    If not found, returns {error: 'NOT_FOUND'}
    """
    q_raw = (query or "").strip()
    if not q_raw:
        return {"error": "NOT_FOUND"}

    idx = _load_registry()

    # 1) Try CAS direct or canonicalized
    cas = q_raw
    if _looks_like_cas(cas):
        rec = idx.by_uid.get(cas)
        if rec:
            uid = rec["uid"]
            return {
                "uid": uid,
                "role": rec["role"],
                "name": rec["name"],
                "aliases": idx.uid_to_aliases.get(uid, []),
                "props": _enrich_props(uid, rec.get("name")),
            }
    can = _maybe_canonicalize_cas(q_raw)
    if can and can in idx.by_uid:
        rec = idx.by_uid[can]
        uid = rec["uid"]
        return {
            "uid": uid,
            "role": rec["role"],
            "name": rec["name"],
            "aliases": idx.uid_to_aliases.get(uid, []),
            "props": _enrich_props(uid, rec.get("name")),
        }

    # 2) Alias-based lookup (normalized)
    q_norm = _normalize_alias(q_raw)
    if q_norm:
        uid = idx.alias_to_uid.get(q_norm)
        if uid and uid in idx.by_uid:
            rec = idx.by_uid[uid]
            return {
                "uid": uid,
                "role": rec["role"],
                "name": rec["name"],
                "aliases": idx.uid_to_aliases.get(uid, []),
                "props": _enrich_props(uid, rec.get("name")),
            }

    return {"error": "NOT_FOUND"}


# Convenience: simple multi-lookup
def resolve_all(queries: List[str]) -> List[Dict[str, Any]]:
    return [resolve(q) for q in queries]


# Simple in-memory search over the loaded registry. Useful for filtering
# by role or by the source "compound_type" field from the JSONL.
def search(
    q: Optional[str] = None,
    *,
    role: Optional[str] = None,
    compound_type: Optional[str] = None,
    limit: int = 50,
) -> List[Dict[str, Any]]:
    """Search registry records.

    - q: optional substring (case-insensitive) matched against name/token/abbreviation/generic_core/uid.
    - role: one of CATALYST|LIGAND|BASE|SOLVENT|ADDITIVE (case-insensitive).
    - compound_type: substring match against the JSONL's compound_type field (case-insensitive).
    - limit: cap number of returned items.

    Returns a list of minimal records: {uid, role, name, compound_type?}
    """
    idx = _load_registry()
    ql = (q or "").strip().lower()
    role_u = (role or "").strip().upper()
    ct_sub = (compound_type or "").strip().lower()

    out: List[Dict[str, Any]] = []
    for uid, rec in idx.by_uid.items():
        if role_u and rec.get("role") != role_u:
            continue
        if ct_sub and ct_sub not in str(rec.get("compound_type") or "").lower():
            continue
        if ql:
            hay = " ".join(
                [
                    str(rec.get("name") or ""),
                    str(rec.get("token") or ""),
                    str(rec.get("abbreviation") or ""),
                    str(rec.get("generic_core") or ""),
                    str(uid or ""),
                ]
            ).lower()
            if ql not in hay:
                continue
        out.append(
            {
                "uid": uid,
                "role": rec.get("role"),
                "name": rec.get("name"),
                "compound_type": rec.get("compound_type"),
            }
        )
        if len(out) >= int(limit):
            break
    return out
