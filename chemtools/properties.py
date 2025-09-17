from __future__ import annotations
"""
chemtools.properties

Minimal properties/role lookup backed by the compound taxonomy (data/compound_taxonomy).
Records include uid (CAS), role, token, optional name/abbreviation, compound_type, and generic_core for catalyst precursors.

- Default source: taxonomy JSON files
- Override taxonomy dir: CHEMTOOLS_TAXONOMY_DIR
- Optional extra overrides: CHEMTOOLS_PROPERTIES_PATH (JSON/JSONL)

No default dependency on a merged CAS registry; lookups prefer taxonomy and only optionally fall back to registry.resolve.
"""

from typing import Dict, Any, Optional, Tuple, List
import json
import os
import re


# Built-in seed properties. These are minimal and can be overridden/extended
# via an external JSON/JSONL pointed to by CHEMTOOLS_PROPERTIES_PATH.
_SEED: Dict[str, Dict[str, Any]] = {
    # Bases
    "7778-53-2": {"role": "BASE", "token": "K3PO4", "pKa_DMSO": 30.0},
    "1310-58-3": {"role": "BASE", "token": "KOH", "pKa_water": 15.7},
    # Solvents
    "7732-18-5": {"role": "SOLVENT", "token": "Water", "KT": {"alpha": 1.17, "beta": 0.47, "pi*": 1.09}},
    # Catalysts/cores
    "7681-65-4": {"role": "CATALYST", "token": "CuI"},
    # Ligands
    "72-52-8": {"role": "LIGAND", "token": "Phenanthroline"},
}


_CACHE: Optional[Dict[str, Dict[str, Any]]] = None
_ALIAS_TO_UID: Optional[Dict[str, str]] = None


def _norm_alias(s: str) -> str:
    s0 = (s or "").strip().lower()
    if not s0:
        return ""
    # unify dashes and basic punctuation; keep alnum only for robust matching
    s0 = s0.replace("—", "-").replace("–", "-")
    return re.sub(r"[^a-z0-9]+", "", s0)


def _taxonomy_dir() -> Optional[str]:
    # Allow override via env
    env = os.environ.get("CHEMTOOLS_TAXONOMY_DIR")
    if env and os.path.isdir(env):
        return env
    try:
        base = os.path.dirname(os.path.dirname(__file__))
        d = os.path.join(base, "data", "compound_taxonomy")
        return d if os.path.isdir(d) else None
    except Exception:
        return None


def _load_taxonomy_props() -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    """Build a properties index from compound taxonomy files.

    Returns (props_by_uid, alias_to_uid), where uid is CAS when available.
    """
    props: Dict[str, Dict[str, Any]] = {}
    alias_to_uid: Dict[str, str] = {}
    d = _taxonomy_dir()
    if not d:
        return props, alias_to_uid

    def add_alias(uid: str, text: Optional[str]) -> None:
        na = _norm_alias(text or "")
        if na:
            alias_to_uid.setdefault(na, uid)

    def add_member(role: str, compound_type: str, fam: Dict[str, Any], em: Dict[str, Any], extra: Dict[str, Any] | None = None) -> None:
        uid = (em.get("cas") or "").strip()
        if not uid:
            # skip entries without CAS for uid; can still be resolved by alias later
            return
        name = (em.get("name") or "").strip()
        abbr = (em.get("abbr") or "").strip()
        token = (abbr or name or uid).strip()
        rec: Dict[str, Any] = {
            "uid": uid,
            "role": role,
            "compound_type": compound_type,
            "name": name or uid,
            "abbreviation": abbr or None,
            "token": token,
        }
        if extra:
            rec.update({k: v for k, v in extra.items() if v is not None})
        # Merge numeric_features minimally if present
        if isinstance(em.get("numeric_features"), dict):
            rec.update({f"num__{k}": v for k, v in em["numeric_features"].items()})
        props[uid] = {**props.get(uid, {}), **rec}
        # Aliases
        add_alias(uid, uid)
        add_alias(uid, name)
        add_alias(uid, abbr)
        for s in (em.get("synonyms") or []):
            add_alias(uid, s)
        # Also index family label/id as weak aliases
        add_alias(uid, fam.get("label"))
        add_alias(uid, fam.get("family_id"))

    def load_json(path: str) -> Dict[str, Any] | None:
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return None

    # Ligands
    lig = load_json(os.path.join(d, "taxonomy_ligand.json"))
    if lig:
        for fam in lig.get("families", []) or []:
            for em in fam.get("example_members", []) or []:
                add_member("LIGAND", "ligand", fam, em, extra={})

    # Bases
    bas = load_json(os.path.join(d, "taxonomy_base.json"))
    if bas:
        for fam in bas.get("families", []) or []:
            for em in fam.get("example_members", []) or []:
                add_member("BASE", "base", fam, em, extra={})

    # Solvents
    solv = load_json(os.path.join(d, "taxonomy_solvent.json"))
    if solv:
        for fam in solv.get("families", []) or []:
            for em in fam.get("example_members", []) or []:
                add_member("SOLVENT", "solvent", fam, em, extra={})

    # Coupling reagents (treat as ADDITIVE for role)
    cou = load_json(os.path.join(d, "taxonomy_coupling_reagent.json"))
    if cou:
        for fam in cou.get("families", []) or []:
            for em in fam.get("example_members", []) or []:
                add_member("ADDITIVE", "coupling_reagent", fam, em, extra={})

    # Reductants (treated as ADDITIVE role)
    red = load_json(os.path.join(d, "taxonomy_reductant.json"))
    if red:
        for fam in red.get("families", []) or []:
            extra: Dict[str, Any] = {}
            fam_class = fam.get("class")
            if fam_class:
                extra["reductant_class"] = fam_class
            mechanism = fam.get("mechanism")
            if mechanism:
                extra["mechanism"] = mechanism
            strength = fam.get("strength_index")
            if strength is not None:
                extra["strength_index"] = strength
            recommended = fam.get("recommended_pairs")
            if recommended:
                extra["family_recommended_pairs"] = recommended
            for em in fam.get("example_members", []) or []:
                uid = (em.get("cas") or "").strip()
                add_member("ADDITIVE", "reductant", fam, em, extra=extra)
                if not uid or uid not in props:
                    continue
                entry = props[uid]
                if em.get("properties"):
                    entry["properties"] = em["properties"]
                if em.get("recommended_pairs"):
                    entry["recommended_pairs"] = em["recommended_pairs"]
                if em.get("embedding_text") and "embedding_text" not in entry:
                    entry["embedding_text"] = em["embedding_text"]

    # Catalyst precursors (metal) → provide generic_core from family metal
    cat = load_json(os.path.join(d, "taxonomy_catalysts_precursor.json"))
    if cat:
        for fam in cat.get("families", []) or []:
            metal = (fam.get("metal") or "").strip()
            for em in fam.get("example_members", []) or []:
                add_member("CATALYST", "catalyst_core", fam, em, extra={"generic_core": metal})

    return props, alias_to_uid


def _load_external() -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    """Load properties from taxonomy dir first. If CHEMTOOLS_PROPERTIES_PATH points
    to a JSON/JSONL, merge that as overrides. Avoid cas_registry_merged.jsonl by default.
    """
    props, alias = _load_taxonomy_props()
    # Optional additional overrides via CHEMTOOLS_PROPERTIES_PATH (JSON/JSONL)
    path = os.environ.get("CHEMTOOLS_PROPERTIES_PATH")
    if path and os.path.exists(path):
        try:
            if path.lower().endswith(".json"):
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    if isinstance(data, dict):
                        for k, v in data.items():
                            if isinstance(v, dict):
                                uid = str(k)
                                props[uid] = {**props.get(uid, {}), **v}
                                # Extend alias by token/name if provided
                                for ak in [uid, v.get("name"), v.get("abbreviation"), v.get("token")]:
                                    na = _norm_alias(ak or "")
                                    if na:
                                        alias.setdefault(na, uid)
            else:
                with open(path, "r", encoding="utf-8") as f:
                    for line in f:
                        s = line.strip()
                        if not s:
                            continue
                        try:
                            obj = json.loads(s)
                        except Exception:
                            continue
                        uid = (obj.get("uid") or obj.get("cas") or "").strip()
                        if not uid:
                            continue
                        rec = dict(obj)
                        props[uid] = {**props.get(uid, {}), **rec}
                        for ak in [uid, rec.get("name"), rec.get("abbreviation"), rec.get("token")]:
                            na = _norm_alias(ak or "")
                            if na:
                                alias.setdefault(na, uid)
        except Exception:
            pass
    return props, alias


def _props() -> Dict[str, Dict[str, Any]]:
    global _CACHE
    global _ALIAS_TO_UID
    if _CACHE is not None:
        return _CACHE
    merged = dict(_SEED)
    ext_props, alias = _load_external()
    for uid, rec in ext_props.items():
        merged[uid] = {**merged.get(uid, {}), **rec}
    _CACHE = merged
    _ALIAS_TO_UID = alias
    return merged


def lookup(query: str, *, allow_registry: bool = True) -> Dict[str, Any]:
    """Lookup properties by CAS or alias.

    - allow_registry: when True, may consult `chemtools.registry.resolve` to
      map names/aliases to a CAS before checking properties. Set to False to
      avoid recursion when `registry.resolve` is enriching results.

    Returns {found, record}. When found, record contains at least {uid, role?, token?, name?, ...props}
    """
    q = (query or "").strip()
    if not q:
        return {"found": False, "record": None}

    props = _props()
    # Direct CAS hit
    if q in props:
        rec = {"uid": q, **props[q]}
        # Best-effort role from compound_type if missing
        if "role" not in rec and rec.get("compound_type"):
            ct = str(rec.get("compound_type") or "").strip().lower()
            role_map = {
                "ligand": "LIGAND",
                "base": "BASE",
                "solvent": "SOLVENT",
                "catalyst_core": "CATALYST",
                "metal": "CATALYST",
                "catalyst": "CATALYST",
                # default other types to ADDITIVE
            }
            rec["role"] = role_map.get(ct, "ADDITIVE")
        return {"found": True, "record": rec}

    # Try case-insensitive match on token/name via registry, unless disabled to avoid recursion
    # First try taxonomy alias map
    global _ALIAS_TO_UID
    if _ALIAS_TO_UID is None:
        _props()  # initializes alias map
    na = _norm_alias(q)
    if na and _ALIAS_TO_UID and na in _ALIAS_TO_UID:
        uid = _ALIAS_TO_UID[na]
        base = props.get(uid, {})
        if base:
            rec = {"uid": uid, **base}
            return {"found": True, "record": rec}

    # Fallback to registry for non-taxonomy names if allowed
    if allow_registry:
        try:
            from .registry import resolve as registry_resolve  # lazy import to avoid cycles at import time
            res = registry_resolve(q)
            if isinstance(res, dict) and res.get("uid") and res.get("uid") in props:
                uid = str(res["uid"])
                base = props[uid]
                rec = {"uid": uid, **base}
                if res.get("name"):
                    rec.setdefault("name", res["name"])  # keep ext/seed if explicitly provided
                if res.get("role"):
                    rec.setdefault("role", res["role"])  # do not overwrite
                return {"found": True, "record": rec}
        except Exception:
            pass

    # Fallback: token/name match across seeds
    ql = q.lower()
    for uid, rec in props.items():
        tok = str(rec.get("token") or "").lower()
        nm = str(rec.get("name") or "").lower()
        if ql == tok or ql == nm:
            return {"found": True, "record": {"uid": uid, **rec}}

    return {"found": False, "record": None}
