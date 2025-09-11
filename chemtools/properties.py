from __future__ import annotations

from typing import Dict, Any, Optional
import json
import os


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


def _load_external() -> Dict[str, Dict[str, Any]]:
    # Prefer explicit override, else fall back to bundled registry JSONL
    path = os.environ.get("CHEMTOOLS_PROPERTIES_PATH")
    if not path:
        # Default to the same merged CAS registry used by chemtools.registry
        try:
            base = os.path.dirname(os.path.dirname(__file__))
            path = os.path.join(base, "data", "cas_registry_merged.jsonl")
        except Exception:
            path = None  # type: ignore
    if not path:
        return {}
    if not os.path.exists(path):
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    try:
        if path.lower().endswith(".json"):
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
                if isinstance(data, dict):
                    # Expect {cas: {props}}
                    for k, v in data.items():
                        if isinstance(v, dict):
                            out[str(k)] = v
        else:
            # JSONL: one object per line with at least {uid/cas, ...}
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    s = line.strip()
                    if not s:
                        continue
                    try:
                        obj = json.loads(s)
                    except Exception:
                        continue
                    uid = obj.get("uid") or obj.get("cas")
                    if not uid:
                        continue
                    rec = dict(obj)
                    out[str(uid)] = rec
    except Exception:
        # Be tolerant of malformed external data
        return {}
    return out


def _props() -> Dict[str, Dict[str, Any]]:
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    merged = dict(_SEED)
    ext = _load_external()
    # Merge with external, allowing overrides
    for uid, rec in ext.items():
        merged[uid] = {**merged.get(uid, {}), **rec}
    _CACHE = merged
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
    res = None
    if allow_registry:
        try:
            from .registry import resolve as registry_resolve  # lazy import to avoid cycles at import time
            res = registry_resolve(q)
        except Exception:
            res = None

    if isinstance(res, dict) and res.get("uid") and res.get("uid") in props:
        uid = str(res["uid"])
        base = props[uid]
        rec = {"uid": uid, **base}
        # Prefer token/name from registry when available
        if res.get("name"):
            rec.setdefault("name", res["name"])  # keep ext/seed value if explicitly provided
        if res.get("role"):
            rec.setdefault("role", res["role"])  # seed may already have role
        return {"found": True, "record": rec}

    # Fallback: token match across seeds
    ql = q.lower()
    for uid, rec in props.items():
        tok = str(rec.get("token") or "").lower()
        nm = str(rec.get("name") or "").lower()
        if ql == tok or ql == nm:
            return {"found": True, "record": {"uid": uid, **rec}}

    return {"found": False, "record": None}
