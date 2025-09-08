from __future__ import annotations

from typing import Dict, Any, List, Optional, Tuple
import os

try:  # avoid hard dependency at import time
    from .contracts import Reagent  # type: ignore
except Exception:  # pragma: no cover
    Reagent = Any  # type: ignore

from .registry import resolve as registry_resolve


def _text_norm(s: str) -> str:
    return (s or "").strip().lower()


def _best_str(*vals: Optional[str]) -> str:
    for v in vals:
        if v:
            return str(v)
    return ""


def _as_bool_env(name: str, default: bool = False) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    v = v.strip().lower()
    return v in ("1", "true", "yes", "on")


def _read_dataset_aliases() -> Tuple[Dict[str, str], Dict[str, str]]:
    """Derive ligand alias maps from Ullmann dataset, if present.

    Returns (by_cas, by_name_norm) where values are canonical ligand names.
    """
    import json, os

    # Allow tests or fast-start modes to skip dataset scanning entirely
    if _as_bool_env("CHEMTOOLS_SKIP_DATASET_ALIASES") or _as_bool_env("CHEMTOOLS_FAST_TEST"):
        return {}, {}

    by_cas: Dict[str, str] = {}
    by_name: Dict[str, str] = {}
    root = os.path.dirname(os.path.dirname(__file__))
    # Allow overriding dataset path via env for flexibility; accept file or directory and scan *.jsonl
    env_path = os.environ.get("CHEMTOOLS_LIGAND_ALIAS_PATH")
    paths: List[str] = []
    if env_path:
        if os.path.isdir(env_path):
            for fn in os.listdir(env_path):
                if fn.lower().endswith(".jsonl"):
                    paths.append(os.path.join(env_path, fn))
        elif os.path.isfile(env_path):
            paths.append(env_path)
    else:
        ds_dir = os.path.join(root, "data", "reaction_dataset")
        if os.path.isdir(ds_dir):
            for fn in os.listdir(ds_dir):
                if fn.lower().endswith(".jsonl"):
                    paths.append(os.path.join(ds_dir, fn))
    if not paths:
        return by_cas, by_name
    try:
        for path in paths:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except Exception:
                        continue
                    core = rec.get("condition_core") or ""
                    if "/" not in core:
                        continue
                    metal, ligand = core.split("/", 1)
                    ligand = (ligand or "").strip()
                    if not ligand or ligand.lower() == "none":
                        continue
                    # Normalize canonical ligand token capitalization
                    canonical = ligand
                    canonical_norm = _text_norm(canonical)
                    metal_norm = _text_norm(metal)

                    def _is_metal_name(n: str) -> bool:
                        t = _text_norm(n)
                        if not t:
                            return False
                        return (
                            (metal_norm and metal_norm in t)
                            or ("copper" in t)
                            or ("palladium" in t)
                            or ("nickel" in t)
                            or ("ruthenium" in t)
                        )

                    # Gather aliases from catalyst.core and enriched names if present
                    cat = rec.get("catalyst") or {}
                    for lst_key in ("core", "full_system"):
                        for item in (cat.get(lst_key) or []):
                            cas = (item.get("cas") or "").strip()
                            nm = (item.get("name") or "").strip()
                            nm_norm = _text_norm(nm)
                            looks_like_ligand = (
                                canonical_norm
                                and (nm_norm == canonical_norm or canonical_norm in nm_norm or nm_norm in canonical_norm)
                                and not _is_metal_name(nm)
                            )
                            if looks_like_ligand:
                                if cas:
                                    by_cas.setdefault(cas, canonical)
                                if nm:
                                    by_name.setdefault(nm_norm, canonical)
                    enr = (rec.get("raw_data") or {}).get("enriched_names") or {}
                    for nm in (enr.get("catalysts") or []):
                        if not _is_metal_name(str(nm)):
                            by_name.setdefault(_text_norm(str(nm)), canonical)
    except Exception:
        # Be resilient if dataset is malformed
        pass
    # Seed some common aliases for safety
    by_name.setdefault("phen", "Phenanthroline")
    by_name.setdefault("1,10-phenanthroline", "Phenanthroline")
    by_name.setdefault("dmeda", "DMEDA")
    by_name.setdefault("tmeda", "TMEDA")
    return by_cas, by_name


_LIG_BY_CAS, _LIG_BY_NAME = _read_dataset_aliases()


def _get_attr(r: Any, key: str) -> str:
    if hasattr(r, key):
        v = getattr(r, key)
        return v if isinstance(v, str) else str(v) if v is not None else ""
    if isinstance(r, dict):
        v = r.get(key)
        return v if isinstance(v, str) else str(v) if v is not None else ""
    return ""


def _pick_metal_symbol(name: str, uid: str) -> Optional[str]:
    # Prefer registry generic_core if available
    res = registry_resolve(uid or name)
    gen = (res.get("record") or {}).get("generic_core") if isinstance(res, dict) else None
    if not gen:
        gen = (res or {}).get("generic_core") if isinstance(res, dict) else None
    if gen:
        return str(gen)
    t = _text_norm(name + " " + uid)
    if "cu" in t or "copper" in t:
        return "Cu"
    if "pd" in t or "palladium" in t:
        return "Pd"
    if "ni" in t or "nickel" in t:
        return "Ni"
    if "ru" in t or "ruthenium" in t:
        return "Ru"
    return None


def _pick_ligand_canonical(name: str, uid: str) -> Optional[str]:
    # Dataset-derived CAS mapping
    if uid and uid in _LIG_BY_CAS:
        return _LIG_BY_CAS[uid]
    # Registry token/name hints
    res = registry_resolve(uid or name)
    if isinstance(res, dict):
        nm = res.get("name") or ""
        tok = (res.get("record") or {}).get("token") if res.get("record") else None
        if tok:
            t = _text_norm(tok)
            if t in _LIG_BY_NAME:
                return _LIG_BY_NAME[t]
            # Normalize common ligands
            if "phen" in t:
                return "Phenanthroline"
        if nm:
            t = _text_norm(nm)
            if t in _LIG_BY_NAME:
                return _LIG_BY_NAME[t]
    # Text-based fallback
    t = _text_norm(name)
    for key, canon in _LIG_BY_NAME.items():
        if key and key in t:
            return canon
    if "xphos" in t:
        return "XPhos"
    if "dmeda" in t:
        return "DMEDA"
    if "tmeda" in t:
        return "TMEDA"
    if "phen" in t:
        return "Phenanthroline"
    return None


def parse(reagents: List[Reagent], text: str | None = None) -> Dict[str, Any]:
    """Parse reagents into a normalized condition core.

    Returns: {core, metal_source_uid, ligand_uid?, ratio?:float, precatalyst?:bool}
    """
    metal: Optional[str] = None
    ligand: Optional[str] = None
    metal_uid: Optional[str] = None
    ligand_uid: Optional[str] = None
    precatalyst = False

    tnorm = _text_norm(text or "")

    # First pass: explicit roles
    for r in reagents:
        role = (_get_attr(r, "role") or "").upper()
        uid = _get_attr(r, "uid")
        nm = _best_str(_get_attr(r, "name"), _get_attr(r, "token"), uid)
        if role in {"CATALYST", "METAL", "METAL_SOURCE", "CATALYST_CORE"}:
            sym = _pick_metal_symbol(nm, uid)
            if sym and not metal:
                metal = sym
                metal_uid = uid
            # pre-catalyst detection via registry compound_type
            res = registry_resolve(uid or nm)
            ctype = None
            if isinstance(res, dict):
                ctype = (res.get("record") or {}).get("compound_type") or res.get("compound_type")
            if ctype and isinstance(ctype, str) and "preformed" in ctype.lower():
                precatalyst = True
            # Some datasets list the ligand under catalyst/core; attempt to pick ligand too
            if not ligand:
                lig = _pick_ligand_canonical(nm, uid)
                if lig:
                    ligand = lig
                    ligand_uid = uid
        if role == "LIGAND" and not ligand:
            lig = _pick_ligand_canonical(nm, uid)
            if lig:
                ligand = lig
                ligand_uid = uid

    # Second pass: if ligand missing, look at additives possibly acting as ligands
    if not ligand:
        for r in reagents:
            role = (_get_attr(r, "role") or "").upper()
            if role in {"ADDITIVE", "ACTIVATOR", "OXIDANT", "BASE"}:
                uid = _get_attr(r, "uid")
                nm = _best_str(_get_attr(r, "name"), _get_attr(r, "token"), uid)
                lig = _pick_ligand_canonical(nm, uid)
                if lig:
                    ligand = lig
                    ligand_uid = uid
                    break

    # If metal still missing, infer from text
    if not metal and tnorm:
        if "copper" in tnorm or "cui" in tnorm or "cu(" in tnorm or " cu " in tnorm:
            metal = "Cu"
        elif "palladium" in tnorm or "pd(" in tnorm or " pd " in tnorm:
            metal = "Pd"

    # If ligand still missing, infer from text
    if not ligand and tnorm:
        ligand = _pick_ligand_canonical(tnorm, "")

    # Compose core string
    if metal and ligand:
        core = f"{metal}/{ligand}"
    elif metal:
        core = metal
    else:
        core = "Unknown"

    out: Dict[str, Any] = {
        "core": core,
        "metal_source_uid": metal_uid,
        "ligand_uid": ligand_uid,
        "precatalyst": precatalyst,
    }
    return out


# Backward compatibility
def parse_core(reagents: List[Reagent], text: str = "") -> Dict[str, Any]:
    return parse(reagents, text)
