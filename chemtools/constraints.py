from typing import Dict, Any, List, Tuple


"""
Constraints filter

Deterministic filters to enforce inventory, environmental and simple operational rules
over a list of candidate reagent/solvent IDs (CAS or tokens).

Inputs:
- candidates: list of strings (UID/CAS or token/name)
- rules: dict with optional keys:
  - inventory: [uids] — whitelist. If present, anything not in it is blocked.
  - blacklist: [uids] — explicit blocklist.
  - no_HMPA: bool — block Hexamethylphosphoramide (CAS 680-31-9, token contains "HMPA").
  - no_chlorinated: bool — block common chlorinated solvents (DCM, chloroform, chlorobenzene, etc.).
  - aqueous_only: bool — allow only water as solvent (non-solvent items pass).
  - min_bp_C: number — require solvent boiling point >= threshold when known.
  - allow_unknown: bool (default True) — unknown properties do not block unless inventory/blacklist says so.

Returns: {allowed: [ids], blocked: [{id, reason}]}
"""


# Small built‑in knowledge for environmental flags (kept tiny; extend via properties DB)
_HMPA_CAS = {"680-31-9"}
_CHLORINATED_SOLVENT_CAS = {
    "75-09-2",   # dichloromethane (DCM)
    "67-66-3",   # chloroform
    "108-90-7",  # chlorobenzene
    "127-18-4",  # tetrachloroethylene
    "79-01-6",   # trichloroethylene
}


def _normalize_id(x: str) -> str:
    return str(x or "").strip().lower()


def _lookup(uid_or_token: str) -> Dict[str, Any]:
    try:
        from .properties import lookup  # lazy import to avoid cycles
        res = lookup(uid_or_token)
        if isinstance(res, dict) and res.get("found") and isinstance(res.get("record"), dict):
            return res["record"]  # contains at least {uid?, role?, token?, name?}
    except Exception:
        pass
    # Fallback minimal record with just the provided id/token
    return {"uid": uid_or_token}


def _is_hmpa(rec: Dict[str, Any]) -> bool:
    uid = str(rec.get("uid") or "").strip()
    tok = str(rec.get("token") or rec.get("name") or rec.get("uid") or "").lower()
    return uid in _HMPA_CAS or "hmpa" in tok


def _is_chlorinated_solvent(rec: Dict[str, Any]) -> bool:
    role = str(rec.get("role") or "").upper()
    uid = str(rec.get("uid") or "").strip()
    tok = str(rec.get("token") or rec.get("name") or rec.get("uid") or "").lower()
    if uid in _CHLORINATED_SOLVENT_CAS:
        return True
    if role == "SOLVENT":
        # Heuristic match on common names
        return any(s in tok for s in ["chloroform", "dichloromethane", "methylene chloride", "chlorobenzene"]) or (
            tok.startswith("chloro") or " dichloro" in tok or " trichloro" in tok
        )
    return False


def _is_water(rec: Dict[str, Any]) -> bool:
    tok = str(rec.get("token") or rec.get("name") or "").lower()
    uid = str(rec.get("uid") or "").strip()
    return tok in {"water", "h2o"} or uid == "7732-18-5"


def _solvent_bp(rec: Dict[str, Any]) -> float | None:
    v = rec.get("bp_C") or rec.get("bp_c")
    try:
        return float(v) if v is not None else None
    except Exception:
        return None


def apply_filter(candidates: List[str], rules: Dict[str, Any]) -> Dict[str, Any]:
    rules = rules or {}
    inventory = {_normalize_id(x) for x in rules.get("inventory", [])}
    blacklist = {_normalize_id(x) for x in rules.get("blacklist", [])}
    no_hmpa = bool(rules.get("no_HMPA") or rules.get("no_hmpa"))
    no_chloro = bool(rules.get("no_chlorinated"))
    aqueous_only = bool(rules.get("aqueous_only"))
    min_bp = rules.get("min_bp_C")
    allow_unknown = True if rules.get("allow_unknown") is None else bool(rules.get("allow_unknown"))

    allowed: List[str] = []
    blocked: List[Dict[str, str]] = []

    for c in candidates or []:
        cid_norm = _normalize_id(c)

        # 1) Explicit blacklist
        if blacklist and cid_norm in blacklist:
            blocked.append({"id": c, "reason": "blacklisted"})
            continue

        # 2) Inventory whitelist
        if inventory and cid_norm not in inventory:
            blocked.append({"id": c, "reason": "not_in_inventory"})
            continue

        # Properties/heuristics
        rec = _lookup(c)

        # 3) Environmental rules
        if no_hmpa and _is_hmpa(rec):
            blocked.append({"id": c, "reason": "banned_hmpa"})
            continue
        if no_chloro and _is_chlorinated_solvent(rec):
            blocked.append({"id": c, "reason": "chlorinated_solvent"})
            continue

        # 4) Aqueous-only rule (solvents only)
        if aqueous_only and str(rec.get("role") or "").upper() == "SOLVENT" and not _is_water(rec):
            blocked.append({"id": c, "reason": "non_aqueous_solvent"})
            continue

        # 5) Boiling point threshold (only when known)
        if min_bp is not None and str(rec.get("role") or "").upper() == "SOLVENT":
            bp = _solvent_bp(rec)
            if bp is None and not allow_unknown:
                blocked.append({"id": c, "reason": "bp_unknown"})
                continue
            if bp is not None and bp < float(min_bp):
                blocked.append({"id": c, "reason": "bp_below_min"})
                continue

        # Pass
        allowed.append(c)

    return {"allowed": allowed, "blocked": blocked}
