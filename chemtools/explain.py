from __future__ import annotations

from typing import Dict, Any, List, Tuple
from collections import Counter


def _fmt_bin(features: Dict[str, Any]) -> str:
    bin_str = features.get("bin")
    if not bin_str:
        lg = features.get("LG", "?")
        nuc = features.get("nuc_class", "?")
        bin_str = f"LG:{lg}|NUC:{nuc}"
    extras: List[str] = []
    # Add simple boolean/ordinal decorations when present
    if features.get("ortho_count") is not None:
        try:
            oc = int(features["ortho_count"])  # type: ignore[index]
            extras.append(f"ortho={oc}")
        except Exception:
            pass
    if isinstance(features.get("para_EWG"), bool):
        extras.append("para_EWG" if features.get("para_EWG") else "para_no_EWG")
    if isinstance(features.get("heteroaryl"), bool) and features.get("heteroaryl"):
        extras.append("heteroaryl")
    return bin_str + ("; " + ", ".join(extras) if extras else "")


def _top_k_counts(items: List[str], k: int = 2) -> List[Tuple[str, int]]:
    c = Counter([str(x) for x in items if x])
    return c.most_common(k)


def _name_from_uid(uid: str) -> str:
    """Best-effort human-friendly token or fall back to uid."""
    if not uid:
        return ""
    try:
        from .properties import lookup

        r = lookup(uid)
        if isinstance(r, dict) and r.get("found") and isinstance(r.get("record"), dict):
            rec = r["record"]
            token = rec.get("token") or rec.get("name") or rec.get("uid")
            return str(token)
    except Exception:
        pass
    return str(uid)


def _fmt_list(entries: List[Tuple[str, int]]) -> str:
    parts: List[str] = []
    for val, cnt in entries:
        name = _name_from_uid(val) if isinstance(val, str) else str(val)
        parts.append(f"{name} ({cnt})")
    return ", ".join(parts)


def _summarize_numeric(values: List[float]) -> str:
    vals = [float(v) for v in values if isinstance(v, (int, float))]
    if not vals:
        return "unknown"
    vals.sort()
    n = len(vals)
    mid = vals[n // 2]
    if n >= 3:
        lo = vals[max(0, n // 4)]
        hi = vals[min(n - 1, (3 * n) // 4)]
        if lo == hi:
            return f"~{int(round(mid))} C" if abs(mid - int(round(mid))) < 1e-9 else f"~{mid:.1f}"
        return f"{int(round(lo))}-{int(round(hi))}"
    return f"~{int(round(mid))}"


def for_pack(pack: Dict[str, Any], features: Dict[str, Any]) -> Dict[str, Any]:
    """Generate short, human-readable reasons that explain the precedent pack.

    Inputs:
    - pack: output from precedent.knn, expected to contain {precedents:[{core, base_uid, solvent_uid, T_C, time_h, ...}], support}
    - features: substrate/features dict including at least bin/LG/nuc_class

    Returns {reasons: [str], precedents: top 3 precedents from pack}
    """
    precedents = list(pack.get("precedents", []) or [])
    support = int(pack.get("support") or len(precedents))

    reasons: List[str] = []
    reasons.append(f"Substrate bin: {_fmt_bin(features)}.")

    if precedents:
        cores = _top_k_counts([str(p.get("core") or "") for p in precedents])
        bases = _top_k_counts([str(p.get("base_uid") or "") for p in precedents])
        solvents = _top_k_counts([str(p.get("solvent_uid") or "") for p in precedents])
        tcs = [p.get("T_C") for p in precedents if p.get("T_C") is not None]
        ths = [p.get("time_h") for p in precedents if p.get("time_h") is not None]

        if cores:
            reasons.append(f"Core(s): {_fmt_list(cores)}.")
        if bases:
            reasons.append(f"Common base(s): {_fmt_list(bases)}.")
        if solvents:
            reasons.append(f"Common solvent(s): {_fmt_list(solvents)}.")
        if tcs:
            reasons.append(f"Typical T (C): {_summarize_numeric(tcs)}.")
        if ths:
            reasons.append(f"Typical time (h): {_summarize_numeric(ths)}.")

    reasons.append(f"Support: {support} precedent(s).")

    return {"reasons": reasons, "precedents": precedents[:3]}
