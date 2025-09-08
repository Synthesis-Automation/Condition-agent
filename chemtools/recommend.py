from __future__ import annotations

from typing import Dict, Any, List, Tuple
from collections import Counter

from .smiles import normalize_reaction
from .router import detect_family
from .featurizers import ullmann as feat_ullmann
from . import precedent, constraints, explain


def _pick_electrophile_nucleophile(reactants: List[str]) -> Tuple[str, str]:
    def is_electrophile(s: str) -> bool:
        t = (s or "").lower()
        return (
            ("br" in t) or ("cl" in t) or (" i" in t)
            or ("os(=o)(=o)c(f)(f)f" in t) or ("otf" in t)
        )
    if not reactants:
        return "", ""
    if len(reactants) == 1:
        return reactants[0], ""
    r0, r1 = reactants[0], reactants[1]
    if is_electrophile(r0):
        return r0, r1
    if is_electrophile(r1):
        return r1, r0
    return r0, r1


def _median(vals: List[float]) -> float | None:
    xs = [float(v) for v in vals if isinstance(v, (int, float))]
    if not xs:
        return None
    xs.sort()
    n = len(xs)
    mid = n // 2
    if n % 2 == 1:
        return xs[mid]
    return 0.5 * (xs[mid - 1] + xs[mid])


def _pick_with_constraints(cands: List[str], rules: Dict[str, Any]) -> Tuple[str | None, Dict[str, Any]]:
    if not cands:
        return None, {"allowed": [], "blocked": []}
    if not rules:
        return cands[0], {"allowed": cands, "blocked": []}
    out = constraints.apply_filter(cands, rules)
    allowed = out.get("allowed") or []
    return (allowed[0] if allowed else None), out


def recommend_from_reaction(
    reaction: str,
    k: int = 25,
    relax: Dict[str, Any] | None = None,
    constraint_rules: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Recommend conditions from a reaction SMILES.

    Returns dict with keys: input, family, features, bin, recommendation, alternatives, pack, reasons.
    """
    relax = dict(relax or {})

    # 1) Normalize and extract reactants
    norm = normalize_reaction(reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect family
    fam = detect_family(reactants).get("family") or "Unknown"

    # 3) Featurize substrates (Ullmann featurizer also used as fallback)
    elec, nuc = _pick_electrophile_nucleophile(reactants)
    features: Dict[str, Any] = {}
    if fam == "Ullmann_CN":
        features = feat_ullmann.featurize(elec, nuc)
    else:
        features = feat_ullmann.featurize(elec, nuc)

    # 4) Retrieve precedents (enable DRFP unless explicitly disabled)
    relax.setdefault("reaction_smiles", norm.get("normalized") or reaction)
    relax.setdefault("use_drfp", True)
    relax.setdefault("precompute_drfp", True)
    relax.setdefault("precompute_scope", "candidates")
    pack = precedent.knn(family=fam, features=features, k=int(k), relax=relax)

    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])
    support = int(pack.get("support") or len(precs))

    # 5) Core vote (Laplace smoothing)
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    labels = list(core_counts.keys())
    alpha = 1.0
    denom = sum(core_counts.values()) + alpha * max(1, len(labels))
    scores = {L: (core_counts.get(L, 0) + alpha) / denom for L in labels}
    chosen_core = max(scores, key=scores.get) if scores else None
    core_vote_share = (core_counts.get(chosen_core, 0) / max(1, sum(core_counts.values()))) if chosen_core else 0.0

    # 6) Choose base and solvent (conditioned on chosen core when possible)
    if chosen_core:
        group = [p for p in precs if str(p.get("core") or "") == chosen_core]
    else:
        group = precs
    bases = [str(p.get("base_uid") or "") for p in group if p.get("base_uid")]
    solvents = [str(p.get("solvent_uid") or "") for p in group if p.get("solvent_uid")]
    base_counts = Counter(bases)
    solv_counts = Counter(solvents)
    base_list = [b for b, _ in base_counts.most_common()] or [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
    solv_list = [s for s, _ in solv_counts.most_common()] or [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]

    base_pick, base_filter = _pick_with_constraints(base_list, constraint_rules or {})
    solv_pick, solv_filter = _pick_with_constraints(solv_list, constraint_rules or {})
    # Fallback to overall precedents if filtered out or empty
    if (not base_pick) and precs:
        all_bases = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
        base_pick, base_filter = _pick_with_constraints(list(dict.fromkeys(all_bases)), constraint_rules or {})
    if (not solv_pick) and precs:
        all_solv = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
        solv_pick, solv_filter = _pick_with_constraints(list(dict.fromkeys(all_solv)), constraint_rules or {})

    # 7) T/time median from the same-core group; fallback to all precedents
    def nums(key: str, items: List[Dict[str, Any]]):
        return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]

    T_med = _median(nums("T_C", group) or nums("T_C", precs))
    t_med = _median(nums("time_h", group) or nums("time_h", precs))

    # 8) Confidence (simple): share of top core among all with any core label, capped to [0.3, 0.95]
    conf = 0.95 * core_vote_share if support >= 5 else 0.5 * core_vote_share
    conf = max(0.3, min(0.95, conf))

    # 9) Reasons from precedents
    reasons_pack = explain.for_pack(pack, features)

    recommendation = {
        "core": chosen_core,
        "base_uid": base_pick,
        "solvent_uid": solv_pick,
        "T_C": T_med,
        "time_h": t_med,
        "confidence": round(float(conf), 3),
    }

    alternatives = {
        "cores": core_counts.most_common(3),
        "bases": base_counts.most_common(3),
        "solvents": solv_counts.most_common(3),
    }

    return {
        "input_reaction": reaction,
        "family": fam,
        "features": features,
        "bin": features.get("bin"),
        "recommendation": recommendation,
        "alternatives": alternatives,
        "precedent_pack": pack,
        "reasons": reasons_pack.get("reasons"),
        "filters": {"base": base_filter, "solvent": solv_filter},
    }


def _well_ids(n: int) -> List[str]:
    # Generate common plate IDs. For 24-well, 4 rows (A-D) x 6 columns (1-6).
    # General small helper that creates as close to square as possible grid.
    import math
    rows = int(math.sqrt(n))
    while rows > 1 and n % rows != 0:
        rows -= 1
    if rows <= 1:
        rows = min(8, n)  # cap rows to 8 for small plates
    cols = (n + rows - 1) // rows
    rows = max(1, min(8, rows))
    cols = max(1, cols)
    letters = [chr(ord('A') + i) for i in range(rows)]
    ids: List[str] = []
    for r in letters:
        for c in range(1, cols + 1):
            ids.append(f"{r}{c}")
            if len(ids) >= n:
                return ids
    return ids[:n]


def design_plate_from_reaction(
    reaction: str,
    plate_size: int = 24,
    relax: Dict[str, Any] | None = None,
    constraint_rules: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """Design a diversified plate across cores for a reaction.

    Returns dict with keys: csv (string), rows (list of dict), meta.
    """
    relax = dict(relax or {})

    # Build one precedent pack (larger k for variety)
    rec = recommend_from_reaction(reaction, k=max(plate_size, 50), relax=relax, constraint_rules=constraint_rules or {})
    pack = rec.get("precedent_pack") or {}
    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])

    # Core ranking by frequency
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    core_list = [c for c, _ in core_counts.most_common() if c]
    if not core_list:
        return {"csv": "well_id,core,base_uid,solvent_uid,additive_uids,T_C,time_h\n", "rows": [], "meta": {"error": "NO_PRECEDENTS"}}

    # Build per-core groups for base/solvent extraction
    by_core: Dict[str, List[Dict[str, Any]]] = {}
    for p in precs:
        c = str(p.get("core") or "")
        if not c:
            continue
        by_core.setdefault(c, []).append(p)

    # Ordered core sequence to fill plate, cycling if needed
    seq: List[str] = []
    while len(seq) < int(plate_size):
        for c in core_list:
            seq.append(c)
            if len(seq) >= int(plate_size):
                break

    # For each core, pick base/solvent and T/time from its group
    rows_out: List[Dict[str, Any]] = []
    for i, core in enumerate(seq):
        group = by_core.get(core, [])
        bases = [str(p.get("base_uid") or "") for p in group if p.get("base_uid")]
        solvents = [str(p.get("solvent_uid") or "") for p in group if p.get("solvent_uid")]
        base_counts = Counter(bases)
        solv_counts = Counter(solvents)
        base_list = [b for b, _ in base_counts.most_common()]
        solv_list = [s for s, _ in solv_counts.most_common()]
        b_pick, _bf = _pick_with_constraints(base_list, constraint_rules or {})
        s_pick, _sf = _pick_with_constraints(solv_list, constraint_rules or {})
        # Fallback across all precedents if needed
        if not b_pick:
            all_b = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
            b_pick, _ = _pick_with_constraints(list(dict.fromkeys(all_b)), constraint_rules or {})
        if not s_pick:
            all_s = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
            s_pick, _ = _pick_with_constraints(list(dict.fromkeys(all_s)), constraint_rules or {})
        # T/time medians per core
        def nums(key: str, items: List[Dict[str, Any]]):
            return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]
        T_med = _median(nums("T_C", group) or nums("T_C", precs))
        t_med = _median(nums("time_h", group) or nums("time_h", precs))

        rows_out.append({
            "well_id": _well_ids(int(plate_size))[i],
            "core": core,
            "base_uid": b_pick or "",
            "solvent_uid": s_pick or "",
            "additive_uids": "",
            "T_C": T_med if T_med is not None else "",
            "time_h": t_med if t_med is not None else "",
        })

    # Render CSV
    header = ["well_id", "core", "base_uid", "solvent_uid", "additive_uids", "T_C", "time_h"]
    def _csv_escape(v: Any) -> str:
        s = "" if v is None else str(v)
        return ('"' + s.replace('"', '""') + '"') if ("," in s or '"' in s) else s
    lines = [",".join(header)]
    for row in rows_out:
        lines.append(
            ",".join(_csv_escape(row[k]) for k in header)
        )
    csv_text = "\n".join(lines) + "\n"

    meta = {
        "family": rec.get("family"),
        "bin": rec.get("bin"),
        "cores": core_counts.most_common(8),
        "precedent_support": int(pack.get("support") or len(precs)),
    }
    return {"csv": csv_text, "rows": rows_out, "meta": meta}
