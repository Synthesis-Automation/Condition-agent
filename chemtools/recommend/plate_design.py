"""
Experimental plate design generation.

Creates diversified plate layouts from reaction recommendations.
"""

from __future__ import annotations

from typing import Dict, Any, List
from collections import Counter
import math

from .utils import median, pick_with_constraints


def _well_ids(n: int) -> List[str]:
    """
    Generate well IDs for plate (e.g., A1, A2, ..., B1, B2, ...).
    
    For 24-well: 4 rows (A-D) x 6 columns (1-6)
    Creates as close to square grid as possible.
    
    Args:
        n: Number of wells
        
    Returns:
        List of well IDs
    """
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
    """
    Design a diversified experimental plate from reaction.
    
    Distributes conditions across different catalyst cores with
    most common bases/solvents for each core.
    
    Args:
        reaction: Reaction SMILES
        plate_size: Number of wells (default: 24)
        relax: Relaxation rules for precedent search
        constraint_rules: Constraint rules for reagent selection
        
    Returns:
        Dict with:
            - csv: CSV string of plate layout
            - rows: List of row dicts (well_id, core, base_uid, solvent_uid, T_C, time_h)
            - meta: Metadata (family, bin, cores, precedent_support)
    """
    # Import here to avoid circular dependency
    from .core import recommend_from_reaction as _recommend_from_reaction
    
    relax = dict(relax or {})

    # Build one precedent pack (larger k for variety)
    rec = _recommend_from_reaction(
        reaction,
        k=max(plate_size, 50),
        relax=relax,
        constraint_rules=constraint_rules or {}
    )
    pack = rec.get("precedent_pack") or {}
    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])

    # Core ranking by frequency
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    core_list = [c for c, _ in core_counts.most_common() if c]
    
    if not core_list:
        return {
            "csv": "well_id,core,base_uid,solvent_uid,additive_uids,T_C,time_h\n",
            "rows": [],
            "meta": {"error": "NO_PRECEDENTS"}
        }

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
        
        b_pick, _bf = pick_with_constraints(base_list, constraint_rules or {})
        s_pick, _sf = pick_with_constraints(solv_list, constraint_rules or {})
        
        # Fallback across all precedents if needed
        if not b_pick:
            all_b = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
            b_pick, _ = pick_with_constraints(list(dict.fromkeys(all_b)), constraint_rules or {})
        if not s_pick:
            all_s = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
            s_pick, _ = pick_with_constraints(list(dict.fromkeys(all_s)), constraint_rules or {})
        
        # T/time medians per core
        def nums(key: str, items: List[Dict[str, Any]]):
            return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]
        
        T_med = median(nums("T_C", group) or nums("T_C", precs))
        t_med = median(nums("time_h", group) or nums("time_h", precs))

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
        lines.append(",".join(_csv_escape(row[k]) for k in header))
    csv_text = "\n".join(lines) + "\n"

    meta = {
        "family": rec.get("family"),
        "bin": rec.get("bin"),
        "cores": core_counts.most_common(8),
        "precedent_support": int(pack.get("support") or len(precs)),
    }
    
    return {"csv": csv_text, "rows": rows_out, "meta": meta}
