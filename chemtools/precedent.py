from typing import Dict, Any, List, Tuple
import os, json
from functools import lru_cache

# Path to small demo dataset used for precedents retrieval
DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "reactions_sample.jsonl")


@lru_cache(maxsize=1)
def _load() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if os.path.exists(DATA_PATH):
        with open(DATA_PATH, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except Exception:
                    # Skip malformed lines in demo data
                    pass
    return rows


def _family_text(family: str) -> str:
    # Map API family tokens to dataset labels
    f = (family or "").strip()
    if f.lower() in {"ullmann_cn", "ullmann c–n", "ullmann c-n", "ullmann"}:
        return "Ullmann C–N"
    return f


def _parse_bin(bin_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not bin_str:
        return out
    for part in str(bin_str).split("|"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def _candidate_pool(rows: List[Dict[str, Any]], family_txt: str, feat: Dict[str, Any], k: int, relax: Dict[str, Any]) -> List[Dict[str, Any]]:
    # Filter rows to family first
    fam_rows = [r for r in rows if (r.get("rxn_type") or "") == family_txt]
    if not fam_rows:
        return []

    strict_bin = relax.get("strict_bin", True)
    min_candidates = int(relax.get("min_candidates", k))
    fallback_order: List[str] = relax.get("fallback_order", ["nuc_class", "LG", "any"])  # type: ignore

    target_bin = (feat.get("bin") or "").strip()
    target_bin_map = _parse_bin(target_bin)
    target_nuc = (feat.get("nuc_class") or target_bin_map.get("NUC") or "").lower()
    target_lg = feat.get("LG") or target_bin_map.get("LG") or ""

    # Exact bin matches
    cands = [r for r in fam_rows if (r.get("features", {}).get("bin") or "") == target_bin]
    if len(cands) >= min_candidates or strict_bin:
        return cands

    # Fallbacks
    remaining = [r for r in fam_rows if r not in cands]
    for fb in fallback_order:
        if fb == "nuc_class" and target_nuc:
            subset = [r for r in remaining if (r.get("features", {}).get("nuc_class") or "").lower() == target_nuc]
        elif fb == "LG" and target_lg:
            subset = [r for r in remaining if (r.get("features", {}).get("LG") or "") == target_lg]
        elif fb == "any":
            subset = remaining[:]
        else:
            subset = []
        cands.extend(subset)
        remaining = [r for r in remaining if r not in subset]
        if len(cands) >= min_candidates:
            break
    return cands


def _similarity(a: Dict[str, Any], b: Dict[str, Any]) -> float:
    # Exact bin match gets perfect similarity
    if (a.get("bin") or "") == (b.get("bin") or "") and a.get("bin"):
        return 1.0

    # Weighted categorical matching
    weights = {
        "LG": 0.35,
        "nuc_class": 0.35,
        "ortho_count": 0.10,
        "para_EWG": 0.10,
        "heteroaryl": 0.10,
    }
    score = 0.0
    total = sum(weights.values())
    for k, w in weights.items():
        av = a.get(k)
        bv = b.get(k)
        # Normalize bools to exact equality
        if isinstance(av, bool) or isinstance(bv, bool):
            if bool(av) == bool(bv):
                score += w
        else:
            if av is not None and bv is not None and str(av).lower() == str(bv).lower():
                score += w

    # Optional small numeric distances if present in feature dicts
    # Use an exponential decay mapped to <= 0.15 extra credit total
    numeric_keys: List[Tuple[str, float, float]] = [
        ("T_C", 50.0, 0.10),  # (scale, weight)
        ("time_h", 8.0, 0.05),
    ]
    for key, scale, w in numeric_keys:
        if key in a and key in b:
            try:
                da = float(a[key]); db = float(b[key])
                import math
                sim_num = math.exp(-abs(da - db) / max(1e-9, scale))
                score += w * sim_num
                total += w
            except Exception:
                # ignore numeric similarity if non-numeric
                pass

    if total <= 0:
        return 0.0
    return max(0.0, min(1.0, score / total))


def knn(family: str, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """
    Retrieve precedents by coarse-bin candidate selection followed by similarity ranking.

    Returns dict with keys: prototype_id, support, precedents[]. If no candidates, returns
    {prototype_id: str, support: 0, precedents: [], error: "NO_PRECEDENTS"}.
    """
    relax = relax or {}
    family_txt = _family_text(family)
    rows = _load()

    # Build candidate set
    cands = _candidate_pool(rows, family_txt, features, k, relax)
    if not cands:
        proto = f"proto_{family_txt.replace(' ', '_').replace('–','-').replace('/','_')}_none_0"
        return {"prototype_id": proto, "support": 0, "precedents": [], "error": "NO_PRECEDENTS"}

    # Score by similarity and yield-weighting
    target_feat = dict(features)
    # Allow bin-derived fallbacks for similarity keys
    if not target_feat.get("LG") or not target_feat.get("nuc_class"):
        bm = _parse_bin(features.get("bin") or "")
        target_feat.setdefault("LG", bm.get("LG"))
        target_feat.setdefault("nuc_class", bm.get("NUC"))

    scored: List[Tuple[float, Dict[str, Any]]] = []
    for r in cands:
        f = r.get("features", {})
        sim = _similarity(target_feat, f)
        if sim <= 0:
            continue
        y = r.get("yield_value")
        y_norm = (float(y) / 100.0) if isinstance(y, (int, float)) else 0.0
        neighbor_score = sim * (0.5 + 0.5 * y_norm)
        scored.append((neighbor_score, r))

    if not scored:
        proto = f"proto_{family_txt.replace(' ', '_').replace('–','-').replace('/','_')}_none_0"
        return {"prototype_id": proto, "support": 0, "precedents": [], "error": "NO_PRECEDENTS"}

    scored.sort(key=lambda x: (-(x[0]), -((x[1].get("yield_value") or 0))))
    top = [r for _, r in scored[: max(1, k)]]
    support = len(scored)

    # Prototype id is a stable-ish hash of family+bin
    family_norm = family_txt.replace(" ", "_").replace("–", "-").replace("/", "_")
    bin_key = str(features.get("bin") or f"LG:{target_feat.get('LG','?')}|NUC:{target_feat.get('nuc_class','?')}")
    proto = f"proto_{family_norm}_{abs(hash(bin_key)) % 100000}"

    precedents = [
        {
            "reaction_id": r.get("reaction_id"),
            "yield": r.get("yield_value"),
            "core": r.get("condition_core"),
            "base_uid": r.get("base_uid"),
            "solvent_uid": r.get("solvent_uid"),
            "T_C": r.get("T_C"),
            "time_h": r.get("time_h"),
        }
        for r in top[:10]
    ]
    return {"prototype_id": proto, "support": support, "precedents": precedents}
