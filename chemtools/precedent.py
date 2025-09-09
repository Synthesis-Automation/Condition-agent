from typing import Dict, Any, List, Tuple, Optional
import os, json
from functools import lru_cache
from .featurizers import molecular as feat_molecular
from . import reaction_similarity as rs

# Local helper to pick electrophile vs nucleophile from reactants list
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

# Path to small demo dataset used for precedents retrieval
DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "reactions_sample.jsonl")
DATASET_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "reaction_dataset")


def _iter_dataset_files() -> List[str]:
    files: List[str] = []
    if os.path.isdir(DATASET_DIR):
        for name in os.listdir(DATASET_DIR):
            if name.lower().endswith(".jsonl"):
                files.append(os.path.join(DATASET_DIR, name))
    return sorted(files)


def _dataset_family_map(raw: str) -> str:
    t = (raw or "").strip()
    # Normalize dataset reaction_type to API family text
    tl = t.lower()
    if tl in {"ullman", "ullmann", "ullman-c-n", "ullmann-c-n", "ullmann c-n"}:
        return "Ullmann C–N"
    if tl in {"buchwald", "buchwald-c-n", "buchwald c-n"}:
        return "Buchwald C–N"
    if tl in {"suzuki", "suzuki-miyaura", "suzuki cc", "suzuki_cc"}:
        return "Suzuki_CC"
    if tl in {"amide-formation", "amide", "amide coupling", "amide_coupling"}:
        return "Amide_Coupling"
    return t


def _make_row_from_dataset(rec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    try:
        rxn_id = rec.get("reaction_id")
        rt = rec.get("reaction_type")
        fam_txt = _dataset_family_map(rt)
        cond = rec.get("conditions") or {}
        y = cond.get("yield_pct")
        T_C = cond.get("temperature_c")
        time_h = cond.get("time_h")
        core = rec.get("condition_core")
        # Base/Solvent: take first entries' CAS where present
        base_uid = None
        for rg in rec.get("reagents", []) or []:
            if (rg.get("role") or "").upper() == "BASE":
                base_uid = rg.get("cas") or rg.get("uid") or rg.get("name")
                if base_uid:
                    break
        solvent_uid = None
        sols = rec.get("solvents", []) or []
        if sols:
            solvent_uid = sols[0].get("cas") or sols[0].get("uid") or sols[0].get("name")
        # Reaction SMILES: reactants>>products (agents side empty in this dataset)
        smiles_block = rec.get("smiles") or {}
        rcts = (smiles_block.get("reactants") or "").strip()
        prods = (smiles_block.get("products") or "").strip()
        rxn_smiles = f"{rcts}>>{prods}"
        # Coarse features for candidate binning (Ullmann currently)
        reactants_list = [p for p in (rcts.split('.') if rcts else []) if p]
        elec, nuc = _pick_electrophile_nucleophile(reactants_list)
        features = feat_molecular.featurize(elec, nuc)
        # Build uniform row
        catalyst_obj = rec.get("catalyst") or {}
        full_system = catalyst_obj.get("full_system") if isinstance(catalyst_obj, dict) else None
        return {
            "reaction_id": rxn_id,
            "rxn_type": fam_txt,
            "yield_value": y,
            "T_C": T_C,
            "time_h": time_h,
            "condition_core": core,
            "base_uid": base_uid,
            "solvent_uid": solvent_uid,
            "reagents": rec.get("reagents") or [],
            "solvents": rec.get("solvents") or [],
            "reference": rec.get("reference") or {},
            "conditions": cond,
            "catalyst": catalyst_obj,
            "full_system": full_system,
            "features": features,
            "reaction_smiles": rxn_smiles,
        }
    except Exception:
        return None


@lru_cache(maxsize=1)
def _load() -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    # 1) Small demo file, if present
    if os.path.exists(DATA_PATH):
        try:
            with open(DATA_PATH, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rows.append(json.loads(line))
                    except Exception:
                        pass
        except Exception:
            pass
    # 2) Auto-load local dataset directory (transformed) when available, unless disabled.
    #    - Explicit override via CHEMTOOLS_LOAD_DATASET (1/true/on to enable; 0/false/off to disable)
    #    - During pytest (PYTEST_CURRENT_TEST set), default to disabled to keep unit tests deterministic
    import os as _os
    _flag = str(_os.environ.get("CHEMTOOLS_LOAD_DATASET", "")).strip().lower()
    if _flag in {"0", "false", "no", "off"}:
        use_dataset = False
    elif _flag in {"1", "true", "yes", "on"}:
        use_dataset = True
    else:
        use_dataset = ("PYTEST_CURRENT_TEST" not in _os.environ)
    if not rows and use_dataset and os.path.isdir(DATASET_DIR):
        for path in _iter_dataset_files():
            try:
                with open(path, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            rec = json.loads(line)
                        except Exception:
                            continue
                        row = _make_row_from_dataset(rec)
                        if row is not None:
                            rows.append(row)
            except Exception:
                continue
    return rows


def _family_text(family: str) -> str:
    # Map API family tokens to dataset labels
    f = (family or "").strip()
    if f.lower() in {"ullmann_cn", "ullmann c–n", "ullmann c-n", "ullmann"}:
        return "Ullmann C–N"
    return f


def _proto_family_id(family_txt: str) -> str:
    """Normalize family text for use in prototype_id.
    Replaces spaces with underscores, normalizes dash characters to '-', and '/' to '_'.
    """
    return (
        str(family_txt)
        .replace(" ", "_")
        .replace("–", "-")  # en dash
        .replace("—", "-")  # em dash
        .replace("−", "-")  # minus sign
        .replace("/", "_")
    )


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
    if len(cands) >= min_candidates:
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


def _as_kv(obj: Dict[str, Any] | None) -> Tuple[Tuple[str, Any], ...]:
    if not obj:
        return tuple()
    # Convert dict to a stable, hashable key
    return tuple(sorted((str(k), obj[k]) for k in obj))


@lru_cache(maxsize=512)
def _knn_cached(family: str, features_kv: Tuple[Tuple[str, Any], ...], k: int, relax_kv: Tuple[Tuple[str, Any], ...]) -> Dict[str, Any]:
    features = {k: v for k, v in features_kv}
    relax = {k: v for k, v in relax_kv}
    return _knn_impl(family, features, k, relax)


def knn(family: str, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
    # Public entrypoint; wrap cached implementation. Return a copy to avoid shared-state mutations.
    out = _knn_cached(family, _as_kv(features or {}), int(k), _as_kv(relax or {}))
    return {**out}


def _knn_impl(family: str, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
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
        proto = f"proto_{_proto_family_id(family_txt)}_none_0"
        return {"prototype_id": proto, "support": 0, "precedents": [], "error": "NO_PRECEDENTS"}

    # Score by similarity and yield-weighting
    target_feat = dict(features)
    # Allow bin-derived fallbacks for similarity keys
    if not target_feat.get("LG") or not target_feat.get("nuc_class"):
        bm = _parse_bin(features.get("bin") or "")
        target_feat.setdefault("LG", bm.get("LG"))
        target_feat.setdefault("nuc_class", bm.get("NUC"))

    # Optional DRFP re-ranking (best-effort)
    use_drfp = bool(relax.get("use_drfp", False))
    rsmi_query = str(relax.get("reaction_smiles") or "")
    drfp_w = float(relax.get("drfp_weight", 0.4))
    drfp_bits = int(relax.get("drfp_n_bits", 4096))
    drfp_radius = int(relax.get("drfp_radius", 3))
    q_fp = None
    if use_drfp and rsmi_query and rs.drfp_available():
        # Optionally precompute fingerprints for entire dataset to warm cache
        if bool(relax.get("precompute_drfp", False)):
            try:
                # Touch encode cache for all rows to speed subsequent scoring
                for _r in cands if relax.get("precompute_scope") == "candidates" else rows:
                    rsmi_val = _r.get("reaction_smiles")
                    if rsmi_val:
                        _ = rs.encode_drfp_cached(rsmi_val, n_bits=drfp_bits, radius=drfp_radius)
            except Exception:
                pass
        q_fp = rs.encode_drfp_cached(rsmi_query, n_bits=drfp_bits, radius=drfp_radius)

    scored: List[Tuple[float, Dict[str, Any]]] = []
    for r in cands:
        f = r.get("features", {})
        sim_cat = _similarity(target_feat, f)
        if sim_cat <= 0:
            # still allow DRFP to rescue a bit if enabled
            pass
        sim_total = sim_cat
        # DRFP component when available for both
        if q_fp is not None:
            r_rsmi = r.get("reaction_smiles")
            if r_rsmi:
                r_fp = rs.encode_drfp_cached(r_rsmi, n_bits=drfp_bits, radius=drfp_radius)
                if r_fp is not None:
                    sim_fp = rs.tanimoto(q_fp, r_fp)
                    try:
                        sim_fp = float(sim_fp)  # coerce in case implementation returns non-float
                    except Exception:
                        sim_fp = 0.0
                    # Blend: convex combination of categorical and DRFP
                    sim_total = max(0.0, min(1.0, (1.0 - drfp_w) * sim_cat + drfp_w * sim_fp))
        if sim_total <= 0:
            continue
        y = r.get("yield_value")
        y_norm = (float(y) / 100.0) if isinstance(y, (int, float)) else 0.0
        neighbor_score = sim_total * (0.5 + 0.5 * y_norm)
        scored.append((neighbor_score, r))

    if not scored:
        proto = f"proto_{_proto_family_id(family_txt)}_none_0"
        return {"prototype_id": proto, "support": 0, "precedents": [], "error": "NO_PRECEDENTS"}

    scored.sort(key=lambda x: (-(x[0]), -((x[1].get("yield_value") or 0))))
    top = [r for _, r in scored[: max(1, k)]]
    support = len(scored)

    # Prototype id is a stable-ish hash of family+bin
    family_norm = _proto_family_id(family_txt)
    bin_key = str(features.get("bin") or f"LG:{target_feat.get('LG','?')}|NUC:{target_feat.get('nuc_class','?')}")
    proto = f"proto_{family_norm}_{abs(hash(bin_key)) % 100000}"

    # Include reaction SMILES and parsed sides for UI/consumers
    try:
        from .smiles import _split_reaction_smiles as _split_rxn  # type: ignore
    except Exception:
        def _split_rxn(rsmi: str):
            parts = (rsmi or "").split(">")
            if len(parts) == 2 and ">>" in (rsmi or ""):
                return parts[0], "", parts[1]
            if len(parts) == 3:
                return parts[0], parts[1], parts[2]
            return rsmi, "", ""

    precedents = []
    for r in top[:10]:
        rsmi = r.get("reaction_smiles") or ""
        reactants_smi, _agents_smi, products_smi = _split_rxn(rsmi)
        precedents.append({
            "reaction_id": r.get("reaction_id"),
            "reaction_smiles": rsmi,
            "reactants_smiles": reactants_smi,
            "products_smiles": products_smi,
            "condition_core": r.get("condition_core"),
            "yield": r.get("yield_value"),
            "core": r.get("condition_core"),
            "base_uid": r.get("base_uid"),
            "solvent_uid": r.get("solvent_uid"),
            "reagents": r.get("reagents"),
            "solvents": r.get("solvents"),
            "reference": r.get("reference"),
            "conditions": r.get("conditions"),
            "catalyst": r.get("catalyst"),
            "full_system": r.get("full_system"),
            "T_C": r.get("T_C"),
            "time_h": r.get("time_h"),
        })
    return {"prototype_id": proto, "support": support, "precedents": precedents}
