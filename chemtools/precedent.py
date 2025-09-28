from typing import Dict, Any, List, Tuple, Optional, Mapping, Sequence
import os, json
from functools import lru_cache
from .featurizers import molecular as feat_molecular
from . import reaction_similarity as rs
from collections import Counter

try:
    from .integrations import molpipeline as molpipe_integration
    _HAS_MOLPIPELINE = molpipe_integration.is_available()
except Exception:
    molpipe_integration = None  # type: ignore
    _HAS_MOLPIPELINE = False

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
BASE_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_DIR, "data", "reactions_sample.jsonl")
# Allow overriding the dataset directory so data-processor can save directly there
_ENV_DIR = os.environ.get("CHEMTOOLS_DATASET_DIR", "").strip()
DATASET_DIR = (
    os.path.abspath(_ENV_DIR) if _ENV_DIR else os.path.join(BASE_DIR, "data", "reaction_dataset")
)


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
        return "Ullmann C鈥揘"
    if tl in {"buchwald", "buchwald-c-n", "buchwald c-n"}:
        return "Buchwald C鈥揘"
    if tl in {"suzuki", "suzuki-miyaura", "suzuki cc", "suzuki_cc"}:
        return "Suzuki_CC"
    if tl in {"amide-formation", "amide formation", "amideformation", "amide", "amide coupling", "amide_coupling"}:
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
    if f.lower() in {"ullmann_cn", "ullmann c鈥搉", "ullmann c-n", "ullmann"}:
        return "Ullmann C鈥揘"
    return f


def _proto_family_id(family_txt: str) -> str:
    """Normalize family text for use in prototype_id."""
    txt = str(family_txt).replace(' ', '_')
    for old in ('–', '—', '−'):
        txt = txt.replace(old, '-')
    return txt.replace('/', '_')


def _parse_bin(bin_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not bin_str:
        return out
    for part in str(bin_str).split('|'):
        if ':' in part:
            k, v = part.split(':', 1)
            out[k.strip()] = v.strip()
    return out
    out: Dict[str, str] = {}
    if not bin_str:
        return out
    for part in str(bin_str).split("|"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out

###############################################
# Catalyst class helpers and filtering
###############################################

_METAL_NAME_TO_SYMBOL = {
    # common transition metals and aliases
    "palladium": "Pd", "pd": "Pd",
    "nickel": "Ni", "ni": "Ni",
    "cobalt": "Co", "co": "Co",
    "copper": "Cu", "cu": "Cu",
    "iron": "Fe", "fe": "Fe",
    "ruthenium": "Ru", "ru": "Ru",
    "rhodium": "Rh", "rh": "Rh",
    "iridium": "Ir", "ir": "Ir",
    "gold": "Au", "au": "Au",
    "silver": "Ag", "ag": "Ag",
    "zinc": "Zn", "zn": "Zn",
    "magnesium": "Mg", "mg": "Mg",
    "manganese": "Mn", "mn": "Mn",
    "chromium": "Cr", "cr": "Cr",
    "vanadium": "V", "v": "V",
    "titanium": "Ti", "ti": "Ti",
    "zirconium": "Zr", "zr": "Zr",
    "molybdenum": "Mo", "mo": "Mo",
    "tungsten": "W", "w": "W",
    "rhenium": "Re", "re": "Re",
    "osmium": "Os", "os": "Os",
    "platinum": "Pt", "pt": "Pt",
}

_METAL_SYMBOLS = set(_METAL_NAME_TO_SYMBOL.values())
_ENZYME_KEYWORDS = {"enzyme", "protein", "lipase", "oxidase", "dehydrogenase", "transferase"}

def _normalize_symbol(token: str) -> Optional[str]:
    t = (token or "").strip()
    if not t:
        return None
    lo = t.lower()
    if lo in _METAL_NAME_TO_SYMBOL:
        return _METAL_NAME_TO_SYMBOL[lo]
    # Try exact symbol match case-insensitively
    up = t[0].upper() + t[1:].lower()
    if up in _METAL_SYMBOLS:
        return up
    return None

def _row_catalyst_class(row: Dict[str, Any]) -> str:
    """Heuristically classify a precedent row into a catalyst class.

    Returns one of metal symbols (e.g., 'Pd', 'Ni', ...), 'enzyme', 'organo_catalyst', or 'other'.
    """
    # 1) Use condition_core like 'Pd/XPhos'
    core = (row.get("condition_core") or "").strip()
    if core:
        head = core.split("/", 1)[0].strip()
        sym = _normalize_symbol(head)
        if sym:
            return sym
    # 2) Scan full_system names for explicit metal names
    fs = row.get("full_system")
    if isinstance(fs, list):
        names = [str((it or {}).get("name") or "") for it in fs]
        text = " ".join(names).lower()
        # enzyme detection first
        if any(k in text for k in _ENZYME_KEYWORDS):
            return "enzyme"
        # metal names
        for key, sym in _METAL_NAME_TO_SYMBOL.items():
            if len(key) <= 2:  # skip symbol aliases here; handled below
                continue
            if key in text:
                return sym
        # symbol occurrence as word-ish
        for sym in _METAL_SYMBOLS:
            if f" {sym.lower()}" in text or f"({sym.lower()}" in text or f"[{sym.lower()}" in text:
                return sym
    # 3) Fallback enzyme detection in catalyst dict
    cat = row.get("catalyst") or {}
    if isinstance(cat, dict):
        nm = str(cat.get("name") or "").lower()
        if any(k in nm for k in _ENZYME_KEYWORDS):
            return "enzyme"
    # 4) If no metal detected, assume organocatalyst when there is any catalyst/ligand info
    if core or (isinstance(fs, list) and fs):
        return "organo_catalyst"
    return "other"

def _match_catalyst_class(selected: str, row_cls: str) -> bool:
    sel = (selected or "").strip().lower()
    if not sel:
        return True
    if sel in {"organo_catalyst", "enzyme", "other"}:
        return row_cls == sel
    sym = _normalize_symbol(sel)
    if sym:
        return row_cls == sym
    return True  # unknown filter -> do not exclude


def _candidate_pool(rows: List[Dict[str, Any]], family_txt: str, feat: Dict[str, Any], k: int, relax: Dict[str, Any]) -> List[Dict[str, Any]]:
    # Filter rows to family first
    fam_rows = [r for r in rows if (r.get("rxn_type") or "") == family_txt]
    # Optional catalyst class filter
    cat_filter = None
    try:
        cat_filter = str(relax.get("catalyst_class")) if isinstance(relax, dict) and relax.get("catalyst_class") is not None else None
    except Exception:
        cat_filter = None
    if cat_filter:
        fam_rows = [r for r in fam_rows if _match_catalyst_class(str(cat_filter), _row_catalyst_class(r))]
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
    relax_dict = dict(relax or {})
    molpipeline_cfg = relax_dict.pop('molpipeline', None)
    out = _knn_cached(family, _as_kv(features or {}), int(k), _as_kv(relax_dict))
    pack = {**out}
    if molpipeline_cfg:
        pack = _attach_molpipeline_features(pack, molpipeline_cfg)
    return pack


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
        # DRFP component when available for both; prefer whole-reaction similarity
        if q_fp is not None:
            r_rsmi = r.get("reaction_smiles")
            if r_rsmi:
                r_fp = rs.encode_drfp_cached(r_rsmi, n_bits=drfp_bits, radius=drfp_radius)
                if r_fp is not None:
                    sim_fp = rs.tanimoto(q_fp, r_fp)
                    try:
                        sim_fp = float(sim_fp)
                    except Exception:
                        sim_fp = 0.0
                    sim_total = max(0.0, min(1.0, sim_fp))
        if sim_total <= 0:
            if sim_cat > 0:
                sim_total = sim_cat
            else:
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


def _parse_core_tokens(core_text: str) -> Tuple[str, str]:
    """Parse a condition_core string like 'Pd/XPhos' into (metal, ligand) tokens (lowercased).

    Returns (metal, ligand) where any missing part is an empty string.
    """
    t = (core_text or "").strip()
    if not t:
        return "", ""
    if "/" in t:
        a, b = t.split("/", 1)
        return (a or "").strip().lower(), (b or "").strip().lower()
    return (t.strip().lower(), "")


def _norm_family(fam: Optional[str]) -> Optional[str]:
    if fam is None:
        return None
    f = (fam or "").strip()
    if not f:
        return None
    return _family_text(f)


def find_reactions_by_core(
    core_query: str,
    *,
    family: Optional[str] = None,
    fuzzy: bool = True,
    limit: int = 50,
) -> List[Dict[str, Any]]:
    """Find dataset reactions that use the same or similar condition core.

    - core_query: e.g., 'Pd/XPhos', 'Pd', or 'XPhos'.
    - family: optional reaction family text (e.g., 'Ullmann C鈥揘', 'Suzuki_CC').
    - fuzzy: if True, also match ligand names from catalyst/full_system entries.
    - limit: maximum number of results.

    Returns a list of dataset rows in the normalized precedent format
    (see _make_row_from_dataset), possibly truncated to `limit`.
    """
    q_metal, q_lig = _parse_core_tokens(core_query)
    fam_norm = _norm_family(family)

    rows = _load()
    out: List[Dict[str, Any]] = []

    def match_row(r: Dict[str, Any]) -> bool:
        if fam_norm and (r.get("rxn_type") or "") != fam_norm:
            return False
        rc = r.get("condition_core") or ""
        rm, rl = _parse_core_tokens(rc)
        # Basic matching rules
        m_ok = (not q_metal) or (rm == q_metal or q_metal in rm)
        l_ok = (not q_lig) or (rl == q_lig or (q_lig and q_lig in rl))
        if m_ok and l_ok:
            return True
        if not fuzzy:
            return False
        # Fuzzy ligand check against catalyst/full_system names
        if q_lig:
            fs = r.get("full_system")
            if isinstance(fs, list):
                for it in fs:
                    nm = str((it or {}).get("name") or "").strip().lower()
                    if nm and q_lig in nm:
                        # Respect metal if provided
                        if not q_metal or m_ok:
                            return True
        return False

    for row in rows:
        try:
            if match_row(row):
                out.append(row)
                if len(out) >= int(limit):
                    break
        except Exception:
            continue
    return out


def list_cores(
    *,
    family: Optional[str] = None,
    top_n: Optional[int] = None,
    include_counts: bool = True,
) -> List[Dict[str, Any]] | List[str]:
    """List unique condition cores from the loaded reaction dataset.

    - family: optional reaction family text to filter (e.g., 'Ullmann C鈥揘', 'Suzuki_CC').
    - top_n: optional cap on number of items returned (by frequency desc).
    - include_counts: when True, return list of {core, count}; else list of core strings.
    """
    fam_norm = _norm_family(family)
    rows = _load()
    ctr: Counter[str] = Counter()
    for r in rows:
        if fam_norm and (r.get("rxn_type") or "") != fam_norm:
            continue
        c = (r.get("condition_core") or "").strip()
        if c:
            ctr[c] += 1
    items = ctr.most_common()
    if top_n is not None:
        items = items[: int(top_n)]
    if include_counts:
        return [{"core": c, "count": n} for c, n in items]
    return [c for c, _ in items]


def _attach_molpipeline_features(pack: Dict[str, Any], cfg: Any) -> Dict[str, Any]:
    if not isinstance(pack, dict):
        return pack
    if not isinstance(cfg, Mapping):
        return pack
    if not _HAS_MOLPIPELINE or molpipe_integration is None:
        return pack
    suppress = bool(cfg.get('suppress_errors', True))
    try:
        aggregator = cfg.get('aggregator') if isinstance(cfg, Mapping) else None
        if aggregator is None:
            roles = cfg.get('roles') if isinstance(cfg.get('roles'), Sequence) else None
            aggregator = molpipe_integration.build_default_role_aggregator(
                roles=roles,
                aggregate=str(cfg.get('aggregate', 'mean')),
                missing_strategy=str(cfg.get('missing_strategy', 'zeros')),
                n_jobs=int(cfg.get('n_jobs', 1)),
                ligand_n_bits=int(cfg.get('ligand_n_bits', 512)),
                ligand_radius=int(cfg.get('ligand_radius', 2)),
            )
    except Exception as exc:
        if suppress:
            pack.setdefault('molpipeline_warnings', []).append(str(exc))
            return pack
        raise
    include_role = bool(cfg.get('include_role_features', True))
    include_concat = bool(cfg.get('include_concat', True))
    role_key = str(cfg.get('role_features_key', 'molpipeline_role_features'))
    concat_key = str(cfg.get('concat_key', 'molpipeline_feature_vector'))
    precedents = pack.get('precedents')
    if isinstance(precedents, list):
        for row in precedents:
            if not isinstance(row, dict):
                continue
            try:
                role_feats = aggregator.featurize_roles(reaction=row)
                if include_role:
                    row[role_key] = {r: (vec.tolist() if hasattr(vec, 'tolist') else vec) if vec is not None else None for r, vec in role_feats.items()}
                if include_concat:
                    row[concat_key] = aggregator.concatenate(reaction=row).tolist()
            except Exception as exc:
                if suppress:
                    row.setdefault('molpipeline_error', str(exc))
                else:
                    raise
    try:
        query_role_smiles = cfg.get('query_role_smiles') if isinstance(cfg, Mapping) else None
        if query_role_smiles is None and isinstance(cfg.get('query_reaction'), Mapping):
            query_role_smiles = molpipe_integration.collect_role_smiles(cfg['query_reaction'])
        if query_role_smiles:
            pack['molpipeline_query_vector'] = aggregator.concatenate(role_smiles=query_role_smiles).tolist()
    except Exception as exc:
        if suppress:
            pack.setdefault('molpipeline_warnings', []).append(str(exc))
        else:
            raise
    return pack
