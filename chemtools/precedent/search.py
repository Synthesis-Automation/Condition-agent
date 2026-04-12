#utf-8
"""k-NN precedent search and core-based reaction lookup.

Main search functionality for finding similar reactions based on molecular
features, catalyst class, and reaction conditions.
"""
from typing import Dict, Any, List, Tuple, Optional
from functools import lru_cache
from collections import Counter
import os
import time
import logging

from .loader import _load, _load_selective
from .core_utils import _family_text, _proto_family_id, _parse_bin, _parse_core_tokens, _norm_family
from .catalyst import _row_catalyst_class, _match_catalyst_class
from .similarity import _similarity

# Configure logging for performance tracking
logger = logging.getLogger(__name__)

# Import reaction_similarity for DRFP support
try:
    from .. import reaction_similarity as rs
except ImportError:
    rs = None  # type: ignore

# Import DRFP binary storage utilities
try:
    from ..util.drfp_storage import DRFPLoader, get_drfp_path_for_family, get_unified_drfp_path
    _drfp_storage_available = True
except ImportError:
    DRFPLoader = None  # type: ignore
    get_drfp_path_for_family = None  # type: ignore
    get_unified_drfp_path = None  # type: ignore
    _drfp_storage_available = False

# Global cache for DRFP loaders (one per family)
_DRFP_LOADER_CACHE: Dict[str, Any] = {}

# Import MixFP + KMN + FAISS routing utilities (MOSAIC-style)
try:
    from ..util.mix_fingerprint import create_mix_fp as _create_mix_fp
    from ..util.kernel_metric_net import get_default_kmn as _get_default_kmn
    from ..util.faiss_router import get_default_router as _get_default_router, is_index_built as _is_index_built
    _MIXFP_AVAILABLE = True
except Exception:
    _create_mix_fp = None   # type: ignore
    _get_default_kmn = None  # type: ignore
    _get_default_router = None  # type: ignore
    _is_index_built = None  # type: ignore
    _MIXFP_AVAILABLE = False

# Singleton KMN / FAISS router (loaded lazily on first use)
_KMN_INSTANCE: Any = None
_FAISS_ROUTER_INSTANCE: Any = None
_DEFAULT_ON_DEMAND_DRFP_RERANK_LIMIT = 400


def _candidate_pool(rows: List[Dict[str, Any]], family_txt: str, feat: Dict[str, Any], k: int, relax: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build candidate pool of precedents matching the query family and features.
    
    Uses relaxation strategies to ensure sufficient candidates are found.
    Optimized with set operations for O(n) performance instead of O(n²).
    
    Args:
        rows: All loaded precedent rows
        family_txt: Canonical family name (or None/"All" for cross-family search)
        feat: Query feature dict
        k: Target number of results
        relax: Relaxation configuration dict
        
    Returns:
        List of candidate precedent rows
    """
    # Filter rows to family first (or use all if family is None/"All")
    if family_txt is None or family_txt.lower() in ["all", "none", ""]:
        fam_rows = rows  # Use all reactions regardless of family
    else:
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
    # Global cap on the candidate pool to prevent brute-force scoring of huge sets.
    # Applied after all fallbacks; 0 means unlimited.  Default is 2000 rows.
    _max_cands = int(relax.get("max_any_candidates", 2000))
    fallback_order: List[str] = relax.get("fallback_order", ["any"])  # type: ignore

    target_bin = (feat.get("bin") or "").strip()
    target_bin_map = _parse_bin(target_bin)
    target_nuc = (feat.get("nuc_class") or target_bin_map.get("NUC") or "").lower()
    target_lg = feat.get("LG") or target_bin_map.get("LG") or ""

    # Exact bin matches
    cands = [r for r in fam_rows if (r.get("features", {}).get("bin") or "") == target_bin]
    if len(cands) >= min_candidates:
        if _max_cands > 0 and len(cands) > _max_cands:
            cands.sort(key=lambda r: float(r.get("yield_value") or 0), reverse=True)
            return cands[:_max_cands]
        return cands

    # Use sets for O(1) lookups instead of O(n) list operations
    # Track reaction IDs to avoid duplicates efficiently
    added_ids = {id(r) for r in cands}

    # Fallbacks with optimized set-based tracking
    for fb in fallback_order:
        if len(cands) >= min_candidates:
            break

        subset = []
        if fb == "nuc_class" and target_nuc:
            subset = [r for r in fam_rows
                     if id(r) not in added_ids
                     and (r.get("features", {}).get("nuc_class") or "").lower() == target_nuc]
        elif fb == "LG" and target_lg:
            subset = [r for r in fam_rows
                     if id(r) not in added_ids
                     and (r.get("features", {}).get("LG") or "") == target_lg]
        elif fb == "any":
            subset = [r for r in fam_rows if id(r) not in added_ids]

        # Add new candidates and update tracking set
        cands.extend(subset)
        added_ids.update(id(r) for r in subset)

    # Apply global cap: sort by yield desc so the highest-quality precedents are kept
    if _max_cands > 0 and len(cands) > _max_cands:
        cands.sort(key=lambda r: float(r.get("yield_value") or 0), reverse=True)
        return cands[:_max_cands]
    return cands


def _resolve_drfp_loader(family_txt: Optional[str]) -> Any:
    if not _drfp_storage_available or DRFPLoader is None:
        return None
    if family_txt is None:
        if "__UNIFIED__" not in _DRFP_LOADER_CACHE and get_unified_drfp_path is not None:
            unified_path = get_unified_drfp_path()
            if os.path.exists(unified_path):
                try:
                    _DRFP_LOADER_CACHE["__UNIFIED__"] = DRFPLoader(unified_path)
                except Exception:
                    _DRFP_LOADER_CACHE["__UNIFIED__"] = None
            else:
                _DRFP_LOADER_CACHE["__UNIFIED__"] = None
        return _DRFP_LOADER_CACHE.get("__UNIFIED__")

    if get_drfp_path_for_family is None:
        return None
    if family_txt not in _DRFP_LOADER_CACHE:
        npz_path = get_drfp_path_for_family(family_txt)
        if os.path.exists(npz_path):
            try:
                _DRFP_LOADER_CACHE[family_txt] = DRFPLoader(npz_path)
            except Exception:
                _DRFP_LOADER_CACHE[family_txt] = None
        else:
            _DRFP_LOADER_CACHE[family_txt] = None
    return _DRFP_LOADER_CACHE.get(family_txt)


def _as_kv(obj: Dict[str, Any] | None) -> Tuple[Tuple[str, Any], ...]:
    """Convert dict to hashable tuple for caching."""
    if not obj:
        return tuple()
    # Convert dict to a stable, hashable key
    return tuple(sorted((str(k), obj[k]) for k in obj))


@lru_cache(maxsize=512)
def _knn_cached(family: str, features_kv: Tuple[Tuple[str, Any], ...], k: int, relax_kv: Tuple[Tuple[str, Any], ...]) -> Dict[str, Any]:
    """Cached k-NN implementation. Converts tuples back to dicts and calls _knn_impl."""
    features = {k: v for k, v in features_kv}
    relax = {k: v for k, v in relax_kv}
    return _knn_impl(family, features, k, relax)


def knn(family: str | None, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """Find k nearest neighbor precedents for a reaction.
    
    Main public API for precedent search. Uses cached implementation for performance.
    
    Args:
        family: Reaction family name (e.g., "C_N_Coupling", "Suzuki"), or None/"All" 
                to search across all reaction families (cross-family search)
        features: Query feature dict (bin, LG, nuc_class, etc.)
        k: Number of precedents to return (default: 50)
        relax: Optional relaxation/configuration dict with options:
            - filter_by_reagent_database (bool): Only return precedents where all 
              reagents are found in the reagent database (default: True)
            - use_drfp (bool): Use DRFP-based similarity (default: False)
            - selective_loading (bool): Load only requested family (default: True)
              Note: Automatically set to False when family=None/"All"
            - debug_timing (bool): Log timing information (default: False)
            - precedent_limit (int): Max precedents to return (default: 10)
        
    Returns:
        Dict with keys: prototype_id, support, precedents, (optional) error
        
    Example:
        >>> # Get precedents with database filtering (default)
        >>> pack = knn("C_N_Coupling", features={}, k=25)
        >>> 
        >>> # Search across all reaction families (cross-family search)
        >>> pack = knn(None, features={}, k=25, 
        ...            relax={"use_drfp": True, "reaction_smiles": "..."})
        >>> 
        >>> # Disable database filtering to get all precedents
        >>> pack = knn("C_N_Coupling", features={}, k=25, 
        ...            relax={"filter_by_reagent_database": False})
    """
    relax_dict = dict(relax or {})
    
    # Set default: enable reagent database filtering
    relax_dict.setdefault("filter_by_reagent_database", True)
    
    # For cross-family search, disable selective loading (must load all families)
    if family is None or (isinstance(family, str) and family.lower() in ["all", "none", ""]):
        relax_dict["selective_loading"] = False
    
    family_key = family if family is not None else "__ALL__"  # Use special key for caching
    out = _knn_cached(family_key, _as_kv(features or {}), int(k), _as_kv(relax_dict))
    return {**out}


def _knn_impl(family: str | None, features: Dict[str, Any], k: int = 50, relax: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """
    Retrieve precedents by coarse-bin candidate selection followed by similarity ranking.

    Returns dict with keys: prototype_id, support, precedents[]. If no candidates, returns
    {prototype_id: str, support: 0, precedents: [], error: "NO_PRECEDENTS"}.
    
    Args:
        family: Reaction family name, or None/"All"/"__ALL__" for cross-family search
    """
    # ⏱️ START TIMING
    t_start_total = time.time()
    timing = {}
    
    relax = relax or {}
    
    # Handle cross-family search (family=None or special markers)
    if family is None or family in ["__ALL__", "All", "all", "None", "none", ""]:
        family_txt = None  # Signal for cross-family search
    else:
        family_txt = _family_text(family)
    
    # Selective loading: only load the requested family to improve performance
    # Check if selective loading is enabled (default: True for better performance)
    use_selective_loading = relax.get("selective_loading", True)
    
    # ⏱️ TIME: Data loading
    t_start_load = time.time()
    if use_selective_loading and family_txt is not None:
        # Load only the requested family
        rows = _load_selective(families=[family_txt])
    else:
        # Load all datasets (legacy behavior OR cross-family search)
        rows = _load()
    timing['load_data'] = time.time() - t_start_load
    
    # ⏱️ TIME: Candidate pool building
    t_start_candidates = time.time()

    # ── MixFP + KMN + FAISS routing (MOSAIC-style) ────────────────────────
    # When `use_mixfp=True`, route the query fingerprint through the KMN
    # and FAISS index to get an ordered list of similar reaction IDs. These
    # are then used as an additional similarity signal (or, when
    # `mixfp_routing_only=True`, as the sole candidate pool).
    #
    # Auto-detection: if the caller did not explicitly set use_mixfp, enable
    # it automatically whenever (a) the KMN index files exist on disk, and
    # (b) a reaction_smiles was provided.  Pass use_mixfp=False to opt out.
    _use_mixfp_explicit = relax.get("use_mixfp")
    if _use_mixfp_explicit is None:
        # Auto-detect: enable if index is built + smiles is available
        _has_smiles = bool(relax.get("reaction_smiles"))
        use_mixfp = (
            _MIXFP_AVAILABLE
            and _has_smiles
            and _is_index_built is not None
            and _is_index_built()
        )
    else:
        use_mixfp = bool(_use_mixfp_explicit)
    mixfp_scores: Dict[str, float] = {}   # reaction_id → cosine similarity
    mixfp_routing_only = bool(relax.get("mixfp_routing_only", False))
    mixfp_w = float(relax.get("mixfp_weight", 0.6))
    if mixfp_w < 0.0:
        mixfp_w = 0.0
    elif mixfp_w > 1.0:
        mixfp_w = 1.0

    if use_mixfp and _MIXFP_AVAILABLE and _create_mix_fp is not None:
        global _KMN_INSTANCE, _FAISS_ROUTER_INSTANCE
        try:
            # Lazy-load singletons
            if _KMN_INSTANCE is None:
                _KMN_INSTANCE = _get_default_kmn()  # type: ignore[misc]
            if _FAISS_ROUTER_INSTANCE is None:
                _FAISS_ROUTER_INSTANCE = _get_default_router()  # type: ignore[misc]

            _rsmi_query = str(relax.get("reaction_smiles") or "")
            _fp_size = int(relax.get("mixfp_fp_size", 1024))
            _k_route = max(int(k * 3), 150)  # fetch more candidates for re-ranking

            if _rsmi_query and _FAISS_ROUTER_INSTANCE.is_ready:
                _fp = _create_mix_fp(_rsmi_query, fp_size=_fp_size)
                if _fp is not None:
                    _emb = _KMN_INSTANCE.get_embeddings(_fp.reshape(1, -1))[0]
                    _ids, _sims = _FAISS_ROUTER_INSTANCE.search(_emb, k=_k_route)
                    mixfp_scores = {
                        rid: float(sim)
                        for rid, sim in zip(_ids, _sims)
                        if rid
                    }
                    if relax.get("debug_timing", False):
                        logger.info(
                            "   MixFP routing: %d candidates (top sim=%.3f)",
                            len(mixfp_scores),
                            max(mixfp_scores.values(), default=0.0),
                        )
        except Exception as _mixfp_exc:
            logger.debug("MixFP routing error (non-fatal): %s", _mixfp_exc)
            mixfp_scores = {}

    # ── Build candidate pool ───────────────────────────────────────────────
    if mixfp_routing_only and mixfp_scores:
        # Use FAISS routing results directly: keep only rows whose reaction_id
        # appears in the FAISS top-k (preserving cross-family results when
        # family=None and mixfp_routing_only=True).
        _rows_by_id = {r.get("reaction_id"): r for r in rows if r.get("reaction_id")}
        cands = [
            _rows_by_id[rid]
            for rid in mixfp_scores
            if rid in _rows_by_id
        ]
        # Apply any family filter if not routing globally
        if family_txt is not None:
            cands = [r for r in cands if (r.get("rxn_type") or "") == family_txt]
    else:
        # Standard categorical candidate pool (existing logic)
        cands = _candidate_pool(rows, family_txt, features, k, relax)
    timing['build_candidates'] = time.time() - t_start_candidates

    if not cands:
        proto_family = _proto_family_id(family_txt) if family_txt else "all_families"
        proto = f"proto_{proto_family}_none_0"
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
    drfp_w = float(relax.get("drfp_weight", 0.7))
    if drfp_w < 0.0:
        drfp_w = 0.0
    elif drfp_w > 1.0:
        drfp_w = 1.0
    drfp_bits = int(relax.get("drfp_n_bits", 4096))
    drfp_radius = int(relax.get("drfp_radius", 3))
    q_fp = None
    drfp_loader = None
    on_demand_drfp_limit = int(
        relax.get(
            "drfp_rerank_limit",
            max(_DEFAULT_ON_DEMAND_DRFP_RERANK_LIMIT, max(int(k) * 8, 0)),
        )
    )
    if on_demand_drfp_limit <= 0:
        on_demand_drfp_limit = len(cands)
    
    # ⏱️ TIME: Query DRFP encoding
    t_start_drfp_query = time.time()
    if use_drfp and rsmi_query and rs and rs.drfp_available():
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
        drfp_loader = _resolve_drfp_loader(family_txt)
    timing['drfp_query_encode'] = time.time() - t_start_drfp_query

    # ⏱️ TIME: Similarity scoring
    t_start_scoring = time.time()
    drfp_load_count = {'binary': 0, 'jsonl': 0, 'computed': 0}
    candidate_rows: List[Tuple[Dict[str, Any], Dict[str, Any], float, float]] = []
    for r in cands:
        f = r.get("features", {})
        sim_cat = _similarity(target_feat, f)
        y = r.get("yield_value")
        y_norm = (float(y) / 100.0) if isinstance(y, (int, float)) else 0.0
        candidate_rows.append((r, f, sim_cat, y_norm))

    allowed_on_demand_ids = None
    if q_fp is not None and drfp_loader is None and len(candidate_rows) > on_demand_drfp_limit:
        ranked_for_drfp = sorted(
            candidate_rows,
            key=lambda item: (item[2], item[3]),
            reverse=True,
        )
        allowed_on_demand_ids = {
            id(item[0])
            for item in ranked_for_drfp[:on_demand_drfp_limit]
        }

    scored: List[Tuple[float, Dict[str, Any]]] = []
    for r, f, sim_cat, y_norm in candidate_rows:
        if sim_cat <= 0:
            # still allow DRFP to rescue a bit if enabled
            pass
        sim_total = sim_cat
        r_fp = None    # initialise so MixFP blend can always reference safely
        sim_fp = 0.0   # initialise so MixFP blend can always reference safely
        # DRFP component when available for both; prefer whole-reaction similarity
        if q_fp is not None:
            r_fp = None
            
            # STRATEGY 1: Try to load from binary NPZ file (fastest, most space-efficient)
            if drfp_loader is not None:
                reaction_id = r.get("reaction_id")
                if reaction_id:
                    try:
                        r_fp = drfp_loader.get_fingerprint(reaction_id)
                        if r_fp is not None:
                            drfp_load_count['binary'] += 1
                    except Exception:
                        pass
            
            # STRATEGY 2: Try to load precomputed DRFP from JSONL (legacy fallback)
            if r_fp is None:
                precomp = r.get("precomputed", {})
                if isinstance(precomp, dict):
                    drfp_fp_list = precomp.get("drfp_fp")
                    precomp_bits = precomp.get("drfp_n_bits", 4096)
                    precomp_radius = precomp.get("drfp_radius", 3)
                    
                    # Use precomputed FP if parameters match
                    if (drfp_fp_list is not None 
                        and precomp_bits == drfp_bits 
                        and precomp_radius == drfp_radius):
                        try:
                            import numpy as np
                            r_fp = np.array(drfp_fp_list, dtype='uint8')
                            if r_fp is not None:
                                drfp_load_count['jsonl'] += 1
                        except Exception:
                            pass  # Fall back to computing
            
            # STRATEGY 3: Fall back to computing on-demand if not precomputed
            if r_fp is None:
                allow_compute = allowed_on_demand_ids is None or id(r) in allowed_on_demand_ids
                if allow_compute:
                    r_rsmi = r.get("reaction_smiles")
                    if r_rsmi:
                        r_fp = rs.encode_drfp_cached(r_rsmi, n_bits=drfp_bits, radius=drfp_radius)
                        if r_fp is not None:
                            drfp_load_count['computed'] += 1
            
            if r_fp is not None:
                sim_fp = rs.tanimoto(q_fp, r_fp)
                try:
                    sim_fp = float(sim_fp)
                except Exception:
                    sim_fp = 0.0
                sim_fp = max(0.0, min(1.0, sim_fp))
                if sim_cat > 0:
                    sim_total = (drfp_w * sim_fp) + ((1.0 - drfp_w) * sim_cat)
                else:
                    sim_total = sim_fp

        # ── MixFP cosine score blend ──────────────────────────────────────
        if mixfp_scores:
            r_id = r.get("reaction_id")
            sim_mixfp = mixfp_scores.get(r_id, 0.0) if r_id else 0.0
            sim_mixfp = max(0.0, min(1.0, sim_mixfp))
            if sim_mixfp > 0.0:
                # Three-way blend: mixfp | drfp | categorical
                # mixfp_w controls the MixFP share; remainder split evenly
                _rem = (1.0 - mixfp_w) / 2.0
                sim_total = (
                    mixfp_w * sim_mixfp
                    + _rem * (sim_fp if q_fp is not None and r_fp is not None else sim_cat)
                    + _rem * sim_cat
                )

        if sim_total <= 0:
            if sim_cat > 0:
                sim_total = sim_cat
            else:
                continue
        neighbor_score = sim_total * (0.5 + 0.5 * y_norm)
        scored.append((neighbor_score, r))

    if not scored:
        proto = f"proto_{_proto_family_id(family_txt)}_none_0"
        return {"prototype_id": proto, "support": 0, "precedents": [], "error": "NO_PRECEDENTS"}

    scored.sort(key=lambda x: (-(x[0]), -((x[1].get("yield_value") or 0))))
    timing['scoring'] = time.time() - t_start_scoring
    
    # ⏱️ TIME: Result preparation
    t_start_prep = time.time()
    return_limit = int(relax.get("precedent_limit", 10))
    if return_limit <= 0:
        return_limit = max(1, int(k))

    # Keep enough scored rows so output limiting can exceed k when requested.
    top = scored[: max(1, max(int(k), return_limit))]
    support = len(scored)

    # Prototype id is a stable-ish hash of family+bin
    family_norm = _proto_family_id(family_txt) if family_txt else "all_families"
    bin_key = str(features.get("bin") or f"LG:{target_feat.get('LG','?')}|NUC:{target_feat.get('nuc_class','?')}")
    proto = f"proto_{family_norm}_{abs(hash(bin_key)) % 100000}"

    # Include reaction SMILES and parsed sides for UI/consumers
    try:
        from ..smiles import _split_reaction_smiles as _split_rxn  # type: ignore
    except Exception:
        def _split_rxn(rsmi: str):
            parts = (rsmi or "").split(">")
            if len(parts) == 2 and ">>" in (rsmi or ""):
                return parts[0], "", parts[1]
            if len(parts) == 3:
                return parts[0], parts[1], parts[2]
            return rsmi, "", ""

    precedents = []
    for score, r in top[:return_limit]:
        rsmi = r.get("reaction_smiles") or ""
        reactants_smi, _agents_smi, products_smi = _split_rxn(rsmi)
        precedents.append({
            "reaction_id": r.get("reaction_id"),
            "dataset_reaction_id": r.get("dataset_reaction_id"),
            "reaction_smiles": rsmi,
            "reactants_smiles": reactants_smi,
            "products_smiles": products_smi,
            "similarity": round(score, 4),  # Add similarity score from scored tuple
            "condition_core": r.get("condition_core"),
            "yield": r.get("yield_value"),
            "core": r.get("condition_core"),
            "base_uid": r.get("base_uid"),
            "solvent_uid": r.get("solvent_uid"),
            "rxn_type": r.get("rxn_type"),  # Include reaction family/type for dataset identification
            "reagents": r.get("reagents"),
            "solvents": r.get("solvents"),
            "source_file": r.get("source_file"),
            "source_group": r.get("source_group"),
            "reference": r.get("reference"),
            "conditions": r.get("conditions"),
            "catalyst": r.get("catalyst"),
            "full_system": r.get("full_system"),
            "catalytic_system": r.get("catalytic_system"),
            "T_C": r.get("T_C"),
            "time_h": r.get("time_h"),
        })
    
    # Optional filtering: only keep precedents where all reagents are in database
    filter_by_database = relax.get("filter_by_reagent_database", True)  # Default: True
    if filter_by_database:
        try:
            from ..reagent import filter_precedents_by_database_availability
            
            # Get full precedent data for filtering (top contains (score, row) tuples)
            precedents_for_filter = [r for score, r in top]
            filtered = filter_precedents_by_database_availability(
                precedents_for_filter,
                require_all_in_database=True
            )
            
            # Create a mapping of reaction_id -> score for filtered results
            filtered_ids = {r.get("reaction_id") for r in filtered}
            filtered_with_scores = [(score, r) for score, r in top if r.get("reaction_id") in filtered_ids]
            
            # Rebuild precedents list from filtered results with scores
            precedents = []
            for score, r in filtered_with_scores[:return_limit]:
                rsmi = r.get("reaction_smiles") or ""
                reactants_smi, _agents_smi, products_smi = _split_rxn(rsmi)
                precedents.append({
                    "reaction_id": r.get("reaction_id"),
                    "dataset_reaction_id": r.get("dataset_reaction_id"),
                    "reaction_smiles": rsmi,
                    "reactants_smiles": reactants_smi,
                    "products_smiles": products_smi,
                    "similarity": round(score, 4),  # Add similarity score
                    "condition_core": r.get("condition_core"),
                    "yield": r.get("yield_value"),
                    "core": r.get("condition_core"),
                    "base_uid": r.get("base_uid"),
                    "solvent_uid": r.get("solvent_uid"),
                    "rxn_type": r.get("rxn_type"),  # Include reaction family/type for dataset identification
                    "reagents": r.get("reagents"),
                    "solvents": r.get("solvents"),
                    "source_file": r.get("source_file"),
                    "source_group": r.get("source_group"),
                    "reference": r.get("reference"),
                    "conditions": r.get("conditions"),
                    "catalyst": r.get("catalyst"),
                    "full_system": r.get("full_system"),
                    "catalytic_system": r.get("catalytic_system"),
                    "T_C": r.get("T_C"),
                    "time_h": r.get("time_h"),
                })
            
            # Update support to reflect filtered count
            support = len(filtered)
            
            # Log filtering results if debugging
            if relax.get("debug_timing", False):
                original_count = len(top)
                filtered_count = len(filtered)
                logger.info(f"   Reagent database filtering: {original_count} → {filtered_count} precedents")
        except ImportError:
            # If reagent_lookup module not available, skip filtering
            pass
    
    timing['result_prep'] = time.time() - t_start_prep
    
    # ⏱️ TOTAL TIME
    timing['total'] = time.time() - t_start_total
    
    # Log timing information if enabled
    if relax.get("debug_timing", False):
        logger.info(f"⏱️  Precedent search timing breakdown:")
        logger.info(f"   - Data loading: {timing['load_data']:.3f}s")
        logger.info(f"   - Candidate pool: {timing['build_candidates']:.3f}s")
        logger.info(f"   - DRFP query encode: {timing['drfp_query_encode']:.3f}s")
        logger.info(f"   - Similarity scoring: {timing['scoring']:.3f}s")
        logger.info(f"   - Result preparation: {timing['result_prep']:.3f}s")
        logger.info(f"   - Total: {timing['total']:.3f}s")
        if use_drfp:
            logger.info(f"   DRFP loading strategy usage:")
            logger.info(f"   - Binary NPZ: {drfp_load_count['binary']} fingerprints")
            logger.info(f"   - JSONL embedded: {drfp_load_count['jsonl']} fingerprints")
            logger.info(f"   - On-demand compute: {drfp_load_count['computed']} fingerprints")
    
    # Always include timing in result for monitoring
    result = {
        "prototype_id": proto, 
        "support": support, 
        "precedents": precedents,
        "timing": timing
    }
    
    if use_drfp:
        result["drfp_load_strategy"] = drfp_load_count
        result["drfp_rerank_limit"] = on_demand_drfp_limit
        result["drfp_rerank_candidates"] = len(candidate_rows)
        result["drfp_loader_available"] = drfp_loader is not None

    if use_mixfp:
        result["mixfp_routing"] = {
            "candidates": len(mixfp_scores),
            "routing_only": mixfp_routing_only,
            "weight": mixfp_w,
        }

    return result


def find_reactions_by_core(
    core_query: str,
    *,
    family: Optional[str] = None,
    fuzzy: bool = True,
    limit: int = 50,
) -> List[Dict[str, Any]]:
    """Find dataset reactions that use the same or similar condition core.

    - core_query: e.g., 'Pd/XPhos', 'Pd', or 'XPhos'.
    - family: optional reaction family text (e.g., 'Ullmann C–N', 'Suzuki_CC').
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

    - family: optional reaction family text to filter (e.g., 'Ullmann C–N', 'Suzuki_CC').
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
