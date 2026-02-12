"""
Aggregate feature extraction across reaction components.

Calculates reaction-level statistics from reactant/product features.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
import re

from ..spectator_rank import rank_spectator_groups

from .utils import extract_scores, group_id_from_motif_id, normalize_motif_id


_DEFAULT_SPECTATOR_GROUP_STOPLIST = {
    "Ar", "R", "Any", "Alkyl", "Alkenyl", "Alkynyl", "H",
}


@lru_cache(maxsize=1)
def _load_group_sets() -> Dict[str, Any]:
    try:
        from chemtools.taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_group_logic()
    if not payload:
        return {}
    group_sets = payload.get("group_sets", {}) or {}
    return group_sets if isinstance(group_sets, dict) else {}


@lru_cache(maxsize=1)
def _load_spectator_group_stoplist() -> Set[str]:
    try:
        from chemtools.taxonomy import loader as taxonomy_loader
    except Exception:
        return set(_DEFAULT_SPECTATOR_GROUP_STOPLIST)
    payload = taxonomy_loader.load_featurizer_logic()
    if not payload:
        return set(_DEFAULT_SPECTATOR_GROUP_STOPLIST)
    aggregation = payload.get("aggregation", {}) or {}
    if not isinstance(aggregation, dict):
        return set(_DEFAULT_SPECTATOR_GROUP_STOPLIST)
    stoplist = aggregation.get("spectator_group_stoplist", [])
    if not isinstance(stoplist, list):
        return set(_DEFAULT_SPECTATOR_GROUP_STOPLIST)
    configured = {str(v) for v in stoplist if isinstance(v, str) and v.strip()}
    return configured or set(_DEFAULT_SPECTATOR_GROUP_STOPLIST)


def _expand_group_set(
    set_id: str,
    group_sets: Dict[str, Any],
    *,
    _seen: Optional[Set[str]] = None,
) -> Set[str]:
    if not group_sets or not set_id:
        return set()
    if _seen is None:
        _seen = set()
    if set_id in _seen:
        return set()
    _seen.add(set_id)
    set_data = group_sets.get(set_id, {}) or {}
    members = set_data.get("members", []) or []
    expanded: Set[str] = set()
    for member in members:
        if not member:
            continue
        member_str = str(member)
        if member_str in group_sets:
            expanded.update(_expand_group_set(member_str, group_sets, _seen=_seen))
        else:
            expanded.add(member_str)
    return expanded


@lru_cache(maxsize=1)
def _load_carbonyl_groups() -> Set[str]:
    group_sets = _load_group_sets()
    if not group_sets:
        return set()
    return _expand_group_set("Carbonyl_Group", group_sets)


@lru_cache(maxsize=1)
def _load_broad_definition_groups() -> Set[str]:
    """Load group IDs flagged as broad_definition from organic group taxonomy."""
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "organic_groups.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    groups = payload.get("groups", []) or []
    if not isinstance(groups, list):
        return set()
    broad: Set[str] = set()
    for entry in groups:
        if not isinstance(entry, dict):
            continue
        if not bool(entry.get("broad_definition")):
            continue
        gid = str(entry.get("id") or "").strip()
        if gid:
            broad.add(gid)
    return broad


def _is_broad_definition_motif_id(motif_id: str) -> bool:
    motif_text = str(motif_id or "").strip()
    if not motif_text:
        return False
    for group_id in _load_broad_definition_groups():
        if motif_text.endswith(group_id):
            return True
    return False


def _broad_fingerprint_changed_ids(
    reactant_motifs: List[Dict[str, Any]],
    product_motifs: List[Dict[str, Any]],
) -> Set[str]:
    """
    Return motif IDs with broad_definition=true whose fingerprint profile changed.

    This catches site-state changes where motif ID counts stay constant across
    reactant/product sides (e.g., Ar-Hydrazine still present but transformed).
    """
    reactant_fp: Dict[str, Counter[str]] = {}
    product_fp: Dict[str, Counter[str]] = {}

    for motif in reactant_motifs:
        cid = normalize_motif_id(str(motif.get("compound_id") or motif.get("id") or ""))
        if not cid or not _is_broad_definition_motif_id(cid):
            continue
        fp = str(motif.get("fingerprint", "") or "")
        reactant_fp.setdefault(cid, Counter())[fp] += 1

    for motif in product_motifs:
        cid = normalize_motif_id(str(motif.get("compound_id") or motif.get("id") or ""))
        if not cid or not _is_broad_definition_motif_id(cid):
            continue
        fp = str(motif.get("fingerprint", "") or "")
        product_fp.setdefault(cid, Counter())[fp] += 1

    changed: Set[str] = set()
    for cid in set(reactant_fp) & set(product_fp):
        if reactant_fp[cid] != product_fp[cid]:
            changed.add(cid)
    return changed


@lru_cache(maxsize=1)
def load_transformation_patterns() -> Dict[str, Any]:
    """Load transformation patterns from taxonomy."""
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "transformation_patterns.json"
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return json.load(handle)
    except Exception:
        return {}


def get_substituent(motif_id: str) -> str:
    """Extract substituent part from motif ID (e.g., 'Ar-Br' -> '-Br')."""
    if "-" in motif_id:
        parts = motif_id.split("-", 1)
        return "-" + parts[1] if len(parts) > 1 else ""
    return ""


def get_scaffold(motif_id: str) -> str:
    """Extract scaffold part from motif ID (e.g., 'Ar-Br' -> 'Ar')."""
    if "-" in motif_id:
        return motif_id.split("-", 1)[0]


def make_motif_key(motif: Dict[str, Any]) -> str:
    """
    Create a unique key for a motif that includes fingerprint info.
    
    This allows distinguishing between the same compound_id with different
    substitution states (e.g., Ar-Hydrazide with NH2 vs Ar-Hydrazide with N=).
    
    Key format: "{compound_id}|{fingerprint}" or just "{compound_id}" if no fingerprint.
    """
    compound_id = motif.get("compound_id") or motif.get("id", "")
    fingerprint = motif.get("fingerprint", "")
    if fingerprint:
        return f"{compound_id}|{fingerprint}"
    return compound_id


def analyze_motif_changes_with_fingerprints(
    reactant_motifs: List[Dict[str, Any]],
    product_motifs: List[Dict[str, Any]],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Analyze motif changes using fingerprints to detect substitution state changes.
    
    This is the fingerprint-aware version of analyze_substituent_changes.
    It can detect when a functional group has changed (e.g., hydrazide NH2 -> hydrazone)
    even though the core SMARTS pattern still matches.
    
    Args:
        reactant_motifs: List of motif dicts from reactants (with fingerprint field)
        product_motifs: List of motif dicts from products (with fingerprint field)
        
    Returns:
        Tuple of (reacted_motifs, formed_motifs, spectator_motifs) as compound_id sets
    """
    # Build keyed representations with fingerprints
    reactant_keys = [make_motif_key(m) for m in reactant_motifs]
    product_keys = [make_motif_key(m) for m in product_motifs]
    
    reactant_key_counts = Counter(reactant_keys)
    product_key_counts = Counter(product_keys)
    
    # Map keys back to compound_ids
    key_to_id: Dict[str, str] = {}
    for m in reactant_motifs + product_motifs:
        key = make_motif_key(m)
        compound_id = m.get("compound_id") or m.get("id", "")
        key_to_id[key] = compound_id
    
    reacted_set: Set[str] = set()
    formed_set: Set[str] = set()
    spectator_set: Set[str] = set()
    
    # Compare by fingerprinted keys
    all_keys = set(reactant_key_counts) | set(product_key_counts)
    for key in all_keys:
        rc = reactant_key_counts.get(key, 0)
        pc = product_key_counts.get(key, 0)
        compound_id = key_to_id.get(key, key.split("|")[0])
        
        if pc > rc:
            formed_set.add(compound_id)
            if rc > 0:
                spectator_set.add(compound_id)
        elif pc < rc:
            reacted_set.add(compound_id)
            if pc > 0:
                spectator_set.add(compound_id)
        else:
            if rc > 0:
                spectator_set.add(compound_id)
    
    # Additional check: same compound_id but different fingerprints = change
    # Group by base compound_id
    reactant_by_id: Dict[str, List[str]] = {}
    for m in reactant_motifs:
        cid = m.get("compound_id") or m.get("id", "")
        fp = m.get("fingerprint", "")
        reactant_by_id.setdefault(cid, []).append(fp)
    
    product_by_id: Dict[str, List[str]] = {}
    for m in product_motifs:
        cid = m.get("compound_id") or m.get("id", "")
        fp = m.get("fingerprint", "")
        product_by_id.setdefault(cid, []).append(fp)
    
    # Check for fingerprint changes within same compound_id
    for cid in set(reactant_by_id) & set(product_by_id):
        r_fps = set(reactant_by_id[cid])
        p_fps = set(product_by_id[cid])
        
        # If fingerprints differ, this group participated in reaction
        if r_fps != p_fps:
            # Remove from spectator, add to both reacted and formed
            spectator_set.discard(cid)
            if r_fps - p_fps:  # Some fingerprints disappeared
                reacted_set.add(cid)
            if p_fps - r_fps:  # Some fingerprints appeared
                formed_set.add(cid)

    # Stable carbonyl groups: if counts match, treat as spectators even if fingerprints differ.
    carbonyl_groups = _load_carbonyl_groups()
    if carbonyl_groups:
        reactant_id_counts = Counter(
            m.get("compound_id") or m.get("id", "") for m in reactant_motifs
        )
        product_id_counts = Counter(
            m.get("compound_id") or m.get("id", "") for m in product_motifs
        )
        for cid in set(reactant_id_counts) & set(product_id_counts):
            if reactant_id_counts[cid] != product_id_counts[cid]:
                continue
            if any(str(cid).endswith(group_id) for group_id in carbonyl_groups):
                reacted_set.discard(cid)
                formed_set.discard(cid)
                spectator_set.add(cid)

    # Substituent-centric override (scaffold changes): if a substituent disappears
    # globally, ensure motifs with that substituent are marked as reacted.
    reactant_ids = [m.get("compound_id") or m.get("id", "") for m in reactant_motifs]
    product_ids = [m.get("compound_id") or m.get("id", "") for m in product_motifs]
    reactant_subs: Dict[str, List[str]] = {}
    for motif_id in reactant_ids:
        if not motif_id:
            continue
        sub = get_substituent(motif_id)
        if sub:
            reactant_subs.setdefault(sub, []).append(motif_id)
    product_subs = {get_substituent(m) for m in product_ids if m}
    consumed_subs = set(reactant_subs.keys()) - product_subs
    if consumed_subs:
        for sub in consumed_subs:
            for motif_id in reactant_subs.get(sub, []):
                reacted_set.add(motif_id)
                spectator_set.discard(motif_id)

    return reacted_set, formed_set, spectator_set


def select_primary_motifs_by_atom(motif_dicts: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Select primary motif for each attachment atom (a_atom_idx).
    
    When multiple motifs share the same attachment atom (e.g., R_acidic-CN and 
    R_acidic-OH on same carbon), only one can be the primary reactive site.
    Uses reactivity_weight to select, falls back to priority.
    
    Args:
        motif_dicts: List of motif dictionaries with a_idx (or a_atom_idx), reactivity_weight, etc.
        
    Returns:
        Filtered list with one primary motif per attachment atom
    """
    if not motif_dicts:
        return []
    
    # Group by attachment atom (supports both a_idx and a_atom_idx keys)
    by_atom: Dict[int, List[Dict[str, Any]]] = {}
    for m in motif_dicts:
        a_idx = m.get("a_idx") or m.get("a_atom_idx")
        if a_idx is None:
            # No atom info, keep as-is
            continue
        by_atom.setdefault(a_idx, []).append(m)
    
    # Select best motif per atom
    selected: List[Dict[str, Any]] = []
    for a_idx, motifs in by_atom.items():
        if len(motifs) == 1:
            selected.append(motifs[0])
        else:
            # Sort by reactivity_weight (desc), then priority (desc)
            best = max(
                motifs, 
                key=lambda m: (m.get("reactivity_weight", 0), m.get("priority", 0))
            )
            selected.append(best)
    
    # Add any motifs without atom info
    for m in motif_dicts:
        a_idx = m.get("a_idx") or m.get("a_atom_idx")
        if a_idx is None:
            selected.append(m)
    
    return selected


def extract_motif_with_bond_info(motif: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract motif info preserving bond-level data and fingerprint.
    
    Returns dict with:
        - id: Normalized compound ID
        - a_idx: Attachment atom index
        - b_idx: Substituent atom index  
        - bond: Site bond tuple
        - reactivity_weight: For primary site selection
        - fingerprint: Heteroatom H-count/connectivity signature (for change detection)
        - h_count: Total H count on heteroatoms (quick comparison)
    """
    cid = motif.get("compound_id") or motif.get("id")
    if not cid:
        return {}
    return {
        "id": normalize_motif_id(str(cid)),
        "compound_id": normalize_motif_id(str(cid)),  # Alias for fingerprint functions
        "a_idx": motif.get("a_atom_idx"),
        "b_idx": motif.get("b_atom_idx"),
        "bond": motif.get("bond"),
        "reactivity_weight": motif.get("reactivity_weight", 0),
        "priority": motif.get("priority", 0),
        "fingerprint": motif.get("fingerprint", ""),
        "h_count": motif.get("h_count", 0),
    }


def analyze_motif_changes_with_bonds(
    reactant_motifs: List[Dict[str, Any]],
    product_motifs: List[Dict[str, Any]],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Analyze motif changes using bond-level tracking.
    
    Compares motifs by both ID and bond position when available.
    Falls back to ID-only comparison when bond info missing.
    
    Args:
        reactant_motifs: List of motif dicts with id, a_idx, b_idx, bond
        product_motifs: List of motif dicts with id, a_idx, b_idx, bond
        
    Returns:
        Tuple of (reacted_ids, formed_ids, spectator_ids)
    """
    # Build reactant index: (id, bond) -> motif
    reactant_by_id: Dict[str, List[Dict]] = {}
    reactant_by_bond: Dict[Tuple, Dict] = {}
    for m in reactant_motifs:
        mid = m.get("id", "")
        if mid:
            reactant_by_id.setdefault(mid, []).append(m)
        bond = m.get("bond")
        if bond:
            reactant_by_bond[tuple(sorted(bond))] = m
    
    # Build product index
    product_by_id: Dict[str, List[Dict]] = {}
    product_by_bond: Dict[Tuple, Dict] = {}
    for m in product_motifs:
        mid = m.get("id", "")
        if mid:
            product_by_id.setdefault(mid, []).append(m)
        bond = m.get("bond")
        if bond:
            product_by_bond[tuple(sorted(bond))] = m
    
    # Analyze changes
    all_ids = set(reactant_by_id.keys()) | set(product_by_id.keys())
    reacted_ids: Set[str] = set()
    formed_ids: Set[str] = set()
    spectator_ids: Set[str] = set()
    
    for mid in all_ids:
        r_count = len(reactant_by_id.get(mid, []))
        p_count = len(product_by_id.get(mid, []))
        
        if p_count > r_count:
            formed_ids.add(mid)
            if r_count > 0:
                spectator_ids.add(mid)
        elif p_count < r_count:
            reacted_ids.add(mid)
            if p_count > 0:
                spectator_ids.add(mid)
        else:
            if r_count > 0:
                spectator_ids.add(mid)
    
    # Bond-level refinement: if bond still exists in product, motif didn't react
    for mid in list(reacted_ids):
        r_motifs = reactant_by_id.get(mid, [])
        for rm in r_motifs:
            bond = rm.get("bond")
            if bond and tuple(sorted(bond)) in product_by_bond:
                # Bond still exists in product - this specific instance is spectator
                # But if ANY instance reacted, keep in reacted
                pass  # Keep current classification for now
    
    # Apply substituent-centric refinement
    return analyze_substituent_changes(
        [m.get("id", "") for m in reactant_motifs if m.get("id")],
        [m.get("id", "") for m in product_motifs if m.get("id")]
    )


def analyze_substituent_changes(
    reactant_motif_ids: List[str],
    product_motif_ids: List[str],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Analyze motif changes using substituent-centric comparison.
    
    This approach focuses on which substituents disappear/appear,
    working across different scaffold types (Ar, HeteroAr, Alkenyl, etc.).
    
    Args:
        reactant_motif_ids: List of motif IDs from reactants
        product_motif_ids: List of motif IDs from products
        
    Returns:
        Tuple of (reacted_motifs, formed_motifs, spectator_motifs)
    """
    reactant_counts = Counter(reactant_motif_ids)
    product_counts = Counter(product_motif_ids)
    
    # Build substituent -> motifs mapping for both sides
    reactant_subs: Dict[str, List[str]] = {}
    for motif_id in reactant_motif_ids:
        sub = get_substituent(motif_id)
        if sub:
            reactant_subs.setdefault(sub, []).append(motif_id)
    
    product_subs: Dict[str, List[str]] = {}
    for motif_id in product_motif_ids:
        sub = get_substituent(motif_id)
        if sub:
            product_subs.setdefault(sub, []).append(motif_id)
    
    reacted_set: Set[str] = set()
    formed_set: Set[str] = set()
    spectator_set: Set[str] = set()
    
    # Analyze each motif using count comparison
    for motif_id in set(reactant_counts) | set(product_counts):
        rc = reactant_counts.get(motif_id, 0)
        pc = product_counts.get(motif_id, 0)
        
        if pc > rc:
            formed_set.add(motif_id)
            if rc > 0:
                spectator_set.add(motif_id)
        elif pc < rc:
            reacted_set.add(motif_id)
            if pc > 0:
                spectator_set.add(motif_id)
        else:
            if rc > 0:
                spectator_set.add(motif_id)
    
    # Cross-scaffold substituent matching:
    # If substituent appears in products on ANY scaffold, motifs with that 
    # substituent might be spectators, not truly reacted
    consumed_subs = set(reactant_subs.keys()) - set(product_subs.keys())
    appeared_subs = set(product_subs.keys()) - set(reactant_subs.keys())
    
    # Refine reacted: prioritize motifs with truly consumed substituents
    refined_reacted: Set[str] = set()
    for motif_id in reacted_set:
        sub = get_substituent(motif_id)
        # If substituent is in consumed list, definitely reacted
        if sub in consumed_subs:
            refined_reacted.add(motif_id)
        # If substituent still exists on other scaffolds, might be spectator
        elif sub in product_subs:
            # Check count - if product count is lower overall, still reacted
            sub_rc = len(reactant_subs.get(sub, []))
            sub_pc = len(product_subs.get(sub, []))
            if sub_pc < sub_rc:
                refined_reacted.add(motif_id)
            else:
                spectator_set.add(motif_id)
        else:
            # No substituent info, keep original classification
            refined_reacted.add(motif_id)
    
    # Refine formed: prioritize motifs with truly new substituents
    refined_formed: Set[str] = set()
    for motif_id in formed_set:
        sub = get_substituent(motif_id)
        if sub in appeared_subs:
            refined_formed.add(motif_id)
        elif sub in reactant_subs:
            sub_rc = len(reactant_subs.get(sub, []))
            sub_pc = len(product_subs.get(sub, []))
            if sub_pc > sub_rc:
                refined_formed.add(motif_id)
            else:
                spectator_set.add(motif_id)
        else:
            refined_formed.add(motif_id)
    
    return refined_reacted, refined_formed, spectator_set


def filter_reacted_by_pattern(
    reacted_motifs: List[str],
    reaction_type: Optional[str],
) -> List[str]:
    """
    Filter reacted motifs using transformation pattern rules.
    
    Uses the unified 'leaving_groups' field from transformation_patterns.json.
    Only keeps motifs whose substituent matches expected leaving groups.
    
    Args:
        reacted_motifs: Raw list of reacted motif IDs
        reaction_type: Detected reaction type ID
        
    Returns:
        Filtered list of reacted motifs relevant to the reaction
    """
    if not reaction_type or not reacted_motifs:
        return reacted_motifs
    
    patterns = load_transformation_patterns()
    reaction_patterns = patterns.get("reaction_patterns", {})
    
    pattern_info = reaction_patterns.get(reaction_type)
    if not pattern_info:
        return reacted_motifs
    
    # All patterns now use unified 'leaving_groups' field
    leaving_groups: Set[str] = set(pattern_info.get("leaving_groups", []))
    
    if not leaving_groups:
        # No leaving groups defined (e.g., addition reactions) - return original
        return reacted_motifs
    
    # Filter: keep only motifs with substituents in leaving_groups
    filtered = []
    for motif in reacted_motifs:
        sub = get_substituent(motif)
        if sub in leaving_groups:
            filtered.append(motif)
    
    # If filtering removed everything, return original (fallback)
    return filtered if filtered else reacted_motifs


@lru_cache(maxsize=1)
def load_scaffold_motif_ids() -> set[str]:
    """Load scaffold motif IDs from taxonomy."""
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "scaffold_motifs.v1.3.json"
    if not path.exists():
        return set()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return set()
    motifs: set[str] = set()
    for entry in payload.get("compounds", []) or []:
        if not isinstance(entry, dict):
            continue
        motif_id = str(entry.get("id") or "").strip()
        if motif_id:
            motifs.add(motif_id)
    return motifs


def aggregate_reaction_features(
    reactants: Iterable[Dict[str, Any]],
    *,
    product_motif_ids: Optional[List[str]] = None,
    product_motifs: Optional[List[Dict[str, Any]]] = None,
    reaction_type: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Aggregate features across all reactants and products.
    
    Args:
        reactants: List of reactant feature bundles
        product_motif_ids: Optional list of product motif IDs for change analysis
        product_motifs: Optional list of product motif dicts with fingerprints
                       (for fingerprint-aware change detection)
        reaction_type: Optional reaction type ID for pattern-based filtering
        
    Returns:
        Aggregated features including motifs, spectators, steric/electronic stats
    """
    reactant_list = list(reactants)
    aryl_scores: List[float] = []
    alkyl_scores: List[float] = []
    electronic_scores: List[float] = []
    motifs: set[str] = set()
    reactant_motif_ids: List[str] = []
    reactant_motifs_per_mol: List[List[Dict[str, Any]]] = []  # Keep separate per reactant to avoid atom index collisions
    spectator_groups: List[str] = []
    spectator_seen: set[str] = set()
    scaffold_ids = load_scaffold_motif_ids()
    reacted_motifs: List[str] = []
    formed_motifs: List[str] = []
    spectator_motifs: List[str] = []
    reacted_motif_counts: Dict[str, int] = {}
    formed_motif_counts: Dict[str, int] = {}
    spectator_motif_counts: Dict[str, int] = {}
    broad_fp_changed_ids: Set[str] = set()

    # Extract features from each reactant
    for reactant in reactant_list:
        for entry in reactant.get("steric", {}).get("aryl", []):
            aryl_scores.extend(extract_scores(entry.get("result")))
        for entry in reactant.get("steric", {}).get("alkyl", []):
            alkyl_scores.extend(extract_scores(entry.get("result")))
        for entry in reactant.get("electronics", {}).get("aryl", []):
            electronic_scores.extend(extract_scores(entry.get("result")))
        
        motif_entries = reactant.get("motifs", [])
        context_entries = reactant.get("context_motifs", [])
        for motif in motif_entries:
            compound_id = normalize_motif_id(motif.get("compound_id") or motif.get("id") or "")
            if compound_id:
                motifs.add(str(compound_id))
        
        # Phase 3: Collect full motif dicts with bond info (per reactant to avoid atom index collisions)
        motifs_for_this_reactant: List[Dict[str, Any]] = []
        for motif in motif_entries:
            if isinstance(motif, dict):
                cid = motif.get("compound_id") or motif.get("id")
                if cid:
                    reactant_motif_ids.append(normalize_motif_id(str(cid)))
                    motif_info = extract_motif_with_bond_info(motif)
                    if motif_info:
                        motifs_for_this_reactant.append(motif_info)
        if context_entries:
            for motif in context_entries:
                if isinstance(motif, dict):
                    cid = motif.get("compound_id") or motif.get("id")
                    if cid:
                        reactant_motif_ids.append(normalize_motif_id(str(cid)))
                        motif_info = extract_motif_with_bond_info(motif)
                        if motif_info:
                            motifs_for_this_reactant.append(motif_info)
        
        if motifs_for_this_reactant:
            reactant_motifs_per_mol.append(motifs_for_this_reactant)
    
    # Phase 3: Select primary motif per attachment atom WITHIN EACH REACTANT
    # This prevents atom index collisions between different reactants (e.g., piperazine atom 3 vs aniline atom 3)
    reactant_motifs_full: List[Dict[str, Any]] = []
    reactant_motifs_for_changes: List[Dict[str, Any]] = []
    for motifs_for_mol in reactant_motifs_per_mol:
        primary_for_mol = select_primary_motifs_by_atom(motifs_for_mol)
        reactant_motifs_full.extend(primary_for_mol)
        reactant_motifs_for_changes.extend(motifs_for_mol)
    
    # Extract IDs from primary motifs for use in change analysis
    primary_motif_ids_list = [m.get("id", "") for m in reactant_motifs_full if m.get("id")]

    # Analyze motif changes if products are provided
    if product_motif_ids or product_motifs:
        # Use fingerprint-aware comparison if full motif dicts are available
        if product_motifs:
            # Extract product motifs with fingerprints
            product_motifs_with_fp = [extract_motif_with_bond_info(m) for m in product_motifs if isinstance(m, dict)]
            product_motifs_with_fp = [m for m in product_motifs_with_fp if m.get("id")]
            product_motifs_primary = select_primary_motifs_by_atom(product_motifs_with_fp)
            
            reacted_set, formed_set, spectator_motifs_set = analyze_motif_changes_with_fingerprints(
                reactant_motifs_for_changes, product_motifs_with_fp
            )
            broad_fp_changed_ids = _broad_fingerprint_changed_ids(
                reactant_motifs_for_changes,
                product_motifs_with_fp,
            )
            if broad_fp_changed_ids:
                reacted_set.update(broad_fp_changed_ids)
                formed_set.update(broad_fp_changed_ids)
                spectator_motifs_set.difference_update(broad_fp_changed_ids)
            reactant_ids_for_counts = [
                m.get("compound_id") or m.get("id", "") for m in reactant_motifs_full
            ]
            product_ids_for_counts = [
                m.get("compound_id") or m.get("id", "") for m in product_motifs_primary
            ]
        else:
            # Fallback to ID-only comparison
            product_ids = product_motif_ids or [m.get("compound_id") or m.get("id", "") for m in (product_motifs or [])]
            reacted_set, formed_set, spectator_motifs_set = analyze_substituent_changes(
                primary_motif_ids_list, product_ids
            )
            reactant_ids_for_counts = primary_motif_ids_list
            product_ids_for_counts = product_ids
        
        # Collect spectator groups
        for motif_id in primary_motif_ids_list:
            if motif_id not in spectator_motifs_set:
                continue
            group_id = group_id_from_motif_id(motif_id)
            if not group_id or group_id in _load_spectator_group_stoplist():
                continue
            if group_id in spectator_seen:
                continue
            spectator_seen.add(group_id)
            spectator_groups.append(group_id)

        # Add scaffold motifs as spectators
        for motif_id in spectator_motifs_set:
            if motif_id in scaffold_ids and motif_id not in spectator_seen:
                spectator_seen.add(motif_id)
                spectator_groups.append(motif_id)
        
        # Apply transformation pattern filtering to reacted motifs
        raw_reacted = sorted(reacted_set)
        reacted_motifs = filter_reacted_by_pattern(raw_reacted, reaction_type)
        formed_motifs = sorted(formed_set)
        spectator_motifs = sorted(spectator_motifs_set)

        # Build count-aware views (supports cases where some motifs react and some remain)
        reactant_id_counts = Counter(
            normalize_motif_id(str(mid)) for mid in reactant_ids_for_counts if mid
        )
        product_id_counts = Counter(
            normalize_motif_id(str(mid)) for mid in product_ids_for_counts if mid
        )
        for motif_id in set(reactant_id_counts) | set(product_id_counts):
            r_count = reactant_id_counts.get(motif_id, 0)
            p_count = product_id_counts.get(motif_id, 0)
            if motif_id in reacted_set:
                reacted_motif_counts[motif_id] = max(r_count - p_count, 0)
            if motif_id in formed_set:
                formed_motif_counts[motif_id] = max(p_count - r_count, 0)
            if (
                motif_id in broad_fp_changed_ids
                and r_count > 0
                and p_count > 0
                and r_count == p_count
            ):
                reacted_motif_counts[motif_id] = max(reacted_motif_counts.get(motif_id, 0), 1)
            if motif_id in spectator_motifs_set:
                spectator_motif_counts[motif_id] = min(r_count, p_count)

        # Reconcile lists with counts (e.g., motifs that are effectively unchanged).
        def _count_positive(value: int | None) -> bool:
            return value is None or value > 0

        reacted_motifs = [
            m for m in reacted_motifs
            if _count_positive(reacted_motif_counts.get(m))
        ]
        formed_motifs = [
            m for m in formed_motifs
            if _count_positive(formed_motif_counts.get(m))
        ]
        # If counts indicate unchanged motifs, ensure they appear as spectators.
        for motif_id, r_count in reactant_id_counts.items():
            p_count = product_id_counts.get(motif_id, 0)
            if r_count <= 0 or p_count <= 0 or r_count != p_count:
                continue
            if reacted_motif_counts.get(motif_id) == 0 and formed_motif_counts.get(motif_id, 0) == 0:
                if motif_id not in spectator_motif_counts:
                    spectator_motif_counts[motif_id] = r_count
                if motif_id not in spectator_motifs_set:
                    spectator_motifs_set.add(motif_id)
        spectator_motifs = sorted(spectator_motifs_set)

    avg_electronic = None
    if electronic_scores:
        avg_electronic = round(sum(electronic_scores) / len(electronic_scores), 2)
    
    # Phase 3: Extract primary motif IDs (after per-atom selection)
    primary_motif_ids = [m.get("id", "") for m in reactant_motifs_full if m.get("id")]

    return {
        "reactant_count": len(reactant_list),
        "motif_ids": sorted(motifs),
        "primary_motif_ids": sorted(set(primary_motif_ids)),  # Phase 3: After per-atom selection
        "reacted_motifs": reacted_motifs,
        "reacted_motif_counts": reacted_motif_counts,
        "formed_motifs": formed_motifs,
        "formed_motif_counts": formed_motif_counts,
        "spectator_motifs": spectator_motifs,
        "spectator_motif_counts": spectator_motif_counts,
        "spectator_groups_combined": spectator_groups,
        "spectator_groups_ranked": rank_spectator_groups(spectator_groups),
        "max_aryl_steric": max(aryl_scores) if aryl_scores else 0.0,
        "max_alkyl_steric": max(alkyl_scores) if alkyl_scores else 0.0,
        "avg_aryl_electronic": avg_electronic if avg_electronic is not None else 5.0,
        "electron_poor_aryl": any(score > 6.5 for score in electronic_scores),
    }


def count_elements(smiles_list: Iterable[str]) -> Optional[Counter[str]]:
    """Count element occurrences across multiple SMILES strings."""
    from chemtools.util.rdkit_helpers import parse_smiles
    
    combined = Counter()
    for smiles in smiles_list or []:
        mol = parse_smiles(smiles)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            combined[atom.GetSymbol()] += 1
    return combined if combined else None


def calculate_product_atom_gain(
    reactant_smiles: Iterable[str],
    product_smiles: Iterable[str],
) -> Optional[int]:
    """Calculate net atom gain/loss from reactants to products."""
    r_counts = count_elements(reactant_smiles)
    p_counts = count_elements(product_smiles)
    if r_counts is None or p_counts is None:
        return None
    
    r_total = sum(r_counts.values())
    p_total = sum(p_counts.values())
    return p_total - r_total


def infer_intramolecular(
    reactant_smiles: Iterable[str],
    product_smiles: Iterable[str],
    roles_summary: Optional[Dict[str, Any]],
) -> Optional[bool]:
    """
    Infer if a reaction is intramolecular.
    
    Uses multiple heuristics:
    1. Single reactant with roles that interact
    2. Product atom count similar to reactant count
    """
    reactants = list(reactant_smiles or [])
    products = list(product_smiles or [])
    
    if not reactants or not products:
        return None
    
    # Heuristic 1: Single reactant
    if len(reactants) == 1:
        if roles_summary:
            electrophile_count = roles_summary.get("electrophile_count", 0)
            nucleophile_count = roles_summary.get("nucleophile_count", 0)
            if electrophile_count > 0 and nucleophile_count > 0:
                return True
    
    # Heuristic 2: Atom count analysis
    atom_gain = calculate_product_atom_gain(reactants, products)
    if atom_gain is not None:
        # If we lose atoms, likely intramolecular (condensation, etc.)
        if atom_gain < 0:
            return True
        # If atom gain is small relative to reactant size, might be intramolecular
        r_counts = count_elements(reactants)
        if r_counts and abs(atom_gain) <= 5 and len(reactants) == 1:
            return True
    
    return False
