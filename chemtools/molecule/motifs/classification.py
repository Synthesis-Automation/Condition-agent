"""Compound classification using motif registry."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from chemtools.core.rdkit import parse_smiles, rdkit_available

from .models import CompoundPattern
from .registry import _default_registry_paths, build_compound_registry
from .utils import calculate_smarts_complexity


def classify_compound_smiles(
    smiles: str,
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
) -> Dict[str, Any]:
    """Classify a SMILES string into compound motif IDs using SMARTS detection.
    
    Detects all compound motifs that match the SMILES and optionally identifies
    the most specific match based on priority and complexity.
    
    Args:
        smiles: SMILES string to classify
        registry: Pre-compiled registry (if None, will build from paths)
        registry_paths: Paths to taxonomy files (uses defaults if None)
        include_best: Whether to compute the best (most specific) match
        
    Returns:
        Dictionary with keys:
        - smiles: Input SMILES string
        - ok: True if classification succeeded
        - hits: List of matching compound IDs
        - best: Most specific compound ID (if include_best=True)
        - error: Error message (if ok=False)
    """
    result: Dict[str, Any] = {"smiles": smiles, "ok": False, "hits": [], "best": None}
    if not rdkit_available():
        result["error"] = "rdkit_unavailable"
        return result

    mol = parse_smiles(smiles)
    if mol is None:
        result["error"] = "invalid_smiles"
        return result

    if registry is None:
        registry_paths = registry_paths or _default_registry_paths()
        registry = build_compound_registry(registry_paths)

    compiled = registry.get("compiled_compounds") or []
    hits: List[str] = []

    # Use the compiled patterns which now include priority
    for pattern in compiled:
        if mol.HasSubstructMatch(pattern.query):
            hits.append(pattern.compound_id)

    result["ok"] = True
    result["hits"] = list(set(hits))
    if include_best:
        best = choose_best_compound_hit(
            result["hits"],
            compound_map=registry.get("compound_map"),
        )
        result["best"] = best
    return result


def classify_compound_batch(
    smiles_list: Iterable[str],
    *,
    registry: Optional[Mapping[str, Any]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
    include_best: bool = True,
) -> List[Dict[str, Any]]:
    """Classify a batch of SMILES strings into compound motif IDs.
    
    More efficient than repeated single calls because registry is built once.
    
    Args:
        smiles_list: Iterable of SMILES strings
        registry: Pre-compiled registry (if None, will build from paths)
        registry_paths: Paths to taxonomy files (uses defaults if None)
        include_best: Whether to compute the best match for each SMILES
        
    Returns:
        List of classification result dictionaries (one per SMILES)
    """
    if registry is None:
        registry_paths = registry_paths or _default_registry_paths()
        registry = build_compound_registry(registry_paths)

    return [
        classify_compound_smiles(
            smiles,
            registry=registry,
            include_best=include_best,
        )
        for smiles in smiles_list
    ]


def choose_best_compound_hit(
    hits: Iterable[str],
    *,
    compound_map: Optional[Mapping[str, CompoundPattern]] = None,
) -> Optional[str]:
    """Choose the most specific motif from multiple matches.
    
    Selection criteria (in order):
    1. Higher priority (scaffold/substituent importance)
    2. Higher complexity (more specific SMARTS pattern)
    3. Longer compound ID (more detailed name)
    4. Alphabetical order (deterministic tie-breaker)
    
    Args:
        hits: List of compound IDs that matched
        compound_map: Mapping of compound_id -> CompoundPattern
        
    Returns:
        Most specific compound ID, or None if hits is empty
    """
    hits_list = [h for h in hits if h]
    if not hits_list:
        return None

    if not compound_map:
        # Fallback to alphabetical if no map provided
        return sorted(hits_list)[0]

    def score_rank(hit: str) -> tuple[int, int, int, str]:
        pattern = compound_map.get(hit)
        priority = pattern.priority if pattern else 0
        complexity = calculate_smarts_complexity(pattern.query) if pattern else 0
        # Higher priority first, then higher complexity (narrower),
        # then longer ID as tie-breaker, then alphabetical
        return (-priority, -complexity, -len(hit), hit)

    return sorted(hits_list, key=score_rank)[0]


__all__ = [
    "classify_compound_smiles",
    "classify_compound_batch",
    "choose_best_compound_hit",
]
