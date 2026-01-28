"""
Aggregate feature extraction across reaction components.

Calculates reaction-level statistics from reactant/product features.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
import re

from ..spectator_rank import rank_spectator_groups

from .utils import extract_scores, group_id_from_motif_id, normalize_motif_id


_SPECTATOR_GROUP_STOPLIST = {
    "Ar", "R", "Any", "Alkyl", "Alkenyl", "Alkynyl", "H",
}


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
) -> Dict[str, Any]:
    """
    Aggregate features across all reactants and products.
    
    Args:
        reactants: List of reactant feature bundles
        product_motif_ids: Optional list of product motif IDs for change analysis
        
    Returns:
        Aggregated features including motifs, spectators, steric/electronic stats
    """
    reactant_list = list(reactants)
    aryl_scores: List[float] = []
    alkyl_scores: List[float] = []
    electronic_scores: List[float] = []
    motifs: set[str] = set()
    reactant_motif_ids: List[str] = []
    spectator_groups: List[str] = []
    spectator_seen: set[str] = set()
    scaffold_ids = load_scaffold_motif_ids()
    reacted_motifs: List[str] = []
    formed_motifs: List[str] = []
    spectator_motifs: List[str] = []

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
            compound_id = normalize_motif_id(motif.get("compound_id") or "")
            if compound_id:
                motifs.add(str(compound_id))
        
        for motif in motif_entries:
            if isinstance(motif, dict):
                cid = motif.get("compound_id")
                if cid:
                    reactant_motif_ids.append(normalize_motif_id(str(cid)))
        if context_entries:
            for motif in context_entries:
                if isinstance(motif, dict):
                    cid = motif.get("compound_id")
                    if cid:
                        reactant_motif_ids.append(normalize_motif_id(str(cid)))

    # Analyze motif changes if products are provided
    if product_motif_ids:
        reactant_counts = Counter(reactant_motif_ids)
        product_counts = Counter(product_motif_ids)
        reacted_set: set[str] = set()
        formed_set: set[str] = set()
        spectator_motifs_set = {
            motif_id
            for motif_id in reactant_counts
            if product_counts.get(motif_id, 0) > 0
        }
        
        for motif_id in set(reactant_counts) | set(product_counts):
            rc = reactant_counts.get(motif_id, 0)
            pc = product_counts.get(motif_id, 0)
            if pc > rc:
                formed_set.add(motif_id)
                if rc > 0:
                    spectator_motifs_set.add(motif_id)
            elif pc < rc:
                reacted_set.add(motif_id)
                if pc > 0:
                    spectator_motifs_set.add(motif_id)
            else:
                if rc > 0:
                    spectator_motifs_set.add(motif_id)
        
        # Collect spectator groups
        for motif_id in reactant_motif_ids:
            if motif_id not in spectator_motifs_set:
                continue
            group_id = group_id_from_motif_id(motif_id)
            if not group_id or group_id in _SPECTATOR_GROUP_STOPLIST:
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
        
        reacted_motifs = sorted(reacted_set)
        formed_motifs = sorted(formed_set)
        spectator_motifs = sorted(spectator_motifs_set)

    avg_electronic = None
    if electronic_scores:
        avg_electronic = round(sum(electronic_scores) / len(electronic_scores), 2)

    return {
        "reactant_count": len(reactant_list),
        "motif_ids": sorted(motifs),
        "reacted_motifs": reacted_motifs,
        "formed_motifs": formed_motifs,
        "spectator_motifs": spectator_motifs,
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
