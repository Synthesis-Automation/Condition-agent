"""
Calculable Feature Detection based on calculable_features.json

Detects structural motifs and properties defined in the feature specification:
- Boolean features (e.g., sp2_chloride_present) via SMARTS matching
- Integer counts (e.g., sp2_halide_site_count) via substructure counting
- Heuristic features (e.g., polarity_high) via descriptor calculations
- Derived features (e.g., internal_alkyne_present) via boolean logic

Usage:
    from chemtools.featurizers.calculable import detect_all_features, detect_feature
    
    # Get all features for a molecule
    features = detect_all_features("c1ccc(Br)cc1")
    # {'sp2_bromide_present': True, 'aryl_halide_present': True, 'ArBr_present': True, ...}
    
    # Get a single feature
    has_aryl_bromide = detect_feature("c1ccc(Br)cc1", "sp2_bromide_present")
    # True
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ..util.rdkit_helpers import rdkit_available, parse_smiles

# Path to the feature specification JSON
_SPEC_PATH = Path(__file__).parent / "calculable_features.json"


# ============================================================================
# Feature Specification Loading
# ============================================================================

@lru_cache(maxsize=1)
def _load_feature_spec() -> Dict[str, Any]:
    """Load and cache the feature specification JSON."""
    with open(_SPEC_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def get_feature_spec() -> Dict[str, Any]:
    """
    Return the loaded feature specification JSON.
    
    Returns:
        Dictionary containing version, schema_notes, features, and derived_shortcuts
    """
    return _load_feature_spec()


# ============================================================================
# SMARTS Pattern Compilation and Caching
# ============================================================================

@lru_cache(maxsize=512)
def _compile_smarts(smarts: str):
    """Compile a SMARTS pattern with caching. Returns None if RDKit unavailable or invalid."""
    if not rdkit_available():
        return None
    
    try:
        from rdkit import Chem
        pattern = Chem.MolFromSmarts(smarts)
        return pattern
    except Exception:
        return None


# ============================================================================
# Boolean Feature Detection (SMARTS-based)
# ============================================================================

def _detect_smarts_feature(mol, smarts_list: List[str]) -> bool:
    """
    Return True if ANY SMARTS pattern matches the molecule.
    
    Args:
        mol: RDKit molecule object
        smarts_list: List of SMARTS patterns to check
        
    Returns:
        True if at least one pattern matches, False otherwise
    """
    if mol is None:
        return False
    
    for smarts in smarts_list:
        pattern = _compile_smarts(smarts)
        if pattern is not None:
            try:
                if mol.HasSubstructMatch(pattern):
                    return True
            except Exception:
                continue
    
    return False


# ============================================================================
# Integer Count Features (Substructure Counting)
# ============================================================================

def _count_substructure_matches(mol, smarts_list: List[str]) -> int:
    """
    Count total matches across all SMARTS patterns.
    
    Args:
        mol: RDKit molecule object
        smarts_list: List of SMARTS patterns to count
        
    Returns:
        Total number of substructure matches
    """
    if mol is None:
        return 0
    
    total = 0
    for smarts in smarts_list:
        pattern = _compile_smarts(smarts)
        if pattern is not None:
            try:
                matches = mol.GetSubstructMatches(pattern)
                total += len(matches)
            except Exception:
                continue
    
    return total


# ============================================================================
# Heuristic Feature Detection
# ============================================================================

def _calculate_polarity_features(mol) -> Dict[str, bool]:
    """
    Calculate polarity-based features using molecular descriptors.
    
    Features:
        - polarity_high: TPSA≥50 or (HBA+HBD)≥4
        - polarity_low: TPSA≤20 and (HBA+HBD)≤2
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dictionary with polarity_high and polarity_low boolean values
    """
    if mol is None or not rdkit_available():
        return {'polarity_high': False, 'polarity_low': False}
    
    try:
        from rdkit.Chem import Descriptors
        
        tpsa = Descriptors.TPSA(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        
        polarity_high = tpsa >= 50 or (hba + hbd) >= 4
        polarity_low = tpsa <= 20 and (hba + hbd) <= 2
        
        return {
            'polarity_high': polarity_high,
            'polarity_low': polarity_low
        }
    except Exception:
        return {'polarity_high': False, 'polarity_low': False}


def _detect_beta_hydride(mol) -> bool:
    """
    Detect if molecule has sp3 C-LG with β-C bearing ≥1 H.
    
    This checks for potential β-hydride elimination risk in alkyl halides/leaving groups.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        True if β-hydride is possible
    """
    if mol is None or not rdkit_available():
        return False
    
    try:
        from rdkit import Chem
        
        # Pattern for sp3 C with common leaving groups (Cl, Br, I, OTs, OTf, OMs)
        lg_patterns = [
            "[CX4][Cl,Br,I]",  # Alkyl halides
            "[CX4]OS(=O)(=O)",  # Alkyl sulfonates
        ]
        
        for smarts in lg_patterns:
            pattern = _compile_smarts(smarts)
            if pattern is None:
                continue
            
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # match[0] is the sp3 carbon with LG
                c_idx = match[0]
                c_atom = mol.GetAtomWithIdx(c_idx)
                
                # Check neighbors for β-carbons
                for neighbor in c_atom.GetNeighbors():
                    # Skip the leaving group
                    if neighbor.GetAtomicNum() in (17, 35, 53, 8):  # Cl, Br, I, O
                        continue
                    
                    # Check if this is a carbon
                    if neighbor.GetAtomicNum() == 6:
                        # Check if it has hydrogens
                        if neighbor.GetTotalNumHs() >= 1:
                            return True
        
        return False
    except Exception:
        return False


def _detect_heuristic_features(mol, heuristic_desc: str, token: str) -> Union[bool, int]:
    """
    Detect features that require heuristic/descriptor-based logic.
    
    Args:
        mol: RDKit molecule object
        heuristic_desc: Description of the heuristic from JSON
        token: Feature token name
        
    Returns:
        Boolean or integer value depending on feature type
    """
    # Polarity features
    if token in ('polarity_high', 'polarity_low'):
        polarity_features = _calculate_polarity_features(mol)
        return polarity_features.get(token, False)
    
    # β-hydride detection
    if token == 'beta_hydride_possible':
        return _detect_beta_hydride(mol)
    
    # Count features with SMARTS patterns embedded in heuristic description
    if token == 'sp2_halide_site_count':
        # "count [cX3]-[Cl,Br,I,F] plus [C;X2]=[C;X2]-[Cl,Br,I,F]"
        smarts_list = [
            "[cX3][Cl,Br,I,F]",
            "[C;X2]=[C;X2][Cl,Br,I,F]"
        ]
        return _count_substructure_matches(mol, smarts_list)
    
    if token == 'sp2_sulfonate_site_count':
        # "count sp2-OS(=O)(=O)R sites"
        smarts_list = [
            "cOS(=O)(=O)",
            "[C;X2]=[C;X2]OS(=O)(=O)"
        ]
        return _count_substructure_matches(mol, smarts_list)
    
    # Unknown heuristic - return safe default
    return False if "present" in token else 0


# ============================================================================
# Derived Feature Evaluation
# ============================================================================

def _evaluate_derived_feature(derive_expr: str, base_features: Dict[str, Any]) -> bool:
    """
    Evaluate boolean expressions for derived features.
    
    Supports:
        - AND operations: "feature1 AND feature2"
        - OR operations: "feature1 OR feature2"
        - NOT operations: "NOT feature1"
        - Parentheses: "(feature1 OR feature2) AND feature3"
        - Combinations: "feature1 AND NOT feature2"
    
    Args:
        derive_expr: Expression like "alkyne_present AND NOT terminal_alkyne_present"
        base_features: Dictionary of already-computed features
        
    Returns:
        Boolean result of the expression
    """
    # Normalize expression
    expr = derive_expr.strip()
    
    # Handle parentheses recursively
    while '(' in expr:
        # Find innermost parenthesized expression
        start = expr.rfind('(')
        end = expr.find(')', start)
        if end == -1:
            break
        inner_expr = expr[start+1:end]
        inner_result = _evaluate_derived_feature(inner_expr, base_features)
        expr = expr[:start] + str(inner_result) + expr[end+1:]
    
    # Replace True/False strings with actual booleans
    expr = expr.replace('True', 'true_token').replace('False', 'false_token')
    
    # Split on OR first (lower precedence)
    if ' OR ' in expr:
        or_parts = [part.strip() for part in expr.split(' OR ')]
        result = False
        for part in or_parts:
            part_result = _evaluate_derived_feature(part, base_features)
            result = result or part_result
            if result:  # Short-circuit on true
                return True
        return result
    
    # Split on AND
    and_parts = [part.strip() for part in expr.split(' AND ')]
    
    result = True
    for part in and_parts:
        # Handle boolean tokens from parentheses
        if part == 'true_token':
            continue
        if part == 'false_token':
            return False
            
        # Check for NOT
        if part.startswith('NOT '):
            feature_name = part[4:].strip()
            feature_value = base_features.get(feature_name, False)
            result = result and (not feature_value)
        else:
            feature_name = part
            feature_value = base_features.get(feature_name, False)
            result = result and feature_value
        
        # Short-circuit on false
        if not result:
            return False
    
    return result


# ============================================================================
# Main Detection Functions
# ============================================================================

@lru_cache(maxsize=2048)
def detect_all_features(smiles: str) -> Dict[str, Any]:
    """
    Detect all calculable features for a molecule.
    
    Returns dictionary with boolean, int, and derived features:
    {
        'sp2_chloride_present': True,
        'sp2_halide_site_count': 2,
        'ArBr_present': False,
        'polarity_high': True,
        ...
    }
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Dictionary mapping feature token to value (bool or int)
        
    Examples:
        >>> detect_all_features("c1ccc(Br)cc1")
        {'sp2_bromide_present': True, 'aryl_halide_present': True, 'ArBr_present': True, ...}
        
        >>> detect_all_features("CCBr")
        {'sp3_bromide_present': True, 'sp2_bromide_present': False, ...}
    """
    if not rdkit_available():
        # Return all False/0 when RDKit is not available
        spec = _load_feature_spec()
        result = {}
        for feature in spec.get("features", []):
            token = feature.get("token")
            ftype = feature.get("type", "bool")
            result[token] = False if ftype == "bool" else 0
        for derived in spec.get("derived_shortcuts", []):
            result[derived.get("token")] = False
        return result
    
    mol = parse_smiles(smiles)
    if mol is None:
        # Invalid SMILES - return all False/0
        spec = _load_feature_spec()
        result = {}
        for feature in spec.get("features", []):
            token = feature.get("token")
            ftype = feature.get("type", "bool")
            result[token] = False if ftype == "bool" else 0
        for derived in spec.get("derived_shortcuts", []):
            result[derived.get("token")] = False
        return result
    
    spec = _load_feature_spec()
    result = {}
    
    # Process base features
    for feature in spec.get("features", []):
        token = feature.get("token")
        ftype = feature.get("type", "bool")
        detect = feature.get("detect", {})
        
        # SMARTS-based detection
        if "smarts_any" in detect:
            smarts_list = detect["smarts_any"]
            if ftype == "bool":
                result[token] = _detect_smarts_feature(mol, smarts_list)
            elif ftype == "int":
                result[token] = _count_substructure_matches(mol, smarts_list)
        
        # Heuristic-based detection
        elif "heuristic" in detect:
            heuristic_desc = detect["heuristic"]
            result[token] = _detect_heuristic_features(mol, heuristic_desc, token)
        
        else:
            # Unknown detection method - return safe default
            result[token] = False if ftype == "bool" else 0
    
    # Process derived features
    for derived in spec.get("derived_shortcuts", []):
        token = derived.get("token")
        derive_expr = derived.get("derive", "")
        result[token] = _evaluate_derived_feature(derive_expr, result)
    
    return result


def detect_feature(smiles: str, feature_token: str) -> Any:
    """
    Detect a single feature by token name.
    
    This is more efficient than detect_all_features when you only need one feature,
    but for multiple features, detect_all_features with caching is more efficient.
    
    Args:
        smiles: SMILES string of the molecule
        feature_token: Name of the feature to detect (e.g., "sp2_chloride_present")
        
    Returns:
        Feature value (bool or int)
        
    Examples:
        >>> detect_feature("c1ccc(Br)cc1", "sp2_bromide_present")
        True
        
        >>> detect_feature("c1ccc(Br)cc1Br", "sp2_halide_site_count")
        2
    """
    # For now, just use detect_all_features and extract the requested feature
    # This leverages caching and is simple to implement
    all_features = detect_all_features(smiles)
    return all_features.get(feature_token, False)


def detect_features_batch(smiles_list: List[str]) -> List[Dict[str, Any]]:
    """
    Batch detection for multiple molecules.
    
    Note: Currently not parallelized, but leverages per-molecule caching.
    Future versions could use multiprocessing for large batches.
    
    Args:
        smiles_list: List of SMILES strings
        
    Returns:
        List of feature dictionaries, one per molecule
        
    Examples:
        >>> detect_features_batch(["c1ccc(Br)cc1", "CCBr", "c1ccccc1"])
        [
            {'sp2_bromide_present': True, 'aryl_halide_present': True, ...},
            {'sp3_bromide_present': True, 'sp2_bromide_present': False, ...},
            {'sp2_bromide_present': False, 'aryl_halide_present': False, ...}
        ]
    """
    return [detect_all_features(smiles) for smiles in smiles_list]


# ============================================================================
# Utility Functions
# ============================================================================

def get_present_features(smiles: str) -> List[str]:
    """
    Get list of feature tokens that are present (True or >0) in a molecule.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        List of feature token names that are present
        
    Examples:
        >>> get_present_features("c1ccc(Br)cc1")
        ['sp2_bromide_present', 'aryl_halide_present', 'ArBr_present', ...]
    """
    all_features = detect_all_features(smiles)
    present = []
    
    for token, value in all_features.items():
        if isinstance(value, bool) and value:
            present.append(token)
        elif isinstance(value, int) and value > 0:
            present.append(token)
    
    return present


def feature_summary(smiles: str) -> str:
    """
    Get a human-readable summary of detected features.
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Multi-line string summarizing the features
        
    Examples:
        >>> print(feature_summary("c1ccc(Br)cc1"))
        Molecule: c1ccc(Br)cc1
        Detected features:
          - sp2_bromide_present
          - aryl_halide_present
          - ArBr_present
          ...
    """
    present = get_present_features(smiles)
    
    lines = [f"Molecule: {smiles}"]
    if present:
        lines.append("Detected features:")
        for token in sorted(present):
            lines.append(f"  - {token}")
    else:
        lines.append("No features detected.")
    
    return "\n".join(lines)


__all__ = [
    "detect_all_features",
    "detect_feature",
    "detect_features_batch",
    "get_present_features",
    "feature_summary",
    "get_feature_spec",
]
