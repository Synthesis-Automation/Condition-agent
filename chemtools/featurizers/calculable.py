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
from ..util.smarts_cache import compile_smarts as _compile_smarts
from ..util.boolean_expr import evaluate as _eval_bool_expr

# Path to the feature specification JSON (centralized in taxonomy/data)
_SPEC_PATH = Path(__file__).resolve().parents[1] / "taxonomy" / "data" / "calculable_features.json"


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
        Dictionary containing version, schema_notes, features, and optional
        derived_shortcuts.
    """
    return _load_feature_spec()


@lru_cache(maxsize=1)
def _get_derived_rules() -> List[tuple[str, str]]:
    """
    Collect all derived boolean rules from the feature spec.

    The spec currently supports two encodings:
      - `features[].derive` or `features[].derived`
      - `derived_shortcuts[].derive`

    This helper normalizes both into a single ordered list.
    """
    spec = _load_feature_spec()
    rules: List[tuple[str, str]] = []

    for entry in spec.get("features", []):
        token = entry.get("token")
        if not token:
            continue
        expr = entry.get("derive") or entry.get("derived")
        if isinstance(expr, str) and expr.strip():
            rules.append((str(token), expr.strip()))

    for entry in spec.get("derived_shortcuts", []):
        token = entry.get("token")
        expr = entry.get("derive")
        if not token:
            continue
        if isinstance(expr, str) and expr.strip():
            rules.append((str(token), expr.strip()))

    return rules


def _apply_derived_rules(result: Dict[str, Any], rules: List[tuple[str, str]]) -> None:
    """
    Apply derived boolean rules in a small fixed-point loop.

    Some derived tokens depend on other derived tokens. Instead of relying on
    ordering in JSON, iterate until values stabilize (or a small cap).
    """
    max_passes = 8
    for _ in range(max_passes):
        changed = False
        for token, expr in rules:
            value = _eval_bool_expr(expr, result)
            if result.get(token) is not value:
                result[token] = value
                changed = True
        if not changed:
            return


# ============================================================================
# SMARTS Pattern Compilation and Caching (now uses centralized cache)
# ============================================================================
# Note: _compile_smarts is now imported from util.smarts_cache for global caching


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


def _detect_ortho_substitution(mol) -> bool:
    """
    Detect ortho-disubstituted benzene rings.
    
    This looks for benzene rings with two substituents in ortho (1,2) positions,
    which can cause steric hindrance in reactions.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        True if ortho-disubstituted pattern is found
    """
    if mol is None or not rdkit_available():
        return False
    
    try:
        from rdkit import Chem
        
        # SMARTS for ortho-disubstituted benzene
        # Pattern: Two adjacent aromatic carbons in a 6-membered aromatic ring,
        # both with non-ring substituents
        ortho_patterns = [
            # Generic ortho-disubstituted pattern
            "c1ccccc1(*)(*)",  # Two adjacent positions with substituents
            # More specific: both positions have non-H substituents
            "c1c([!H])c([!H])ccc1",
        ]
        
        for smarts in ortho_patterns:
            pattern = _compile_smarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                return True
        
        return False
    except Exception:
        return False


def _detect_molecular_weight(mol) -> float:
    """
    Calculate molecular weight.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Molecular weight in g/mol
    """
    if mol is None or not rdkit_available():
        return 0.0
    
    try:
        from rdkit.Chem import Descriptors
        return Descriptors.MolWt(mol)
    except Exception:
        return 0.0


def _detect_fused_ring_system(mol) -> bool:
    """
    Detect fused ring systems (e.g., naphthalene, quinoline).
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        True if molecule has fused rings
    """
    if mol is None or not rdkit_available():
        return False
    
    try:
        from rdkit import Chem
        
        # Get ring info
        ri = mol.GetRingInfo()
        rings = ri.AtomRings()
        
        if len(rings) < 2:
            return False
        
        # Check if any two rings share at least one edge (two atoms)
        for i, ring1 in enumerate(rings):
            for ring2 in rings[i+1:]:
                shared_atoms = set(ring1) & set(ring2)
                if len(shared_atoms) >= 2:
                    return True
        
        return False
    except Exception:
        return False


def _detect_chiral_centers(mol) -> tuple[bool, int]:
    """
    Detect chiral centers in molecule.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Tuple of (has_chiral_center, chiral_center_count)
    """
    if mol is None or not rdkit_available():
        return (False, 0)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import FindMolChiralCenters
        
        # Find chiral centers
        chiral_centers = FindMolChiralCenters(mol, includeUnassigned=True)
        count = len(chiral_centers)
        
        return (count > 0, count)
    except Exception:
        return (False, 0)


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
    
    # Ortho-substitution detection
    if token == 'ortho_substitution_present':
        return _detect_ortho_substitution(mol)
    
    # Molecular weight features
    if token in ('low_molecular_weight', 'high_molecular_weight'):
        mw = _detect_molecular_weight(mol)
        if token == 'low_molecular_weight':
            return mw < 200
        else:  # high_molecular_weight
            return mw > 500
    
    # Fused ring system detection
    if token == 'fused_ring_system':
        return _detect_fused_ring_system(mol)
    
    # Chiral center detection
    if token == 'chiral_center_present':
        has_chiral, _ = _detect_chiral_centers(mol)
        return has_chiral
    
    if token == 'chiral_center_count':
        _, count = _detect_chiral_centers(mol)
        return count
    
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
        - Comparisons: "halogen_count >= 2", "ring_count > 0"
        - Combinations: "feature1 AND NOT feature2"
    
    Args:
        derive_expr: Expression like "alkyne_present AND NOT terminal_alkyne_present"
        base_features: Dictionary of already-computed features
        
    Returns:
        Boolean result of the expression
    """
    return _eval_bool_expr(derive_expr, base_features)



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
        for token, _expr in _get_derived_rules():
            result[token] = False
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
        for token, _expr in _get_derived_rules():
            result[token] = False
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
        
        # SMARTS count (always returns int)
        elif "smarts_count" in detect:
            smarts_pattern = detect["smarts_count"]
            result[token] = _count_substructure_matches(mol, [smarts_pattern])
        
        # Heuristic-based detection
        elif "heuristic" in detect:
            heuristic_desc = detect["heuristic"]
            result[token] = _detect_heuristic_features(mol, heuristic_desc, token)
        
        else:
            # Unknown detection method - return safe default
            result[token] = False if ftype == "bool" else 0

    # Process all derived features (both features[].derive and derived_shortcuts[]).
    _apply_derived_rules(result, _get_derived_rules())

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


# ============================================================================
# Reactant Type Utility Functions
# ============================================================================

def get_reactant_type_features(smiles: str) -> Dict[str, Any]:
    """
    Extract only reactant type features from a molecule.
    
    Returns dict with both member-level and category-level reactant types:
    {
        'ArBr_reactant': True,
        'ArX_reactant': True,
        'member_types': ['ArBr'],
        'categories': ['ArX*'],
        ...
    }
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Dictionary with reactant type features and metadata
        
    Examples:
        >>> get_reactant_type_features("c1ccc(Br)cc1")
        {'ArBr_reactant': True, 'ArX_reactant': True, 'member_types': ['ArBr'], ...}
    """
    all_features = detect_all_features(smiles)
    
    # Extract reactant features (those with reactant_metadata)
    spec = _load_feature_spec()
    reactant_features = {}
    member_types = []
    categories = []
    
    for feature in spec.get("features", []):
        token = feature.get("token")
        reactant_meta = feature.get("reactant_metadata")
        
        if reactant_meta and all_features.get(token):
            reactant_features[token] = all_features[token]
            
            # Extract member type from metadata
            if "reactant_member" in reactant_meta:
                member_types.append(reactant_meta["reactant_member"])
            if "reactant_category" in reactant_meta:
                cat = reactant_meta["reactant_category"]
                if cat not in categories:
                    categories.append(cat)
    
    # Add derived category features
    for derived in spec.get("derived_shortcuts", []):
        token = derived.get("token")
        reactant_meta = derived.get("reactant_metadata")
        
        if reactant_meta and reactant_features.get(token):
            if "reactant_category" in reactant_meta:
                cat = reactant_meta["reactant_category"]
                if cat not in categories:
                    categories.append(cat)
    
    reactant_features["member_types"] = member_types
    reactant_features["categories"] = categories
    
    return reactant_features


def classify_reactant_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    """
    Backward-compatible wrapper for reactant classification.
    
    Returns a dict similar to the legacy ReactantMatch structure:
    {
        'category': 'ArX*',
        'member_type': 'ArBr',
        'name': 'aryl bromide',
        'smarts': 'c[Br]',
        'coupling_role': 'electrophile',
        ...
    }
    
    Args:
        smiles: SMILES string of the molecule
        
    Returns:
        Dictionary with reactant match info, or None if no match
        
    Examples:
        >>> classify_reactant_smiles("c1ccc(Br)cc1")
        {'category': 'ArX*', 'member_type': 'ArBr', 'name': 'aryl bromide', ...}
    """
    reactant_features = get_reactant_type_features(smiles)
    
    # If no reactant types detected, return None
    if not reactant_features.get("member_types"):
        return None
    
    # Find the most specific member match
    spec = _load_feature_spec()
    
    # General categories that should be deprioritized
    GENERAL_REACTANT_CATEGORIES = {"Alkyl-C-H", "ArH"}
    
    # Priority: prefer non-general categories, then longer SMARTS patterns
    candidates = []
    
    for feature in spec.get("features", []):
        token = feature.get("token")
        reactant_meta = feature.get("reactant_metadata")
        
        if reactant_meta and reactant_features.get(token):
            smarts = feature.get("detect", {}).get("smarts_any", [""])[0]
            category = reactant_meta.get("reactant_category", "")
            member_type = reactant_meta.get("reactant_member", "")
            is_general = category in GENERAL_REACTANT_CATEGORIES
            
            candidates.append({
                "category": category,
                "member_type": member_type,
                "name": reactant_meta.get("compound_type", ""),
                "category_name": category,
                "smarts": smarts,
                "coupling_role": reactant_meta.get("coupling_role", ""),
                "description": "",
                "specificity": len(smarts),
                "is_general": is_general,
                "group": "",
                "category_smarts": "",
            })
    
    if not candidates:
        return None
    
    # Sort by (general? -> False first, specificity descending, member type for determinism)
    candidates.sort(key=lambda m: (m["is_general"], -m["specificity"], m["member_type"]))
    return candidates[0]


def get_reactant_categories(smiles: str) -> List[str]:
    """
    Get list of reactant categories that match a molecule.
    
    Args:
        smiles: SMILES string
        
    Returns:
        List of category IDs (e.g., ['ArX*', 'HetAr-X'])
        
    Examples:
        >>> get_reactant_categories("c1ccc(Br)cc1")
        ['ArX*']
    """
    features = get_reactant_type_features(smiles)
    return features.get("categories", [])


def get_reactant_members(smiles: str) -> List[str]:
    """
    Get list of reactant member types that match a molecule.
    
    Args:
        smiles: SMILES string
        
    Returns:
        List of member IDs (e.g., ['ArBr', 'ArCl'])
        
    Examples:
        >>> get_reactant_members("c1ccc(Br)cc1")
        ['ArBr']
    """
    features = get_reactant_type_features(smiles)
    return features.get("member_types", [])


__all__ = [
    "detect_all_features",
    "detect_feature",
    "detect_features_batch",
    "get_present_features",
    "feature_summary",
    "get_feature_spec",
    # Reactant type functions
    "get_reactant_type_features",
    "classify_reactant_smiles",
    "get_reactant_categories",
    "get_reactant_members",
]
