"""
Functional Group Detection for Organic Molecules.

Comprehensive SMARTS-based detection of common functional groups with fallback
to text-based pattern matching when RDKit is unavailable.

Usage:
    from chemtools.util.functional_groups import detect_all, get_functional_groups
    
    # Get all functional groups as dict
    groups = detect_all("CC(=O)O")  # {"carboxylic_acid": True, "carbonyl": True, ...}
    
    # Get list of present functional groups
    present = get_functional_groups("CC(=O)O")  # ["carboxylic_acid", "carbonyl"]
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple
from .rdkit_helpers import rdkit_available, parse_smiles


# ============================================================================
# SMARTS Pattern Definitions
# ============================================================================

FUNCTIONAL_GROUP_SMARTS = {
    # Oxygen-containing groups
    "alcohol": "[OX2H]",
    "phenol": "c[OX2H]",
    "ether": "[OD2]([#6])[#6]",
    "carbonyl": "[CX3]=[OX1]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6][CX3](=O)[OX2][#6]",
    "acyl_chloride": "[CX3](=O)[Cl]",
    "anhydride": "[CX3](=O)[OX2][CX3](=O)",
    "peroxide": "[OX2][OX2]",
    
    # Nitrogen-containing groups
    "amine_primary": "[NX3;H2;!$(NC=O)]",
    "amine_secondary": "[NX3;H1;!$(NC=O)]",
    "amine_tertiary": "[NX3;H0;!$(NC=O);!$(N=*)]",
    "aniline": "c[NX3;H2,H1;!$(NC=O)]",
    "amide": "[NX3][CX3](=[OX1])",
    "amide_primary": "[NX3;H2][CX3](=[OX1])",
    "amide_secondary": "[NX3;H1][CX3](=[OX1])",
    "amide_tertiary": "[NX3;H0][CX3](=[OX1])",
    "nitrile": "[NX1]#[CX2]",
    "nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
    "imine": "[NX2]=[CX3]",
    "hydrazine": "[NX3][NX3]",
    "hydroxylamine": "[NX3][OX2H]",
    "isocyanate": "N=C=O",
    "azide": "[NX2]=[NX2]=[NX1]",
    
    # Sulfur-containing groups
    "thiol": "[SX2H]",
    "sulfide": "[#6][SX2][#6]",
    "disulfide": "[SX2][SX2]",
    "sulfoxide": "[$([SX3](=O)[#6]),$([SX3+]([O-])[#6])]",
    "sulfone": "[$([SX4](=O)(=O)([#6])[#6]),$([SX4+2]([O-])([O-])([#6])[#6])]",
    "sulfonyl_chloride": "[SX4](=O)(=O)[Cl]",
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H]",
    "sulfonamide": "[SX4](=O)(=O)[NX3]",
    "sulfate": "[SX4](=O)(=O)([OX2H,OX1-])[OX2][#6]",
    "thioester": "[#6][CX3](=O)[SX2][#6]",
    
    # Halides
    "alkyl_fluoride": "[CX4][F]",
    "alkyl_chloride": "[CX4][Cl]",
    "alkyl_bromide": "[CX4][Br]",
    "alkyl_iodide": "[CX4][I]",
    "aryl_fluoride": "c[F]",
    "aryl_chloride": "c[Cl]",
    "aryl_bromide": "c[Br]",
    "aryl_iodide": "c[I]",
    "acyl_halide": "[CX3](=O)[F,Cl,Br,I]",
    
    # Phosphorus-containing groups
    "phosphine": "[PX3]",
    "phosphine_oxide": "[PX4](=O)",
    "phosphate": "[PX4](=O)([OX2H,OX1-])([OX2H,OX1-])[OX2]",
    
    # Boron-containing groups (for cross-coupling)
    "boronic_acid": "[BX3](O)(O)",
    "boronic_ester": "[BX3](O[#6])(O[#6])",
    "boron_reagent": "[B]",
    
    # Vinyl halides (for cross-coupling)
    "vinyl_halide": "[CX3]=[CX3][F,Cl,Br,I]",
    "vinyl_chloride": "[CX3]=[CX3][Cl]",
    "vinyl_bromide": "[CX3]=[CX3][Br]",
    
    # Other common groups
    "epoxide": "[OX2r3]1[#6][#6]1",
    "aziridine": "[NX3r3]1[#6][#6]1",
    "lactone": "[OX2r][CX3](=O)",
    "lactam": "[NX3r][CX3](=O)",
    "urea": "[NX3][CX3](=[OX1])[NX3]",
    "carbamate": "[NX3][CX3](=[OX1])[OX2]",
    "carbonate": "[OX2][CX3](=[OX1])[OX2]",
    "imide": "[CX3](=O)[NX3][CX3](=O)",
    "enol": "[OX2H][#6]=[#6]",
    "enamine": "[NX3][#6]=[#6]",
    "imine_n_oxide": "[NX3+]([O-])=[CX3]",
    "n_oxide": "[$([#7+][O-]),$([#7v5]=[O])]",
    
    # Aromatic heterocycles (simplified)
    "pyridine": "n1ccccc1",
    "pyrrole": "[nH]1cccc1",
    "furan": "o1cccc1",
    "thiophene": "s1cccc1",
    "imidazole": "n1c[nH]cc1",
    "thiazole": "n1ccsc1",
    "oxazole": "n1ccoc1",
    
    # Functional patterns
    "alkene": "[CX3]=[CX3]",
    "alkyne": "[CX2]#[CX2]",
    "aromatic": "a",
    "benzylic": "[#6]a",
    "allylic": "[#6][CX3]=[CX3]",
    "propargylic": "[#6][CX2]#[CX2]",
    
    # Special leaving groups
    "triflate": "[$(S(=O)(=O)(O[#6])C(F)(F)F)]",
    "tosylate": "[$(S(=O)(=O)(O[#6])c1ccc(C)cc1)]",
    "mesylate": "[$(S(=O)(=O)(O[#6])C)]",
    
    # Protecting groups (common)
    "boc": "[$(C(=O)OC(C)(C)C)]",
    "cbz": "[$(C(=O)OCc1ccccc1)]",
    "fmoc": "[$(C(=O)OCC1c2ccccc2-c2ccccc12)]",
    "silyl_ether": "[$(O[Si])]",
}


# ============================================================================
# Text Pattern Fallbacks (when RDKit unavailable)
# ============================================================================

TEXT_PATTERNS = {
    # Simple substring matches (case-insensitive)
    "alcohol": ["oh", "[oh]"],
    "phenol": ["coh", "c[oh]"],
    "ether": ["oc", "co"],
    "carbonyl": ["c(=o)", "c=o"],
    "aldehyde": ["c(=o)h", "c=oh"],
    "carboxylic_acid": ["c(=o)o", "cooh"],
    "ester": ["oc(=o)", "c(=o)o"],
    "amine_primary": ["nh2", "[nh2]"],
    "amine_secondary": ["nh", "[nh]"],
    "nitrile": ["c#n"],
    "nitro": ["n(=o)=o", "[n+](=o)[o-]"],
    "thiol": ["sh", "[sh]"],
    "sulfone": ["s(=o)(=o)", "so2"],
    "sulfonamide": ["s(=o)(=o)n", "so2n"],
    "alkyl_fluoride": ["cf"],
    "alkyl_chloride": ["ccl"],
    "alkyl_bromide": ["cbr"],
    "alkyl_iodide": ["ci"],
    "aryl_chloride": ["c1ccc(cl)"],
    "aryl_bromide": ["c1ccc(br)"],
    "aryl_iodide": ["c1ccc(i)"],
    "alkene": ["c=c"],
    "alkyne": ["c#c"],
    "triflate": ["otf", "oso2cf3"],
    "boronic_acid": ["b(o)o", "b(oh)oh"],
    "boronic_ester": ["b(o"],
    "vinyl_chloride": ["c=ccl"],
    "vinyl_bromide": ["c=cbr"],
}


# ============================================================================
# Detection Functions
# ============================================================================

def _detect_with_rdkit(smiles: str, patterns: Dict[str, str]) -> Dict[str, bool]:
    """Detect functional groups using RDKit SMARTS matching."""
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return {}
    
    mol = parse_smiles(smiles)
    if mol is None:
        return {}
    
    results = {}
    for name, smarts in patterns.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None:
                results[name] = mol.HasSubstructMatch(pattern)
            else:
                results[name] = False
        except Exception:
            results[name] = False
    
    return results


def _detect_with_text(smiles: str, patterns: Dict[str, List[str]]) -> Dict[str, bool]:
    """Detect functional groups using simple text pattern matching."""
    s_lower = smiles.lower()
    results = {}
    
    for name, pattern_list in patterns.items():
        found = any(p in s_lower for p in pattern_list)
        results[name] = found
    
    return results


def detect_all(smiles: Optional[str]) -> Dict[str, bool]:
    """
    Detect all functional groups in a molecule.
    
    Uses RDKit SMARTS matching when available, falls back to text patterns.
    
    Args:
        smiles: SMILES string to analyze
        
    Returns:
        Dictionary mapping functional group name to presence (True/False)
        
    Examples:
        >>> detect_all("CC(=O)O")
        {'carboxylic_acid': True, 'carbonyl': True, ...}
        
        >>> detect_all("c1ccccc1N")
        {'aniline': True, 'amine_primary': True, 'aromatic': True, ...}
    """
    if not smiles:
        return {name: False for name in FUNCTIONAL_GROUP_SMARTS.keys()}
    
    # Try RDKit first
    if rdkit_available():
        results = _detect_with_rdkit(smiles, FUNCTIONAL_GROUP_SMARTS)
        if results:
            return results
    
    # Fallback to text patterns
    text_results = _detect_with_text(smiles, TEXT_PATTERNS)
    
    # Fill in missing groups as False
    all_groups = set(FUNCTIONAL_GROUP_SMARTS.keys())
    for name in all_groups:
        if name not in text_results:
            text_results[name] = False
    
    return text_results


def get_functional_groups(smiles: Optional[str]) -> List[str]:
    """
    Get list of functional groups present in a molecule.
    
    Args:
        smiles: SMILES string to analyze
        
    Returns:
        List of functional group names found in the molecule
        
    Examples:
        >>> get_functional_groups("CC(=O)O")
        ['carboxylic_acid', 'carbonyl']
        
        >>> get_functional_groups("c1ccc(Br)cc1N")
        ['aniline', 'amine_primary', 'aryl_bromide', 'aromatic']
    """
    all_groups = detect_all(smiles)
    return [name for name, present in all_groups.items() if present]


def has_functional_group(smiles: Optional[str], group_name: str) -> bool:
    """
    Check if a molecule has a specific functional group.
    
    Args:
        smiles: SMILES string to analyze
        group_name: Name of functional group to check (see FUNCTIONAL_GROUP_SMARTS keys)
        
    Returns:
        True if functional group is present, False otherwise
        
    Examples:
        >>> has_functional_group("CC(=O)O", "carboxylic_acid")
        True
        
        >>> has_functional_group("CCO", "aldehyde")
        False
    """
    if not smiles or group_name not in FUNCTIONAL_GROUP_SMARTS:
        return False
    
    # Fast path: check single group
    if rdkit_available():
        result = _detect_with_rdkit(smiles, {group_name: FUNCTIONAL_GROUP_SMARTS[group_name]})
        return result.get(group_name, False)
    
    # Fallback
    if group_name in TEXT_PATTERNS:
        s_lower = smiles.lower()
        return any(p in s_lower for p in TEXT_PATTERNS[group_name])
    
    return False


def count_functional_groups(smiles: Optional[str], group_name: str) -> int:
    """
    Count occurrences of a specific functional group.
    
    Args:
        smiles: SMILES string to analyze
        group_name: Name of functional group to count
        
    Returns:
        Number of times the functional group appears
        
    Examples:
        >>> count_functional_groups("O=C(O)CC(=O)O", "carboxylic_acid")
        2
        
        >>> count_functional_groups("c1ccc(Br)cc1Br", "aryl_bromide")
        2
    """
    if not smiles or group_name not in FUNCTIONAL_GROUP_SMARTS:
        return 0
    
    if not rdkit_available():
        # Fallback: just return 1 if present, 0 if not
        return 1 if has_functional_group(smiles, group_name) else 0
    
    try:
        from rdkit import Chem  # type: ignore
        mol = parse_smiles(smiles)
        if mol is None:
            return 0
        
        pattern = Chem.MolFromSmarts(FUNCTIONAL_GROUP_SMARTS[group_name])
        if pattern is None:
            return 0
        
        matches = mol.GetSubstructMatches(pattern)
        return len(matches)
    except Exception:
        return 0


def get_group_categories(smiles: Optional[str]) -> Dict[str, List[str]]:
    """
    Organize detected functional groups by chemical category.
    
    Args:
        smiles: SMILES string to analyze
        
    Returns:
        Dictionary mapping category to list of functional groups
        
    Examples:
        >>> get_group_categories("CC(=O)Oc1ccccc1")
        {
            'oxygen': ['ester', 'carbonyl', 'phenol'],
            'aromatic': ['aromatic', 'phenol'],
            'halides': []
        }
    """
    groups = get_functional_groups(smiles)
    
    categories = {
        "oxygen": [],
        "nitrogen": [],
        "sulfur": [],
        "phosphorus": [],
        "halides": [],
        "aromatic": [],
        "unsaturated": [],
        "protecting_groups": [],
        "leaving_groups": [],
    }
    
    # Categorize groups
    oxygen_groups = {"alcohol", "phenol", "ether", "carbonyl", "aldehyde", "ketone", 
                     "carboxylic_acid", "ester", "anhydride", "epoxide", "lactone", 
                     "peroxide", "enol", "n_oxide", "imine_n_oxide"}
    nitrogen_groups = {"amine_primary", "amine_secondary", "amine_tertiary", "aniline",
                       "amide", "amide_primary", "amide_secondary", "amide_tertiary",
                       "nitrile", "nitro", "imine", "hydrazine", "hydroxylamine",
                       "isocyanate", "azide", "lactam", "urea", "enamine", "imide"}
    sulfur_groups = {"thiol", "sulfide", "disulfide", "sulfoxide", "sulfone",
                     "sulfonyl_chloride", "sulfonic_acid", "sulfonamide", "sulfate",
                     "thioester"}
    phosphorus_groups = {"phosphine", "phosphine_oxide", "phosphate"}
    halide_groups = {"alkyl_fluoride", "alkyl_chloride", "alkyl_bromide", "alkyl_iodide",
                     "aryl_fluoride", "aryl_chloride", "aryl_bromide", "aryl_iodide",
                     "acyl_halide"}
    aromatic_groups = {"aromatic", "pyridine", "pyrrole", "furan", "thiophene",
                       "imidazole", "thiazole", "oxazole", "benzylic"}
    unsaturated_groups = {"alkene", "alkyne", "allylic", "propargylic"}
    protecting_groups = {"boc", "cbz", "fmoc", "silyl_ether"}
    leaving_groups = {"triflate", "tosylate", "mesylate"}
    
    for group in groups:
        if group in oxygen_groups:
            categories["oxygen"].append(group)
        if group in nitrogen_groups:
            categories["nitrogen"].append(group)
        if group in sulfur_groups:
            categories["sulfur"].append(group)
        if group in phosphorus_groups:
            categories["phosphorus"].append(group)
        if group in halide_groups:
            categories["halides"].append(group)
        if group in aromatic_groups:
            categories["aromatic"].append(group)
        if group in unsaturated_groups:
            categories["unsaturated"].append(group)
        if group in protecting_groups:
            categories["protecting_groups"].append(group)
        if group in leaving_groups:
            categories["leaving_groups"].append(group)
    
    return categories


def summarize_functional_groups(smiles: Optional[str]) -> str:
    """
    Get a human-readable summary of functional groups.
    
    Args:
        smiles: SMILES string to analyze
        
    Returns:
        Formatted string summary
        
    Examples:
        >>> summarize_functional_groups("CC(=O)Oc1ccccc1")
        'Oxygen: ester, carbonyl, phenol\\nAromatic: aromatic, phenol'
    """
    categories = get_group_categories(smiles)
    
    lines = []
    for category, groups in categories.items():
        if groups:
            category_name = category.replace("_", " ").title()
            groups_str = ", ".join(sorted(groups))
            lines.append(f"{category_name}: {groups_str}")
    
    return "\n".join(lines) if lines else "No functional groups detected"


# ============================================================================
# Compatibility Functions (for substrate_analysis.py)
# ============================================================================

def has_free_alcohol(smiles: Optional[str]) -> bool:
    """Check if molecule has free alcohol (not part of carboxylic acid)."""
    if not smiles:
        return False
    
    # Detect both alcohol and carboxylic acid
    has_alcohol = has_functional_group(smiles, "alcohol")
    has_carboxylic = has_functional_group(smiles, "carboxylic_acid")
    has_phenol_group = has_functional_group(smiles, "phenol")
    
    # If has carboxylic acid, check if there are additional alcohols
    if has_carboxylic:
        count_oh = count_functional_groups(smiles, "alcohol")
        count_cooh = count_functional_groups(smiles, "carboxylic_acid")
        # Free alcohol if more OH than COOH (COOH has one OH)
        return count_oh > count_cooh
    
    # Has alcohol or phenol but no carboxylic acid
    return has_alcohol or has_phenol_group


def has_phenol(smiles: Optional[str]) -> bool:
    """Check if molecule has phenol group."""
    return has_functional_group(smiles, "phenol")


def has_sulfonamide(smiles: Optional[str]) -> bool:
    """Check if molecule has sulfonamide group."""
    return has_functional_group(smiles, "sulfonamide")


def has_hydroxylamine(smiles: Optional[str]) -> bool:
    """Check if molecule has hydroxylamine group."""
    return has_functional_group(smiles, "hydroxylamine")
