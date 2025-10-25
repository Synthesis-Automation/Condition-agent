"""
Automatic Reactant Type Classifier using SMARTS patterns.

This module provides functions to automatically classify reactants based on their
SMILES structure using the SMARTS patterns defined in reactant_types.json.

Usage:
    from classify_reactant import classify_reactant, classify_batch
    
    # Single molecule
    types = classify_reactant("c1ccccc1Br")
    # Returns: {'category': 'ArX*', 'member_type': 'ArBr', 'name': 'aryl bromide'}
    
    # Batch classification
    results = classify_batch(["c1ccccc1Br", "CCBr", "c1cccnc1Br"])
"""

import json
from typing import Dict, List, Optional, Tuple
from rdkit import Chem


def load_reactant_types() -> Dict:
    """Load reactant types with SMARTS patterns."""
    import os
    
    # Try different paths depending on where the script is run from
    possible_paths = [
        'reactant_types.json',  # Same directory
        'data-processor/other_data/reactant_types.json',  # From root
        '../../data-processor/other_data/reactant_types.json',  # From elsewhere
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as f:
                return json.load(f)
    
    raise FileNotFoundError("Could not find reactant_types.json")


def get_category_matches(smiles: str, reactant_types: Optional[Dict] = None) -> List[str]:
    """
    Get all matching categories (top-level) for a molecule.
    
    Args:
        smiles: SMILES string of the molecule
        reactant_types: Optional pre-loaded reactant types dict
        
    Returns:
        List of matching category names
    """
    if reactant_types is None:
        reactant_types = load_reactant_types()
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    matching_categories = []
    for category, data in reactant_types.items():
        if 'smarts' not in data:
            continue
        
        pattern = Chem.MolFromSmarts(data['smarts'])
        if pattern and mol.HasSubstructMatch(pattern):
            matching_categories.append(category)
    
    return matching_categories


def classify_reactant(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[Dict]:
    """
    Classify a reactant based on SMILES using SMARTS pattern matching.
    
    Uses hierarchical matching:
    1. First checks category-level patterns (broad)
    2. Then finds best specific member within matching categories
    
    Args:
        smiles: SMILES string of the molecule
        reactant_types: Optional pre-loaded reactant types dict
        
    Returns:
        Dict with classification info or None if no match:
        {
            'category': str,  # e.g., 'ArX*'
            'member_type': str,  # e.g., 'ArBr'
            'name': str,  # e.g., 'aryl bromide'
            'group': str,  # e.g., 'Electrophiles'
            'smarts': str,  # The matched SMARTS pattern
            'category_smarts': str  # The category-level SMARTS (if available)
        }
    """
    if reactant_types is None:
        reactant_types = load_reactant_types()
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    
    # Try to match against all patterns
    matches = []
    
    # Define very general categories that should be deprioritized
    # Note: Alkene and Alkyne removed since they now have specific terminal/internal members
    general_categories = {'Alkyl-C-H', 'ArH'}
    
    for category, data in reactant_types.items():
        for member in data['members']:
            if 'smarts' not in member:
                continue
                
            pattern = Chem.MolFromSmarts(member['smarts'])
            if pattern and mol.HasSubstructMatch(pattern):
                is_general = category in general_categories
                matches.append({
                    'category': category,
                    'member_type': member['id'],
                    'name': member['name'],
                    'group': data['group'],
                    'description': data.get('description', ''),
                    'smarts': member['smarts'],
                    'category_smarts': data.get('smarts', ''),  # Include category SMARTS
                    'specificity': len(member['smarts']),  # Longer SMARTS = more specific
                    'is_general': is_general
                })
    
    if not matches:
        return None
    
    # First, try to find non-general matches
    specific_matches = [m for m in matches if not m['is_general']]
    
    # If we have specific matches, use the most specific one
    if specific_matches:
        best_match = max(specific_matches, key=lambda x: x['specificity'])
    else:
        # Otherwise, use the most specific general match
        best_match = max(matches, key=lambda x: x['specificity'])
    
    del best_match['specificity']
    del best_match['is_general']
    return best_match


def classify_batch(smiles_list: List[str], reactant_types: Optional[Dict] = None) -> List[Optional[Dict]]:
    """
    Classify multiple reactants in batch.
    
    Args:
        smiles_list: List of SMILES strings
        reactant_types: Optional pre-loaded reactant types dict
        
    Returns:
        List of classification dicts (or None for unmatched)
    """
    if reactant_types is None:
        reactant_types = load_reactant_types()
    
    return [classify_reactant(smiles, reactant_types) for smiles in smiles_list]


def get_all_matches(smiles: str, reactant_types: Optional[Dict] = None) -> List[Dict]:
    """
    Get all matching reactant types (not just best match).
    
    Useful for molecules that could be classified multiple ways.
    """
    if reactant_types is None:
        reactant_types = load_reactant_types()
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    
    matches = []
    for category, data in reactant_types.items():
        for member in data['members']:
            if 'smarts' not in member:
                continue
                
            pattern = Chem.MolFromSmarts(member['smarts'])
            if pattern and mol.HasSubstructMatch(pattern):
                matches.append({
                    'category': category,
                    'member_type': member['id'],
                    'name': member['name'],
                    'group': data['group'],
                    'smarts': member['smarts']
                })
    
    return matches


def classify_by_category(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[str]:
    """
    Classify and return just the category (e.g., 'ArX*', 'Aliphatic-amine').
    
    Quick function for when you only need the top-level category.
    """
    result = classify_reactant(smiles, reactant_types)
    return result['category'] if result else None


def classify_by_group(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[str]:
    """
    Classify and return just the functional group (e.g., 'Electrophiles', 'Nucleophiles').
    
    Quick function for when you only need the broad group.
    """
    result = classify_reactant(smiles, reactant_types)
    return result['group'] if result else None


# ============================================================================
# Testing and Validation
# ============================================================================

if __name__ == "__main__":
    # Test classification with examples
    test_cases = [
        ("c1ccccc1Br", "ArBr (aryl bromide)"),
        ("CCBr", "Alkyl-Br (alkyl bromide)"),
        ("c1cccnc1Br", "HetArBr or PyridineBr (heteroaryl bromide)"),
        ("C=CCBr", "Allyl-Br (allylic bromide)"),
        ("c1ccccc1CBr", "Bn-Br (benzylic bromide)"),
        ("c1ccccc1B(O)O", "ArB(OH)2 (aryl boronic acid)"),
        ("CCN", "RNH2 (primary aliphatic amine)"),
        ("c1ccccc1N", "ArNH2 (aniline)"),
        ("CCO", "ROH-primary (primary aliphatic alcohol)"),
        ("c1ccccc1O", "ArOH (phenol)"),
        ("CC#N", "R-CN (alkyl nitrile)"),
        ("c1ccccc1C#N", "Ar-CN (aryl nitrile)"),
        ("CC=O", "RCHO? (might be acetaldehyde)"),
        ("C=C", "alkene"),
        ("C#C", "alkyne"),
    ]
    
    print("=" * 80)
    print("REACTANT TYPE CLASSIFIER - Test Suite (with Category-Level Matching)")
    print("=" * 80)
    
    reactant_types = load_reactant_types()
    
    # Test category-level matching
    print("\n--- Testing Category-Level Matching ---")
    for smiles, expected in test_cases[:5]:
        categories = get_category_matches(smiles, reactant_types)
        print(f"\nSMILES: {smiles}")
        print(f"  Expected: {expected}")
        print(f"  Matching Categories: {', '.join(categories) if categories else 'None'}")
    
    print("\n" + "=" * 80)
    print("--- Testing Full Classification (Category + Member) ---")
    print("=" * 80)
    
    reactant_types = load_reactant_types()
    
    for smiles, expected in test_cases:
        result = classify_reactant(smiles, reactant_types)
        all_matches = get_all_matches(smiles, reactant_types)
        
        print(f"\nSMILES: {smiles}")
        print(f"Expected: {expected}")
        
        if result:
            print(f"✅ Best Match: {result['member_type']} ({result['name']})")
            print(f"   Category: {result['category']}")
            print(f"   Group: {result['group']}")
            print(f"   SMARTS: {result['smarts']}")
        else:
            print("❌ No match found")
        
        if len(all_matches) > 1:
            print(f"   ⚠️  {len(all_matches)} total matches:")
            for match in all_matches[:3]:  # Show first 3
                print(f"      - {match['member_type']} ({match['name']})")
    
    print("\n" + "=" * 80)
    print("BATCH CLASSIFICATION TEST")
    print("=" * 80)
    
    batch_smiles = ["c1ccccc1Br", "CCBr", "c1cccnc1Br", "c1ccccc1N", "CCN"]
    results = classify_batch(batch_smiles, reactant_types)
    
    for smiles, result in zip(batch_smiles, results):
        if result:
            print(f"{smiles:20s} → {result['member_type']:20s} ({result['name']})")
        else:
            print(f"{smiles:20s} → No match")
    
    print("\n✅ Classification system ready!")
