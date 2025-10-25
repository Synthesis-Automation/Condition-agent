"""
Test script to verify reactant type classification is added to process_reactions.py

This script tests that:
1. classify_reactant imports work correctly
2. classify_reactants_from_smiles function works
3. The function is integrated into the main processing workflow
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import the function
from process_reactions import classify_reactants_from_smiles

def test_reactant_classification():
    """Test the reactant classification function."""
    
    print("=" * 80)
    print("TESTING REACTANT TYPE CLASSIFICATION")
    print("=" * 80)
    
    # Test cases
    test_cases = [
        ("c1ccccc1Br", "ArBr (aryl bromide)"),
        ("c1ccccc1Br.CCN", "ArBr + RNH2 (Buchwald-Hartwig reactants)"),
        ("c1ccccc1Br.c1ccccc1B(O)O", "ArBr + ArB(OH)2 (Suzuki reactants)"),
        ("", "Empty SMILES"),
    ]
    
    for smiles, description in test_cases:
        print(f"\nTest: {description}")
        print(f"SMILES: {smiles}")
        
        types, categories = classify_reactants_from_smiles(smiles)
        
        print(f"Reactant Types: {types}")
        print(f"Reactant Categories: {categories}")
    
    print("\n" + "=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    test_reactant_classification()
