"""
Minimal test case to debug SMARTS pattern matching issue.

This script tests the most basic case: simple bromobenzene should match ArBr pattern.
"""

import sys
from chemtools.analysis.reactants import iter_reactant_matches, classify_reactant_smiles
from chemtools.analysis._registry import get_registry


def test_simple_molecules():
    """Test SMARTS matching on simple, well-known molecules."""
    
    print("=" * 80)
    print("MINIMAL SMARTS MATCHING DEBUG TEST")
    print("=" * 80)
    print()
    
    # Test cases: (SMILES, Expected Category)
    test_cases = [
        ("Brc1ccccc1", "ArBr / aryl_halide"),
        ("Clc1ccccc1", "ArCl / aryl_halide"),
        ("Ic1ccccc1", "ArI / aryl_halide"),
        ("Nc1ccccc1", "ArNH2 / amine"),
        ("c1ccc(B(O)O)cc1", "ArB(OH)2 / boron"),
        ("c1ccc(C#N)cc1", "ArCN / nitrile"),
        ("CCBr", "Alkyl-Br / alkyl_halide"),
        ("CC(C)N", "RNH2 / amine"),
    ]
    
    print(f"Testing {len(test_cases)} simple molecules...\n")
    
    for smiles, expected in test_cases:
        print("-" * 80)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}")
        print()
        
        # Test with iter_reactant_matches
        matches = list(iter_reactant_matches(smiles))
        print(f"Matches found: {len(matches)}")
        
        if matches:
            print("SUCCESS - Matches:")
            for i, match in enumerate(matches, 1):
                print(f"  {i}. Category: {match.category}")
                print(f"     Member Type: {match.member_type}")
                print(f"     Name: {match.name}")
                print(f"     Description: {match.description}")
        else:
            print("FAILED - No matches returned")
            
            # Try to debug why
            print("\nDebugging:")
            
            # Check if registry is loaded
            registry = get_registry()
            if registry:
                print(f"  ✓ Registry loaded")
                print(f"  ✓ Reactant types: {len(registry.reactant_types)}")
                
                # Try to manually check a few patterns
                if hasattr(registry, 'reactant_types'):
                    sample_types = list(registry.reactant_types.keys())[:5]
                    print(f"  ✓ Sample types: {sample_types}")
            else:
                print(f"  ✗ Registry not loaded!")
        
        print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    total_tests = len(test_cases)
    successful = sum(1 for smiles, _ in test_cases if list(iter_reactant_matches(smiles)))
    failed = total_tests - successful
    
    print(f"Total Tests: {total_tests}")
    print(f"Successful: {successful} ({successful/total_tests*100:.1f}%)")
    print(f"Failed: {failed} ({failed/total_tests*100:.1f}%)")
    print()
    
    if failed == total_tests:
        print("⚠ ALL TESTS FAILED - SMARTS matching is completely broken")
        return 1
    elif failed > 0:
        print("⚠ PARTIAL FAILURE - Some patterns not matching")
        return 1
    else:
        print("✓ ALL TESTS PASSED - SMARTS matching working correctly")
        return 0


def test_registry_inspection():
    """Inspect the registry to understand the data structure."""
    
    print("\n")
    print("=" * 80)
    print("REGISTRY INSPECTION")
    print("=" * 80)
    print()
    
    registry = get_registry()
    
    if not registry:
        print("✗ Registry is None - cannot load taxonomy data")
        return
    
    print(f"✓ Registry loaded")
    print()
    
    # Check reactant types
    if hasattr(registry, 'reactant_types'):
        print(f"Reactant Types Count: {len(registry.reactant_types)}")
        
        # Show first 10 reactant types
        print("\nFirst 10 Reactant Types:")
        for i, (key, value) in enumerate(list(registry.reactant_types.items())[:10], 1):
            print(f"  {i}. {key}: {type(value)}")
            if isinstance(value, dict):
                print(f"     Keys: {list(value.keys())[:5]}")
            
        # Look for specific patterns we expect
        print("\nLooking for expected patterns:")
        expected_patterns = ["ArBr", "ArX", "aryl_halide", "ArNH2", "amine"]
        for pattern in expected_patterns:
            if pattern in registry.reactant_types:
                print(f"  ✓ Found: {pattern}")
                rt = registry.reactant_types[pattern]
                if isinstance(rt, dict) and 'smarts' in rt:
                    print(f"    SMARTS: {rt['smarts'][:100]}")
            else:
                print(f"  ✗ Missing: {pattern}")
    else:
        print("✗ No reactant_types attribute")
    
    # Check reaction types
    if hasattr(registry, 'reaction_types'):
        print(f"\nReaction Types Count: {len(registry.reaction_types)}")
        sample_reactions = list(registry.reaction_types.keys())[:10]
        print(f"Sample reactions: {sample_reactions}")


def test_direct_smarts():
    """Test SMARTS pattern matching directly with RDKit."""
    
    print("\n")
    print("=" * 80)
    print("DIRECT RDKIT SMARTS TEST")
    print("=" * 80)
    print()
    
    try:
        from rdkit import Chem
        
        # Simple test: bromobenzene should match aryl bromide SMARTS
        mol = Chem.MolFromSmiles("Brc1ccccc1")
        
        if mol is None:
            print("✗ Failed to create RDKit mol from 'Brc1ccccc1'")
            return
        
        print("✓ Created RDKit mol from 'Brc1ccccc1'")
        print(f"  Atoms: {mol.GetNumAtoms()}")
        print(f"  Bonds: {mol.GetNumBonds()}")
        print()
        
        # Test common SMARTS patterns
        patterns = [
            ("[Br]", "Bromine atom"),
            ("c", "Aromatic carbon"),
            ("[cH]", "Aromatic CH"),
            ("c1ccccc1", "Benzene ring"),
            ("[Br;$([Br]c)]", "Aryl bromide (simple)"),
            ("[$([Br;$([Br]c)])]", "Aryl bromide (complex)"),
        ]
        
        print("Testing SMARTS patterns:")
        for smarts, description in patterns:
            try:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.GetSubstructMatches(pattern)
                    print(f"  {smarts:30s} ({description:30s}): {len(matches)} matches")
                else:
                    print(f"  {smarts:30s} ({description:30s}): Invalid SMARTS")
            except Exception as e:
                print(f"  {smarts:30s} ({description:30s}): Error - {e}")
        
    except ImportError:
        print("✗ RDKit not available - cannot test SMARTS directly")


if __name__ == "__main__":
    print("\n")
    
    # Run tests
    exit_code = test_simple_molecules()
    test_registry_inspection()
    test_direct_smarts()
    
    print("\n")
    print("=" * 80)
    print("DEBUG TEST COMPLETED")
    print("=" * 80)
    
    sys.exit(exit_code)
