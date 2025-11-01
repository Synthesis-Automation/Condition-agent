"""
Test hybrid bond analysis combining RXNMapper and MCS.

This demonstrates the new hybrid approach that runs both methods
and combines results for higher confidence.
"""

from chemtools import (
    analyze_bond_changes_hybrid,
    analyze_bond_changes,
    analyze_with_mcs,
    rxnmapper_available
)
import json


def test_hybrid_analysis():
    """Test hybrid bond analysis on example reactions."""
    
    print("=" * 80)
    print("Testing Hybrid Bond Analysis (RXNMapper + MCS)")
    print("=" * 80)
    print()
    
    # Check availability
    print(f"RXNMapper available: {rxnmapper_available()}")
    print()
    
    # Test reactions
    test_reactions = [
        {
            "name": "Suzuki-Miyaura Coupling",
            "smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        },
        {
            "name": "Sonogashira Coupling",
            "smiles": "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
        },
        {
            "name": "Amide Formation",
            "smiles": "CC(=O)O.NCC>>CC(=O)NCC"
        },
        {
            "name": "Esterification (Manual Mapped)",
            "smiles": "[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][OH:7]>>[CH3:1][C:2](=[O:3])[O:4][CH2:6][CH3:5]",
            "has_manual": True
        }
    ]
    
    for test in test_reactions:
        print("=" * 80)
        print(f"Reaction: {test['name']}")
        print(f"SMILES: {test['smiles']}")
        if test.get('has_manual'):
            print("  [Manual atom mapping detected]")
        print("-" * 80)
        
        # Run hybrid analysis
        result = analyze_bond_changes_hybrid(
            test['smiles'],
            use_rxnmapper=True,
            use_mcs=True,
            auto_map=True
        )
        
        if result['success']:
            print(f"✓ Analysis successful")
            print(f"  Method: {result['method']}")
            print(f"  Combined confidence: {result['combined_confidence']:.3f}")
            
            # Show agreement details
            if result.get('agreement'):
                agreements = result['agreement']
                print(f"  Agreement:")
                for comparison, agreed in agreements.items():
                    if agreed is not None:
                        status = "✓" if agreed else "✗"
                        print(f"    {status} {comparison.replace('_', ' ')}")
            
            print(f"  Validation: {result.get('validation', 'N/A')}")
            print()
            
            # Manual mapping results
            if result['manual_result']:
                manual = result['manual_result']
                if manual.get('success'):
                    print("  Manual Mapping Results (Ground Truth):")
                    print(f"    Confidence: 1.000 (manual)")
                    print(f"    Broken bonds: {len(manual.get('broken_bonds', []))}")
                    for bond in manual.get('broken_bonds', []):
                        print(f"      - {bond}")
                    print(f"    Formed bonds: {len(manual.get('formed_bonds', []))}")
                    for bond in manual.get('formed_bonds', []):
                        print(f"      - {bond}")
                else:
                    print(f"  Manual mapping failed: {manual.get('error')}")
                print()
            
            # RXNMapper results
            if result['rxnmapper_result']:
                rxn = result['rxnmapper_result']
                if rxn.get('success'):
                    print("  RXNMapper Results:")
                    print(f"    Confidence: {rxn.get('mapping_confidence', 'N/A')}")
                    print(f"    Broken bonds: {len(rxn.get('broken_bonds', []))}")
                    for bond in rxn.get('broken_bonds', []):
                        print(f"      - {bond}")
                    print(f"    Formed bonds: {len(rxn.get('formed_bonds', []))}")
                    for bond in rxn.get('formed_bonds', []):
                        print(f"      - {bond}")
                else:
                    print(f"  RXNMapper failed: {rxn.get('error')}")
                print()
            
            # MCS results
            if result['mcs_result']:
                mcs = result['mcs_result']
                if mcs.get('success'):
                    print("  MCS Results:")
                    print(f"    Confidence: {mcs.get('confidence', 'N/A'):.3f}")
                    print(f"    MCS size: {mcs.get('mcs_size')} atoms")
                    print(f"    MCS coverage: {mcs.get('mcs_coverage', 0):.1%}")
                    print(f"    Estimated broken bonds: {mcs.get('likely_broken_bonds')}")
                    print(f"    Estimated formed bonds: {mcs.get('likely_formed_bonds')}")
                    print(f"    Interpretation: {mcs.get('interpretation')}")
                    if mcs.get('warning'):
                        print(f"    Warning: {mcs['warning']}")
                else:
                    print(f"  MCS failed: {mcs.get('error')}")
                print()
            
            # Recommended result
            if result.get('recommended_result'):
                print(f"  ✓ Recommended: Use {result['method']} result")
                
        else:
            print(f"✗ Analysis failed: {result.get('error')}")
        
        print()


def test_individual_methods():
    """Compare individual methods side-by-side."""
    
    print("\n" + "=" * 80)
    print("Individual Method Comparison")
    print("=" * 80)
    print()
    
    test_smiles = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    
    print(f"Reaction: {test_smiles}")
    print()
    
    # RXNMapper only
    print("-" * 80)
    print("RXNMapper Analysis:")
    print("-" * 80)
    rxn_result = analyze_bond_changes(test_smiles, auto_map=True)
    print(json.dumps({
        'success': rxn_result.get('success'),
        'confidence': rxn_result.get('mapping_confidence'),
        'broken_bonds': rxn_result.get('broken_bonds', []),
        'formed_bonds': rxn_result.get('formed_bonds', [])
    }, indent=2))
    print()
    
    # MCS only
    print("-" * 80)
    print("MCS Analysis:")
    print("-" * 80)
    mcs_result = analyze_with_mcs(test_smiles)
    print(json.dumps({
        'success': mcs_result.get('success'),
        'confidence': mcs_result.get('confidence'),
        'mcs_size': mcs_result.get('mcs_size'),
        'likely_broken_bonds': mcs_result.get('likely_broken_bonds'),
        'likely_formed_bonds': mcs_result.get('likely_formed_bonds'),
        'interpretation': mcs_result.get('interpretation')
    }, indent=2))
    print()
    
    # Hybrid
    print("-" * 80)
    print("Hybrid Analysis:")
    print("-" * 80)
    hybrid_result = analyze_bond_changes_hybrid(test_smiles)
    print(json.dumps({
        'success': hybrid_result.get('success'),
        'method': hybrid_result.get('method'),
        'combined_confidence': hybrid_result.get('combined_confidence'),
        'agreement': hybrid_result.get('agreement'),
        'validation': hybrid_result.get('validation')
    }, indent=2))
    print()


if __name__ == "__main__":
    # Run hybrid tests
    test_hybrid_analysis()
    
    # Compare individual methods
    test_individual_methods()
    
    print("\n" + "=" * 80)
    print("Test Complete!")
    print("=" * 80)
