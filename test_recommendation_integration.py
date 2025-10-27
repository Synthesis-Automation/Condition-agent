#!/usr/bin/env python3
"""
Test script to validate the integration of the new analysis module
into the recommendation pipeline.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from chemtools import chem
import json


def print_separator(char="=", length=80):
    """Print a separator line."""
    print(char * length)


def test_recommendation_with_analysis():
    """Test that the new analysis module is being used in recommendations."""
    
    print_separator()
    print("TESTING NEW ANALYSIS MODULE INTEGRATION")
    print_separator()
    print()
    
    # Test case: Buchwald-Hartwig C-N coupling
    reaction_smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    
    print(f"Input Reaction: {reaction_smiles}")
    print()
    
    print("Running recommendation with new analysis module...")
    print()
    
    try:
        # Get recommendation
        result = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=20,
            search_all_families=False,
            relax={"use_analysis_module": True}  # Explicitly enable new module
        )
        
        # Check detection metadata
        detection = result.get("detection", {})
        
        print_separator("-")
        print("DETECTION INFORMATION")
        print_separator("-")
        print()
        
        print(f"Detected Family: {result.get('family', 'Unknown')}")
        print(f"Auto Family: {detection.get('auto_family', 'N/A')}")
        print(f"Rule Family: {detection.get('rule_family', 'N/A')}")
        print(f"Detection Source: {detection.get('source', 'N/A')}")
        print(f"Analysis Module Used: {detection.get('analysis_module_used', False)}")
        print()
        
        # Check if reactant classification is present
        reactant_classification = detection.get("reactant_classification")
        
        if reactant_classification:
            print_separator("-")
            print("REACTANT CLASSIFICATION (from Analysis Module)")
            print_separator("-")
            print()
            
            print(f"Reaction Type: {reactant_classification.get('reaction_type', 'N/A')}")
            print(f"Confidence: {reactant_classification.get('reaction_confidence', 'N/A')}")
            print(f"Detection Method: {reactant_classification.get('detection_method', 'N/A')}")
            print()
            
            reactants = reactant_classification.get("reactants", [])
            if reactants:
                print(f"Reactants Classified: {len(reactants)}")
                print()
                for i, reactant in enumerate(reactants, 1):
                    print(f"  Reactant {i}:")
                    print(f"    Category: {reactant.get('category', 'N/A')}")
                    print(f"    Type: {reactant.get('member_type', 'N/A')}")
                    print(f"    Name: {reactant.get('name', 'N/A')}")
                    if reactant.get('role'):
                        print(f"    Role: {reactant.get('role')}")
                    print(f"    Expected: {reactant.get('is_expected', False)}")
                    
                    alternatives = reactant.get('alternative_matches', [])
                    if alternatives:
                        print(f"    Alternatives: {len(alternatives)} other functional groups")
                    print()
        else:
            print("WARNING: No reactant classification found - analysis module may not have been used")
            print()
        
        # Show recommendation
        print_separator("-")
        print("RECOMMENDATION")
        print_separator("-")
        print()
        
        recommendation = result.get("recommendation", {})
        print(f"Catalyst: {recommendation.get('core', 'N/A')}")
        print(f"Base: {recommendation.get('base_uid', 'N/A')}")
        print(f"Solvent: {recommendation.get('solvent_uid', 'N/A')}")
        print(f"Temperature: {recommendation.get('T_C', 'N/A')}°C")
        print(f"Time: {recommendation.get('time_h', 'N/A')} h")
        print(f"Confidence: {recommendation.get('confidence', 'N/A')}")
        print()
        
        print_separator()
        print("SUCCESS: TEST COMPLETED")
        print_separator()
        
        return True
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def test_comparison_with_without_analysis():
    """Compare recommendations with and without the analysis module."""
    
    print_separator()
    print("COMPARISON: WITH vs WITHOUT ANALYSIS MODULE")
    print_separator()
    print()
    
    reaction_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    
    print(f"Input Reaction: {reaction_smiles}")
    print()
    
    try:
        # Test WITH analysis module
        print("=" * 40)
        print("WITH ANALYSIS MODULE")
        print("=" * 40)
        result_with = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=20,
            relax={"use_analysis_module": True}
        )
        
        detection_with = result_with.get("detection", {})
        print(f"Detected Family: {result_with.get('family', 'Unknown')}")
        print(f"Analysis Module Used: {detection_with.get('analysis_module_used', False)}")
        if detection_with.get('reactant_classification'):
            reactants = detection_with['reactant_classification'].get('reactants', [])
            print(f"Reactants Classified: {len(reactants)}")
        print()
        
        # Test WITHOUT analysis module
        print("=" * 40)
        print("WITHOUT ANALYSIS MODULE")
        print("=" * 40)
        result_without = chem.recommend.conditions(
            reaction=reaction_smiles,
            k=20,
            relax={"use_analysis_module": False}
        )
        
        detection_without = result_without.get("detection", {})
        print(f"Detected Family: {result_without.get('family', 'Unknown')}")
        print(f"Analysis Module Used: {detection_without.get('analysis_module_used', False)}")
        print()
        
        print_separator()
        print("SUCCESS: COMPARISON COMPLETED")
        print_separator()
        
        return True
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("\n")
    
    # Run tests
    success1 = test_recommendation_with_analysis()
    print("\n\n")
    
    success2 = test_comparison_with_without_analysis()
    print("\n")
    
    # Exit with appropriate code
    if success1 and success2:
        sys.exit(0)
    else:
        sys.exit(1)
