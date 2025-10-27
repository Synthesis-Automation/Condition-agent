"""
Test the analysis module with sample reactions from tests/sample_reactions.py

This tests:
1. Reaction type detection (router)
2. Reactant classification (reactants)
3. Context-aware classification (reaction_context)
"""

import sys
import os

# Add project root to path
ROOT = os.path.abspath(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Import sample reactions
sys.path.insert(0, os.path.join(ROOT, 'tests'))
from sample_reactions import SAMPLE_REACTIONS
from chemtools.router import detect_family_from_reaction
from chemtools.analysis.reaction_context import classify_reactants_with_context, get_reactant_summary
from chemtools.analysis.smiles import normalize_reaction


def extract_reaction_smiles(sample: str) -> tuple[str, str]:
    """
    Extract SMILES and description from sample reaction string.
    
    Format: "SMILES (Description)"
    Returns: (smiles, description)
    """
    if "(" in sample and sample.endswith(")"):
        smiles, desc = sample.rsplit("(", 1)
        smiles = smiles.strip()
        desc = desc.rstrip(")").strip()
        return smiles, desc
    return sample.strip(), "Unknown"


def test_reaction_detection():
    """Test reaction type detection on all sample reactions."""
    
    print("=" * 100)
    print("TESTING REACTION TYPE DETECTION")
    print("=" * 100)
    
    # Skip the first entry (it's a placeholder)
    samples = [s for s in SAMPLE_REACTIONS[1:] if s.strip() and not s.startswith("#")]
    
    print(f"\nTesting {len(samples)} sample reactions...")
    print()
    
    results = {
        "total": 0,
        "detected": 0,
        "high_confidence": 0,  # >= 0.8
        "medium_confidence": 0,  # >= 0.5
        "low_confidence": 0,  # < 0.5
        "unknown": 0,
        "by_family": {}
    }
    
    # Test a subset for detailed output
    detailed_count = 10
    
    for i, sample in enumerate(samples, 1):
        smiles, description = extract_reaction_smiles(sample)
        
        # Skip if no valid SMILES
        if not smiles or ">>" not in smiles:
            continue
        
        results["total"] += 1
        
        try:
            # Detect reaction type
            detection = detect_family_from_reaction(smiles)
            family = detection.get("family", "Unknown")
            confidence = detection.get("confidence", 0.0)
            
            # Count by family
            results["by_family"][family] = results["by_family"].get(family, 0) + 1
            
            if family != "Unknown":
                results["detected"] += 1
                
                # Categorize by confidence
                if confidence >= 0.8:
                    results["high_confidence"] += 1
                elif confidence >= 0.5:
                    results["medium_confidence"] += 1
                else:
                    results["low_confidence"] += 1
            else:
                results["unknown"] += 1
            
            # Show detailed output for first N reactions
            if i <= detailed_count:
                print(f"\n{'─' * 100}")
                print(f"Sample {i}: {description}")
                print(f"{'─' * 100}")
                print(f"SMILES: {smiles[:80]}{'...' if len(smiles) > 80 else ''}")
                print(f"\nDetected Family: {family}")
                print(f"Confidence: {confidence:.2f}")
                
                # Show rule hits
                hits = detection.get("hits", {})
                if hits:
                    print(f"\nSMARTS Hits:")
                    hit_list = [k for k, v in hits.items() if v]
                    for hit in sorted(hit_list)[:5]:  # Show first 5
                        print(f"  ✓ {hit}")
                    if len(hit_list) > 5:
                        print(f"  ... and {len(hit_list) - 5} more")
                
                # Show catalyst info
                catalysts = detection.get("catalysts", {})
                if catalysts:
                    metals = catalysts.get("metals", [])
                    if metals:
                        print(f"\nCatalysts Detected: {', '.join(metals)}")
        
        except Exception as e:
            print(f"\n❌ Error in sample {i}: {e}")
            continue
    
    # Print summary statistics
    print(f"\n\n{'=' * 100}")
    print("SUMMARY STATISTICS")
    print(f"{'=' * 100}\n")
    
    print(f"Total Reactions Tested: {results['total']}")
    print(f"Successfully Detected: {results['detected']} ({results['detected']/results['total']*100:.1f}%)")
    print(f"Unknown/Failed: {results['unknown']} ({results['unknown']/results['total']*100:.1f}%)")
    
    print(f"\nConfidence Distribution:")
    print(f"  High (≥0.8):   {results['high_confidence']} ({results['high_confidence']/results['total']*100:.1f}%)")
    print(f"  Medium (≥0.5): {results['medium_confidence']} ({results['medium_confidence']/results['total']*100:.1f}%)")
    print(f"  Low (<0.5):    {results['low_confidence']} ({results['low_confidence']/results['total']*100:.1f}%)")
    
    print(f"\nDetected Reaction Families:")
    sorted_families = sorted(results['by_family'].items(), key=lambda x: x[1], reverse=True)
    for family, count in sorted_families[:15]:  # Top 15
        print(f"  {family:30s}: {count:3d} ({count/results['total']*100:.1f}%)")
    
    if len(sorted_families) > 15:
        print(f"  ... and {len(sorted_families) - 15} more families")
    
    return results


def test_context_aware_classification():
    """Test context-aware reactant classification."""
    
    print(f"\n\n{'=' * 100}")
    print("TESTING CONTEXT-AWARE REACTANT CLASSIFICATION")
    print(f"{'=' * 100}\n")
    
    # Test a few specific reactions
    test_cases = [
        ("Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1", "suzuki_miyaura", "Suzuki - Simple"),
        ("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1", "buchwald_hartwig_c_n", "BH - Simple"),
        ("Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1", "sonogashira", "Sonogashira"),
        ("Brc1ccccc1.C=C>>C(=Cc1ccccc1)", "heck", "Heck"),
    ]
    
    for smiles, expected_type, desc in test_cases:
        print(f"\n{'─' * 100}")
        print(f"Test: {desc}")
        print(f"{'─' * 100}")
        print(f"SMILES: {smiles}")
        print(f"Expected Type: {expected_type}")
        
        try:
            # Test with user-provided reaction type
            result = classify_reactants_with_context(
                smiles,
                reaction_type=expected_type
            )
            
            summary = get_reactant_summary(result)
            
            print(f"\nReaction Type: {summary['reaction_type']}")
            print(f"Detection Method: {summary['detection_method']}")
            print(f"Confidence: {summary['confidence']:.2f}")
            print(f"Number of Reactants: {summary['num_reactants']}")
            
            if summary['num_reactants'] > 0:
                print(f"\nReactants:")
                for r in summary['reactants']:
                    print(f"  {r['position'] + 1}. {r['category']} ({r['member_type']}) - {r['name']}")
                    print(f"     Role: {r['role']}, Expected: {r['is_expected']}")
            else:
                print("\n⚠ Note: No reactants classified (SMARTS pattern issue)")
        
        except Exception as e:
            print(f"\n❌ Error: {e}")


def test_normalization():
    """Test SMILES normalization on sample reactions."""
    
    print(f"\n\n{'=' * 100}")
    print("TESTING SMILES NORMALIZATION")
    print(f"{'=' * 100}\n")
    
    # Test a few reactions
    samples = SAMPLE_REACTIONS[1:6]  # First 5 real reactions
    
    for i, sample in enumerate(samples, 1):
        smiles, description = extract_reaction_smiles(sample)
        
        if not smiles or ">>" not in smiles:
            continue
        
        print(f"\nSample {i}: {description}")
        print(f"Input:  {smiles[:80]}{'...' if len(smiles) > 80 else ''}")
        
        try:
            normalized = normalize_reaction(smiles)
            
            print(f"Output: {normalized['normalized'][:80]}{'...' if len(normalized['normalized']) > 80 else ''}")
            print(f"Reactants: {len(normalized['reactants'])}")
            print(f"Agents: {len(normalized['agents'])}")
            print(f"Products: {len(normalized['products'])}")
            
            if normalized.get('errors'):
                print(f"⚠ Errors: {normalized['errors']}")
        
        except Exception as e:
            print(f"❌ Error: {e}")


def main():
    """Run all tests."""
    
    print("\n" + "█" * 100)
    print("█" + " " * 98 + "█")
    print("█" + " " * 25 + "TESTING ANALYSIS MODULE WITH SAMPLE REACTIONS" + " " * 27 + "█")
    print("█" + " " * 98 + "█")
    print("█" * 100)
    
    try:
        # Test 1: Reaction type detection
        results = test_reaction_detection()
        
        # Test 2: Context-aware classification
        test_context_aware_classification()
        
        # Test 3: Normalization
        test_normalization()
        
        print(f"\n\n{'=' * 100}")
        print("✓ ALL TESTS COMPLETED")
        print(f"{'=' * 100}\n")
        
        return 0
    
    except Exception as e:
        print(f"\n\n{'=' * 100}")
        print(f"❌ FATAL ERROR: {e}")
        print(f"{'=' * 100}\n")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
