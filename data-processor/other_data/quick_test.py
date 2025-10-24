"""
Quick Test Script - Simple examples to validate the recommender
Run this for quick testing without the full test framework.
"""

from simple_condition_recommender import SimpleConditionRecommender
from pathlib import Path
import json

def print_result(result, show_full=False):
    """Pretty print a recommendation result."""
    print(f"\n{'='*80}")
    print(f"Reaction: {result['reaction_type']}")
    print(f"Substrate: {result['substrate']['electrophile']} + {result['substrate']['nucleophile']}")
    print(f"Match Level: {result['match_level'].upper()}")
    print(f"Total Precedents: {result['total_precedents']}")
    print(f"\n{result['metadata']['match_explanation']}")
    
    if result['recommendations']:
        print(f"\nTop {len(result['recommendations'])} Recommendations:")
        print("-" * 80)
        
        for rec in result['recommendations']:
            print(f"\n#{rec['rank']} (Confidence: {rec['confidence_score']:.3f})")
            print(f"  Catalyst: {rec['catalyst']}")
            print(f"  Ligand:   {rec['ligand']}")
            print(f"  Base:     {rec['base']}")
            print(f"  Solvent:  {rec['solvent']}")
            print(f"  Evidence:")
            print(f"    - Success cases: {rec['evidence']['successful_cases']} / {rec['evidence']['total_precedents']} ({rec['evidence']['success_rate']})")
            print(f"    - Avg z-score: {rec['evidence']['avg_zscore']}")
            print(f"    - Z-score range: {rec['evidence']['zscore_range']}")
            print(f"    - Example ELN IDs: {', '.join(rec['evidence']['example_eln_ids'][:3])}")
    else:
        print("\n⚠ No recommendations available")
    
    if show_full:
        print(f"\n{'='*80}")
        print("FULL JSON OUTPUT:")
        print(json.dumps(result, indent=2))


def main():
    # Initialize recommender
    csv_path = Path(__file__).parent / "z-Score Peaks with FG_STANDARDIZED.csv"
    print("Loading recommender...")
    recommender = SimpleConditionRecommender(csv_path, zscore_threshold=1.0)
    
    print("\n" + "="*80)
    print("QUICK TEST EXAMPLES")
    print("="*80)
    
    # Test 1: Very common reaction (should get exact match with many precedents)
    print("\n\n*** TEST 1: Common Buchwald-Hartwig ***")
    result1 = recommender.recommend(
        reaction_type="Buchwald-Hartwig",
        electrophile="ArBr",
        nucleophile="RNH2",
        top_n=3
    )
    print_result(result1)
    
    # Test 2: Common Suzuki (should get exact match)
    print("\n\n*** TEST 2: Common Suzuki-Miyaura ***")
    result2 = recommender.recommend(
        reaction_type="Suzuki-Miyaura",
        electrophile="ArCl",
        nucleophile="ArB(OH)2",
        top_n=3
    )
    print_result(result2)
    
    # Test 3: Less common substrate (might get category match)
    print("\n\n*** TEST 3: Less Common - Amide Coupling ***")
    result3 = recommender.recommend(
        reaction_type="Amide-coupling",
        electrophile="RCO2H",
        nucleophile="RNH2",
        top_n=3
    )
    print_result(result3)
    
    # Test 4: Rare combination (should fall back to category or reaction type)
    print("\n\n*** TEST 4: Rare Combination - C-O Coupling ***")
    result4 = recommender.recommend(
        reaction_type="CO-Coupling",
        electrophile="ArBr",
        nucleophile="ROH-primary",
        top_n=3
    )
    print_result(result4)
    
    # Test 5: Show what happens with completely unknown substrate
    print("\n\n*** TEST 5: Novel Substrate (Fallback Test) ***")
    result5 = recommender.recommend(
        reaction_type="Buchwald-Hartwig",
        electrophile="ArF",  # Less common
        nucleophile="Lactam",  # Less common
        top_n=3
    )
    print_result(result5)
    
    # Summary
    print("\n\n" + "="*80)
    print("SUMMARY OF TESTS")
    print("="*80)
    
    tests = [
        ("Test 1 (BH ArBr+RNH2)", result1),
        ("Test 2 (Suzuki ArCl+ArB)", result2),
        ("Test 3 (Amide RCO2H+RNH2)", result3),
        ("Test 4 (C-O ArBr+ROH)", result4),
        ("Test 5 (BH ArF+Lactam)", result5),
    ]
    
    for name, result in tests:
        status = "✓" if result['recommendations'] else "✗"
        print(f"{status} {name:30s} | Match: {result['match_level']:15s} | Precedents: {result['total_precedents']:4d} | Recs: {len(result['recommendations'])}")
    
    print("\n" + "="*80)
    print("INTERPRETATION:")
    print("="*80)
    print("✓ = Recommendations found")
    print("✗ = No recommendations (reaction type not in database)")
    print("\nMatch levels:")
    print("  - exact: High confidence (exact substrate match)")
    print("  - category: Medium confidence (similar substrates)")
    print("  - reaction_type: Lower confidence (general conditions)")
    print("\nHigh precedent count (>50) = Very reliable recommendations")
    print("Low precedent count (<10) = Less reliable, consider screening")


if __name__ == "__main__":
    main()
