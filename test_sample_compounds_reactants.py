"""
Test the unified reactant type system with sample_compounds.py
============================================================

This script tests the new unified feature system (chemtools.featurizers.calculable)
against all compounds in sample_compounds.py and validates that reactant type
detection works correctly.
"""

import sys
from pathlib import Path

# Add tests directory to path
tests_dir = Path(__file__).parent / "tests"
sys.path.insert(0, str(tests_dir))

from sample_compounds import ALL_SAMPLE_COMPOUNDS, ARYL_HALIDES, PHASE_3_4_COMPOUNDS
from chemtools.featurizers.calculable import (
    classify_reactant_smiles,
    get_reactant_type_features,
    get_reactant_categories,
    get_reactant_members,
    detect_all_features,
)
from chemtools.analysis.reactants import (
    classify_reactant_smiles as analysis_classify,
    iter_reactant_matches,
)
from typing import Dict, Any, List
import sys


def test_compound(compound: Dict[str, Any]) -> Dict[str, Any]:
    """Test a single compound and return results."""
    smiles = compound["smiles"]
    name = compound["name"]
    role = compound["role"]
    
    # Test with unified feature system
    result = classify_reactant_smiles(smiles)
    reactant_features = get_reactant_type_features(smiles)
    categories = get_reactant_categories(smiles)
    members = get_reactant_members(smiles)
    
    # Test with analysis module (should use same system now)
    analysis_result = analysis_classify(smiles)
    
    # Get all features for comprehensive check
    all_features = detect_all_features(smiles)
    
    return {
        "name": name,
        "smiles": smiles,
        "role": role,
        "unified_result": result,
        "analysis_result": analysis_result,
        "categories": categories,
        "members": members,
        "member_count": len(members),
        "category_count": len(categories),
        "has_reactant_type": result is not None,
        "unified_and_analysis_match": (
            result is not None and analysis_result is not None and
            result.get("member_type") == analysis_result.member_type
        ),
    }


def print_summary_table(results: List[Dict[str, Any]]):
    """Print a summary table of results."""
    print("\n" + "=" * 120)
    print("REACTANT TYPE DETECTION SUMMARY")
    print("=" * 120)
    print(f"{'Name':<40} {'Role':<15} {'Detected?':<12} {'Member':<20} {'Category':<20}")
    print("-" * 120)
    
    for r in results:
        name = r["name"][:38]
        role = r["role"][:13]
        detected = "✓" if r["has_reactant_type"] else "✗"
        member = r["unified_result"]["member_type"][:18] if r["unified_result"] else "N/A"
        category = r["unified_result"]["category"][:18] if r["unified_result"] else "N/A"
        
        print(f"{name:<40} {role:<15} {detected:<12} {member:<20} {category:<20}")


def print_detailed_results(results: List[Dict[str, Any]], limit: int = 10):
    """Print detailed results for first N compounds."""
    print("\n" + "=" * 120)
    print(f"DETAILED RESULTS (First {limit} compounds)")
    print("=" * 120)
    
    for i, r in enumerate(results[:limit]):
        print(f"\n[{i+1}] {r['name']}")
        print(f"    SMILES: {r['smiles']}")
        print(f"    Role: {r['role']}")
        
        if r["unified_result"]:
            result = r["unified_result"]
            print(f"    ✓ Detected: {result['member_type']} ({result['category']})")
            print(f"      Name: {result['name']}")
            print(f"      Coupling role: {result.get('coupling_role', 'N/A')}")
            print(f"      Specificity: {result.get('specificity', 0)}")
        else:
            print(f"    ✗ No reactant type detected")
        
        print(f"    Members detected: {len(r['members'])}")
        if r["members"]:
            print(f"      {', '.join(r['members'][:5])}" + ("..." if len(r['members']) > 5 else ""))
        
        print(f"    Categories detected: {len(r['categories'])}")
        if r["categories"]:
            print(f"      {', '.join(r['categories'][:3])}" + ("..." if len(r['categories']) > 3 else ""))
        
        # Check system consistency
        if r["unified_and_analysis_match"]:
            print(f"    ✓ Unified & Analysis systems match")
        elif r["unified_result"] and r["analysis_result"]:
            print(f"    ⚠ System mismatch!")
            print(f"      Unified: {r['unified_result']['member_type']}")
            print(f"      Analysis: {r['analysis_result'].member_type}")


def analyze_by_role(results: List[Dict[str, Any]]):
    """Analyze detection rates by role."""
    print("\n" + "=" * 120)
    print("DETECTION RATES BY ROLE")
    print("=" * 120)
    
    roles = {}
    for r in results:
        role = r["role"]
        if role not in roles:
            roles[role] = {"total": 0, "detected": 0}
        roles[role]["total"] += 1
        if r["has_reactant_type"]:
            roles[role]["detected"] += 1
    
    print(f"{'Role':<20} {'Total':<10} {'Detected':<10} {'Rate':<10}")
    print("-" * 60)
    
    for role, stats in sorted(roles.items()):
        rate = f"{stats['detected']/stats['total']*100:.1f}%"
        print(f"{role:<20} {stats['total']:<10} {stats['detected']:<10} {rate:<10}")
    
    total = sum(s["total"] for s in roles.values())
    detected = sum(s["detected"] for s in roles.values())
    overall_rate = f"{detected/total*100:.1f}%"
    print("-" * 60)
    print(f"{'OVERALL':<20} {total:<10} {detected:<10} {overall_rate:<10}")


def find_issues(results: List[Dict[str, Any]]):
    """Find and report issues."""
    print("\n" + "=" * 120)
    print("ISSUE DETECTION")
    print("=" * 120)
    
    # Find compounds that should be detected but aren't
    expected_roles = {"electrophile", "nucleophile", "bifunctional"}
    undetected = [r for r in results if r["role"] in expected_roles and not r["has_reactant_type"]]
    
    if undetected:
        print(f"\n⚠ {len(undetected)} compounds not detected (expected to have reactant type):")
        for r in undetected[:10]:
            print(f"  - {r['name']} ({r['role']}): {r['smiles']}")
        if len(undetected) > 10:
            print(f"  ... and {len(undetected) - 10} more")
    else:
        print("\n✓ All expected compounds detected successfully!")
    
    # Find system mismatches
    mismatches = [r for r in results if not r["unified_and_analysis_match"] and r["unified_result"] and r["analysis_result"]]
    
    if mismatches:
        print(f"\n⚠ {len(mismatches)} system mismatches (unified vs analysis):")
        for r in mismatches[:5]:
            print(f"  - {r['name']}")
            print(f"    Unified: {r['unified_result']['member_type']}")
            print(f"    Analysis: {r['analysis_result'].member_type}")
    else:
        print("\n✓ All systems consistent (unified = analysis)!")


def test_specific_categories():
    """Test specific reactant categories mentioned in sample_compounds."""
    print("\n" + "=" * 120)
    print("SPECIFIC CATEGORY TESTS")
    print("=" * 120)
    
    test_cases = [
        ("Brc1ccccc1", "ArBr", "ArX*", "Bromobenzene"),
        ("c1ccc(B(O)O)cc1", "ArB(OH)2", "ArB*", "Phenylboronic acid"),
        ("c1ccc(N)cc1", "ArNH2", "ArNH2/Ar2NH", "Aniline"),
        ("CCN", "RNH2", "RNH2/R2NH", "Ethylamine"),
        ("C#Cc1ccccc1", "terminal-alkyne", "Alkyne", "Phenylacetylene"),
        ("CCBr", "Alkyl-Br", "Alkyl-X", "Ethyl bromide"),
    ]
    
    print(f"{'Name':<25} {'SMILES':<20} {'Expected Member':<20} {'Detected':<20} {'Match':<10}")
    print("-" * 120)
    
    for smiles, expected_member, expected_cat, name in test_cases:
        result = classify_reactant_smiles(smiles)
        if result:
            detected_member = result["member_type"]
            match = "✓" if detected_member == expected_member else "✗"
        else:
            detected_member = "N/A"
            match = "✗"
        
        print(f"{name:<25} {smiles:<20} {expected_member:<20} {detected_member:<20} {match:<10}")


def main():
    """Main test function."""
    print("=" * 120)
    print("TESTING UNIFIED REACTANT TYPE SYSTEM WITH SAMPLE COMPOUNDS")
    print("=" * 120)
    print(f"\nTotal compounds to test: {len(ALL_SAMPLE_COMPOUNDS)}")
    print(f"  - ARYL_HALIDES: {len(ARYL_HALIDES)}")
    print(f"  - PHASE_3_4_COMPOUNDS: {len(PHASE_3_4_COMPOUNDS)}")
    
    # Run tests on all compounds
    print("\nRunning detection tests...")
    results = []
    for compound in ALL_SAMPLE_COMPOUNDS:
        try:
            result = test_compound(compound)
            results.append(result)
        except Exception as e:
            print(f"ERROR testing {compound['name']}: {e}")
            results.append({
                "name": compound["name"],
                "smiles": compound["smiles"],
                "role": compound["role"],
                "unified_result": None,
                "analysis_result": None,
                "categories": [],
                "members": [],
                "member_count": 0,
                "category_count": 0,
                "has_reactant_type": False,
                "unified_and_analysis_match": False,
                "error": str(e),
            })
    
    # Print results
    print_summary_table(results)
    analyze_by_role(results)
    find_issues(results)
    test_specific_categories()
    print_detailed_results(results, limit=15)
    
    # Final summary
    total = len(results)
    detected = sum(1 for r in results if r["has_reactant_type"])
    consistent = sum(1 for r in results if r["unified_and_analysis_match"])
    
    print("\n" + "=" * 120)
    print("FINAL SUMMARY")
    print("=" * 120)
    print(f"Total compounds tested: {total}")
    print(f"Reactant types detected: {detected} ({detected/total*100:.1f}%)")
    print(f"System consistency (unified = analysis): {consistent} ({consistent/total*100:.1f}%)")
    
    # Success criteria
    detection_rate = detected / total
    consistency_rate = consistent / total
    
    if detection_rate > 0.70 and consistency_rate > 0.95:
        print("\n✅ ALL TESTS PASSED!")
        print(f"   - Detection rate: {detection_rate*100:.1f}% (target: >70%)")
        print(f"   - Consistency rate: {consistency_rate*100:.1f}% (target: >95%)")
        return 0
    else:
        print("\n⚠️ SOME TESTS NEED ATTENTION")
        if detection_rate <= 0.70:
            print(f"   - Detection rate: {detection_rate*100:.1f}% (target: >70%) ❌")
        if consistency_rate <= 0.95:
            print(f"   - Consistency rate: {consistency_rate*100:.1f}% (target: >95%) ❌")
        return 1


if __name__ == "__main__":
    sys.exit(main())
