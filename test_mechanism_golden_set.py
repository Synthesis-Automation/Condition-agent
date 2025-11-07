"""
Test the mechanism analyzer against the golden set of reactions.

This test suite validates the new mechanism analyzer tool using curated
reaction examples with known mechanism types.
"""

import json
from pathlib import Path
from typing import Dict, Any, List

from chem_assistant.chemtools_wrapper import analyze_mechanism_tool


def load_golden_set() -> List[Dict[str, Any]]:
    """Load the golden set test cases."""
    golden_file = Path("tests/data/mechanism_golden_set.json")
    if not golden_file.exists():
        golden_file = Path("tests") / "data" / "mechanism_golden_set.json"
    if not golden_file.exists():
        golden_file = Path(__file__).parent / "tests" / "data" / "mechanism_golden_set.json"
    
    with open(golden_file, 'r') as f:
        return json.load(f)


def normalize_mechanism_name(name: str) -> str:
    """Normalize mechanism names for comparison."""
    return name.lower().replace("-", "_").replace(" ", "_")


def test_mechanism_analyzer_golden_set():
    """Test mechanism analyzer against all golden set examples."""
    
    print("=" * 80)
    print("MECHANISM ANALYZER GOLDEN SET TESTS")
    print("=" * 80)
    
    golden_set = load_golden_set()
    
    results = {
        "total": len(golden_set),
        "passed": 0,
        "failed": 0,
        "details": []
    }
    
    for i, test_case in enumerate(golden_set, 1):
        name = test_case["name"]
        reaction_smiles = test_case["reaction_smiles"]
        expected_mechanism = test_case["expected_mechanism"]
        expected_family = test_case.get("expected_family")
        
        print(f"\n[{i}/{len(golden_set)}] Testing: {name}")
        print(f"   Reaction: {reaction_smiles}")
        print(f"   Expected mechanism: {expected_mechanism}")
        
        try:
            # Run the mechanism analyzer
            result = analyze_mechanism_tool.invoke({
                "reaction_smiles": reaction_smiles,
                "detail_level": "standard",
                "include_bond_changes": True,
                "include_electron_flow": True,
                "include_intermediates": True,
                "include_precedents": False  # Skip for speed
            })
            
            if not result.get("success"):
                print(f"   ❌ FAILED: Tool returned error: {result.get('error')}")
                results["failed"] += 1
                results["details"].append({
                    "name": name,
                    "status": "error",
                    "error": result.get("error"),
                    "expected": expected_mechanism
                })
                continue
            
            # Extract predicted mechanism
            predicted_mechanism = result.get("mechanism_type", "unknown")
            confidence = result.get("confidence", 0.0)
            description = result.get("description", "")
            
            # Normalize names for comparison
            predicted_norm = normalize_mechanism_name(predicted_mechanism)
            expected_norm = normalize_mechanism_name(expected_mechanism)
            
            # Check if prediction matches
            match = predicted_norm == expected_norm
            
            if match:
                print(f"   ✅ PASSED: {predicted_mechanism} (confidence: {confidence:.2f})")
                results["passed"] += 1
                status = "pass"
            else:
                print(f"   ⚠️  MISMATCH: Got '{predicted_mechanism}', expected '{expected_mechanism}'")
                print(f"      Confidence: {confidence:.2f}")
                results["failed"] += 1
                status = "mismatch"
            
            # Additional info
            print(f"   Description: {description[:80]}...")
            
            # Check electron flow
            electron_flow = result.get("electron_flow", {})
            arrows = electron_flow.get("arrows", [])
            if arrows:
                print(f"   Electron flow: {len(arrows)} step(s)")
                for arrow in arrows[:2]:  # Show first 2
                    print(f"      • {arrow.get('from')} → {arrow.get('to')}")
            
            # Check intermediates
            intermediates = result.get("intermediates", [])
            if intermediates:
                print(f"   Intermediates: {len(intermediates)} predicted")
                for intermed in intermediates[:2]:  # Show first 2
                    print(f"      • {intermed.get('name', 'unknown')}")
            
            # Store details
            results["details"].append({
                "name": name,
                "status": status,
                "predicted": predicted_mechanism,
                "expected": expected_mechanism,
                "confidence": confidence,
                "num_arrows": len(arrows),
                "num_intermediates": len(intermediates)
            })
            
        except Exception as e:
            print(f"   ❌ EXCEPTION: {e}")
            results["failed"] += 1
            results["details"].append({
                "name": name,
                "status": "exception",
                "error": str(e),
                "expected": expected_mechanism
            })
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total tests: {results['total']}")
    print(f"Passed: {results['passed']} ✅")
    print(f"Failed: {results['failed']} ❌")
    print(f"Success rate: {results['passed'] / results['total'] * 100:.1f}%")
    
    # Detailed breakdown
    print("\n" + "-" * 80)
    print("DETAILED RESULTS")
    print("-" * 80)
    for detail in results["details"]:
        status_emoji = {
            "pass": "✅",
            "mismatch": "⚠️ ",
            "error": "❌",
            "exception": "💥"
        }.get(detail["status"], "❓")
        
        print(f"{status_emoji} {detail['name']}")
        if detail["status"] == "pass":
            print(f"   Predicted: {detail['predicted']} (confidence: {detail.get('confidence', 0):.2f})")
        elif detail["status"] == "mismatch":
            print(f"   Predicted: {detail['predicted']} (confidence: {detail.get('confidence', 0):.2f})")
            print(f"   Expected:  {detail['expected']}")
        else:
            print(f"   Error: {detail.get('error', 'Unknown error')}")
    
    return results


def test_detailed_analysis():
    """Run detailed analysis on one example to show all features."""
    
    print("\n" + "=" * 80)
    print("DETAILED ANALYSIS EXAMPLE")
    print("=" * 80)
    
    # Use Buchwald-Hartwig as example
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    
    print(f"\nAnalyzing: Buchwald-Hartwig C-N coupling")
    print(f"Reaction: {reaction}")
    
    result = analyze_mechanism_tool.invoke({
        "reaction_smiles": reaction,
        "detail_level": "high",  # Get more info
        "include_bond_changes": True,
        "include_electron_flow": True,
        "include_intermediates": True,
        "include_precedents": False
    })
    
    if result.get("success"):
        print("\n✅ Analysis successful!\n")
        
        # Mechanism type
        print(f"Mechanism: {result.get('mechanism_type')}")
        print(f"Confidence: {result.get('confidence', 0):.2f}")
        print(f"Description: {result.get('description')}")
        
        # Steps
        steps = result.get("steps", [])
        if steps:
            print(f"\nMechanism Steps ({len(steps)}):")
            for i, step in enumerate(steps, 1):
                print(f"  {i}. {step.get('description', 'N/A')}")
        
        # Electron flow
        electron_flow = result.get("electron_flow", {})
        arrows = electron_flow.get("arrows", [])
        if arrows:
            print(f"\nElectron Flow ({len(arrows)} arrows):")
            for i, arrow in enumerate(arrows, 1):
                print(f"  {i}. {arrow.get('from')} → {arrow.get('to')}")
                print(f"     {arrow.get('description', '')}")
        
        # Intermediates
        intermediates = result.get("intermediates", [])
        if intermediates:
            print(f"\nIntermediates ({len(intermediates)}):")
            for i, intermed in enumerate(intermediates, 1):
                print(f"  {i}. {intermed.get('name', 'Unknown')}")
                print(f"     {intermed.get('description', '')}")
        
        # Narrative
        narrative = result.get("narrative", "")
        if narrative:
            print(f"\nNarrative:")
            print(f"  {narrative[:300]}...")
        
        # Highlights
        highlights = result.get("highlights", [])
        if highlights:
            print(f"\nKey Highlights:")
            for highlight in highlights:
                print(f"  • {highlight}")
        
        # Warnings
        warnings = result.get("warnings", [])
        if warnings:
            print(f"\nWarnings:")
            for warning in warnings:
                print(f"  ⚠️  {warning}")
        
    else:
        print(f"\n❌ Analysis failed: {result.get('error')}")


def test_electron_flow_coverage():
    """Test that all mechanism types have electron flow data."""
    
    print("\n" + "=" * 80)
    print("ELECTRON FLOW COVERAGE TEST")
    print("=" * 80)
    
    golden_set = load_golden_set()
    
    coverage = {
        "with_flow": 0,
        "without_flow": 0,
        "details": []
    }
    
    for test_case in golden_set:
        name = test_case["name"]
        reaction = test_case["reaction_smiles"]
        
        result = analyze_mechanism_tool.invoke({
            "reaction_smiles": reaction,
            "include_electron_flow": True
        })
        
        if result.get("success"):
            electron_flow = result.get("electron_flow", {})
            arrows = electron_flow.get("arrows", [])
            
            if arrows:
                coverage["with_flow"] += 1
                status = "✅"
            else:
                coverage["without_flow"] += 1
                status = "❌"
            
            coverage["details"].append({
                "name": name,
                "has_flow": bool(arrows),
                "num_arrows": len(arrows)
            })
            
            print(f"{status} {name}: {len(arrows)} arrow(s)")
    
    print(f"\n📊 Coverage: {coverage['with_flow']}/{len(golden_set)} reactions have electron flow data")
    print(f"   Coverage rate: {coverage['with_flow'] / len(golden_set) * 100:.1f}%")


def main():
    """Run all tests."""
    
    print("🧪 Starting Mechanism Analyzer Test Suite\n")
    
    # Test 1: Golden set validation
    results = test_mechanism_analyzer_golden_set()
    
    # Test 2: Detailed analysis example
    test_detailed_analysis()
    
    # Test 3: Electron flow coverage
    test_electron_flow_coverage()
    
    # Final summary
    print("\n" + "=" * 80)
    print("🎯 FINAL SUMMARY")
    print("=" * 80)
    
    success_rate = results["passed"] / results["total"] * 100
    
    if success_rate >= 85:
        grade = "🏆 EXCELLENT"
    elif success_rate >= 70:
        grade = "✅ GOOD"
    elif success_rate >= 50:
        grade = "⚠️  FAIR"
    else:
        grade = "❌ NEEDS WORK"
    
    print(f"\nOverall Grade: {grade}")
    print(f"Success Rate: {success_rate:.1f}%")
    print(f"Tests Passed: {results['passed']}/{results['total']}")
    
    if results["failed"] > 0:
        print(f"\n💡 Tip: Review mismatches to improve mechanism classification rules")
    
    print("\n✨ Test suite complete!\n")
    
    return results


if __name__ == "__main__":
    main()
