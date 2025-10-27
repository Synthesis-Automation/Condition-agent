"""
Test Analysis Module with Sample Reactions
==========================================

This script tests the chemtools.analysis module components with reactions
from sample_reactions.py to validate:
1. Reaction type detection
2. Reactant type classification  
3. Context-aware classification (Two-Pass Approach)
4. Family resolution and alias handling
"""

import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

# Add tests directory to path
sys.path.insert(0, str(Path(__file__).parent / "tests"))

# Import sample reactions
from sample_reactions import SAMPLE_REACTIONS

# Import analysis modules
from chemtools.analysis.reactions import (
    canonical_family_label,
    resolve_reaction_family,
    apply_catalyst_override,
)
from chemtools.analysis.reactants import iter_reactant_matches
from chemtools.analysis.reaction_context import (
    classify_reactants_with_context,
    get_reactant_summary,
)
from chemtools.router import detect_family_from_reaction
from chemtools.smiles import normalize_reaction, split_reaction


def extract_reaction_info(sample: str) -> Tuple[str, str, str]:
    """Extract SMILES, expected type, and description from sample."""
    if " (" in sample:
        smiles_part = sample.split(" (")[0].strip()
        desc_part = sample.split(" (")[1].rstrip(")")
        
        # Try to extract expected reaction type from description
        desc_lower = desc_part.lower()
        expected_type = "unknown"
        
        # Map description keywords to canonical reaction types
        if "suzuki" in desc_lower:
            expected_type = "suzuki_miyaura"
        elif "stille" in desc_lower:
            expected_type = "stille"
        elif "sonogashira" in desc_lower:
            expected_type = "sonogashira"
        elif "heck" in desc_lower:
            expected_type = "heck"
        elif "negishi" in desc_lower:
            expected_type = "negishi"
        elif "kumada" in desc_lower:
            expected_type = "kumada"
        elif "buchwald" in desc_lower or "b-h" in desc_lower:
            expected_type = "buchwald_hartwig_c_n"
        elif "c-n" in desc_lower or "ullmann_cn" in desc_lower:
            expected_type = "cn_coupling"  # Generic
        elif "chan" in desc_lower or "lam" in desc_lower:
            expected_type = "chan_lam"
        elif "ullmann ether" in desc_lower or "c-o" in desc_lower:
            expected_type = "ullmann_ether"
        elif "esterification" in desc_lower:
            expected_type = "esterification"
        elif "amidation" in desc_lower:
            expected_type = "amidation"
        elif "hydrogenation" in desc_lower:
            expected_type = "hydrogenation"
        elif "nabh4" in desc_lower:
            expected_type = "reduction"
        elif "lialh4" in desc_lower:
            expected_type = "reduction"
        elif "birch" in desc_lower:
            expected_type = "birch_reduction"
        elif "transfer hydrogenation" in desc_lower:
            expected_type = "transfer_hydrogenation"
        
        return smiles_part, expected_type, desc_part
    else:
        return sample.strip(), "unknown", ""


def test_reaction_detection():
    """Test 1: Reaction type detection accuracy."""
    print("=" * 80)
    print("TEST 1: Reaction Type Detection")
    print("=" * 80)
    
    stats = {
        "total": 0,
        "detected": 0,
        "correct": 0,
        "by_method": defaultdict(int),
        "by_expected_type": defaultdict(lambda: {"total": 0, "correct": 0}),
    }
    
    failed_detections = []
    mismatches = []
    
    for sample in SAMPLE_REACTIONS[1:]:  # Skip "Select a sample reaction..."
        smiles, expected_type, desc = extract_reaction_info(sample)
        
        if not smiles or expected_type == "unknown":
            continue
        
        stats["total"] += 1
        stats["by_expected_type"][expected_type]["total"] += 1
        
        try:
            # Normalize reaction SMILES
            normalized = normalize_reaction(smiles)
            
            # Detect reaction family
            result = detect_family_from_reaction(normalized)
            
            if result and result.get("family"):
                detected_family = result["family"]
                confidence = result.get("confidence", 0.0)
                method = result.get("method", "unknown")
                
                stats["detected"] += 1
                stats["by_method"][method] += 1
                
                # Check if detection matches expected
                if detected_family == expected_type:
                    stats["correct"] += 1
                    stats["by_expected_type"][expected_type]["correct"] += 1
                else:
                    # Check if it's an alias match
                    canonical = canonical_family_label(detected_family)
                    if canonical == expected_type:
                        stats["correct"] += 1
                        stats["by_expected_type"][expected_type]["correct"] += 1
                    else:
                        mismatches.append({
                            "smiles": smiles[:60] + "..." if len(smiles) > 60 else smiles,
                            "expected": expected_type,
                            "detected": detected_family,
                            "confidence": confidence,
                            "method": method,
                            "desc": desc[:50] + "..." if len(desc) > 50 else desc,
                        })
            else:
                failed_detections.append({
                    "smiles": smiles[:60] + "..." if len(smiles) > 60 else smiles,
                    "expected": expected_type,
                    "desc": desc[:50] + "..." if len(desc) > 50 else desc,
                })
        except Exception as e:
            print(f"  ⚠ Error processing: {smiles[:50]}... - {e}")
    
    # Print summary
    print(f"\nTotal reactions tested: {stats['total']}")
    print(f"Successfully detected: {stats['detected']} ({stats['detected']/stats['total']*100:.1f}%)")
    print(f"Correctly identified: {stats['correct']} ({stats['correct']/stats['total']*100:.1f}%)")
    print(f"\nDetection by method:")
    for method, count in stats["by_method"].items():
        print(f"  {method}: {count}")
    
    print(f"\nAccuracy by reaction type:")
    for rxn_type, counts in sorted(stats["by_expected_type"].items()):
        if counts["total"] > 0:
            accuracy = counts["correct"] / counts["total"] * 100
            print(f"  {rxn_type}: {counts['correct']}/{counts['total']} ({accuracy:.1f}%)")
    
    if failed_detections:
        print(f"\nFailed detections ({len(failed_detections)}):")
        for item in failed_detections[:5]:  # Show first 5
            print(f"  • {item['desc']}")
            print(f"    Expected: {item['expected']}")
            print(f"    SMILES: {item['smiles']}")
    
    if mismatches:
        print(f"\nMismatches ({len(mismatches)}):")
        for item in mismatches[:5]:  # Show first 5
            print(f"  • {item['desc']}")
            print(f"    Expected: {item['expected']}, Got: {item['detected']} ({item['confidence']:.2f}, {item['method']})")
            print(f"    SMILES: {item['smiles']}")
    
    return stats


def test_reactant_classification():
    """Test 2: Reactant type classification."""
    print("\n" + "=" * 80)
    print("TEST 2: Reactant Type Classification")
    print("=" * 80)
    
    stats = {
        "total_reactions": 0,
        "total_reactants": 0,
        "reactants_classified": 0,
        "multi_functional": 0,
    }
    
    # Test a subset of reactions
    test_samples = [
        # Suzuki examples
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Simple Ph-Ph)",
        "Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1 (Suzuki - Electron-poor ArCl)",
        
        # C-N examples
        "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (Buchwald-Hartwig - Diphenylamine)",
        "Clc1ccncc1.NCC>>CCNc1ccncc1 (Buchwald-Hartwig - Pyridine ethylamine)",
        
        # Hydrogenation
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1 (Hydrogenation - Ethylbenzene)",
        "c1ccc([N+](=O)[O-])cc1.[H][H]>>c1ccc(N)cc1 (Hydrogenation - Nitro to amine)",
    ]
    
    print("\nTesting reactant classification on selected reactions:\n")
    
    for sample in test_samples:
        smiles, expected_type, desc = extract_reaction_info(sample)
        
        if not smiles:
            continue
        
        stats["total_reactions"] += 1
        
        try:
            # Split reaction into reactants
            normalized = normalize_reaction(smiles)
            parts = split_reaction(normalized)
            reactants = parts.get("reactants", [])
            
            print(f"  {desc[:70]}")
            print(f"    Reactants: {len(reactants)}")
            
            for i, reactant in enumerate(reactants, 1):
                stats["total_reactants"] += 1
                
                # Classify reactant
                matches = list(iter_reactant_matches(reactant))
                
                if matches:
                    stats["reactants_classified"] += 1
                    match_strs = [f"{m.reactant_type_id}" for m in matches[:3]]
                    print(f"      [{i}] {reactant[:40]}: {', '.join(match_strs)}")
                    
                    if len(matches) > 1:
                        stats["multi_functional"] += 1
                else:
                    print(f"      [{i}] {reactant[:40]}: No matches")
            
            print()
            
        except Exception as e:
            print(f"    ⚠ Error: {e}\n")
    
    # Print summary
    print(f"Summary:")
    print(f"  Reactions tested: {stats['total_reactions']}")
    print(f"  Total reactants: {stats['total_reactants']}")
    print(f"  Reactants classified: {stats['reactants_classified']}")
    print(f"  Multi-functional reactants: {stats['multi_functional']}")
    
    if stats['total_reactants'] > 0:
        classification_rate = stats['reactants_classified'] / stats['total_reactants'] * 100
        print(f"  Classification rate: {classification_rate:.1f}%")
    
    return stats


def test_context_aware_classification():
    """Test 3: Context-aware classification (Two-Pass Approach)."""
    print("\n" + "=" * 80)
    print("TEST 3: Context-Aware Classification (Two-Pass Approach)")
    print("=" * 80)
    
    # Test with user-provided reaction types
    test_cases = [
        {
            "smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
            "reaction_type": "suzuki_miyaura",
            "desc": "Suzuki - Simple Ph-Ph (user-provided type)",
        },
        {
            "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "reaction_type": "buchwald_hartwig_c_n",
            "desc": "Buchwald-Hartwig - Diphenylamine (user-provided type)",
        },
        {
            "smiles": "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
            "reaction_type": None,  # Auto-detect
            "desc": "Sonogashira - Auto-detect",
        },
        {
            "smiles": "C=Cc1ccccc1.[H][H]>>CCc1ccccc1",
            "reaction_type": None,  # Auto-detect
            "desc": "Hydrogenation - Auto-detect",
        },
    ]
    
    print("\nTesting context-aware classification:\n")
    
    for i, case in enumerate(test_cases, 1):
        print(f"  [{i}] {case['desc']}")
        
        try:
            # Classify with context
            result = classify_reactants_with_context(
                case["smiles"],
                reaction_type=case["reaction_type"],
                auto_detect=(case["reaction_type"] is None),
            )
            
            print(f"      Detected: {result.reaction_type} ({result.detection_method}, conf={result.reaction_confidence:.2f})")
            print(f"      Reactants found: {len(result.reactants)}")
            
            for reactant in result.reactants:
                role_str = f" [{reactant.role}]" if reactant.role else ""
                expected_str = " ✓" if reactant.is_expected else " ?"
                print(f"        • {reactant.match.reactant_type_id}{role_str}{expected_str}")
                
                if reactant.alternative_matches:
                    alt_strs = [m.reactant_type_id for m in reactant.alternative_matches[:2]]
                    print(f"          Alternatives: {', '.join(alt_strs)}")
            
            if result.has_multi_functional:
                print(f"      ⚠ Contains multi-functional reactants")
            
            print()
            
        except Exception as e:
            print(f"      ⚠ Error: {e}\n")
    
    print("Context-aware classification test complete.")


def test_family_resolution():
    """Test 4: Family resolution and alias handling."""
    print("\n" + "=" * 80)
    print("TEST 4: Family Resolution and Alias Handling")
    print("=" * 80)
    
    test_aliases = [
        "Suzuki",
        "Suzuki-Miyaura",
        "Suzuki_CC",
        "suzuki_miyaura",
        "Buchwald-Hartwig",
        "Buchwald_CN",
        "buchwald_hartwig_c_n",
        "Chan_Lam_CN",
        "chan_lam_cn",
        "chan_lam",
        "C_N_Coupling",
        "C_O_Coupling",
        "Ullmann_CN",
        "Ullmann_CO",
        "Sonogashira",
        "Sonogashira_CC",
        "Heck",
        "Negishi",
        "Stille",
        "Kumada",
    ]
    
    print("\nTesting alias resolution:\n")
    
    for alias in test_aliases:
        canonical = canonical_family_label(alias)
        resolved = resolve_reaction_family(alias)
        
        if canonical or resolved:
            result = canonical or resolved or "No match"
            print(f"  {alias:25} → {result}")
        else:
            print(f"  {alias:25} → ⚠ No resolution")
    
    print("\nAlias resolution test complete.")


def main():
    """Run all tests."""
    print("\n" + "=" * 80)
    print("TESTING ANALYSIS MODULE WITH SAMPLE REACTIONS")
    print("=" * 80)
    print(f"Total sample reactions available: {len(SAMPLE_REACTIONS) - 1}")
    print()
    
    try:
        # Test 1: Reaction detection
        detection_stats = test_reaction_detection()
        
        # Test 2: Reactant classification
        reactant_stats = test_reactant_classification()
        
        # Test 3: Context-aware classification
        test_context_aware_classification()
        
        # Test 4: Family resolution
        test_family_resolution()
        
        # Final summary
        print("\n" + "=" * 80)
        print("OVERALL SUMMARY")
        print("=" * 80)
        print(f"✓ Reaction detection test complete")
        print(f"  - Accuracy: {detection_stats['correct']}/{detection_stats['total']} ({detection_stats['correct']/detection_stats['total']*100:.1f}%)")
        print(f"✓ Reactant classification test complete")
        print(f"  - Classification rate: {reactant_stats['reactants_classified']}/{reactant_stats['total_reactants']}")
        print(f"✓ Context-aware classification test complete")
        print(f"✓ Family resolution test complete")
        print()
        print("⚠ Note: SMARTS pattern matching issues may affect reactant classification results")
        print("  (See diagnostic check for details on functional group detection)")
        print()
        
    except Exception as e:
        print(f"\n⚠ Fatal error during testing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
