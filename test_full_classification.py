"""
Comprehensive test for both reaction type and reactant type identification.

Tests all 320 reactions from sample_reactions.py to evaluate:
1. Reaction type detection accuracy
2. Reactant type classification accuracy
3. Context-aware classification with Two-Pass Approach
"""

import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Any

# Add tests directory to path
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, 'tests'))

from sample_reactions import SAMPLE_REACTIONS
from chemtools.analysis.reaction_context import (
    classify_reactants_with_context,
    get_reactant_summary,
)
from chemtools.analysis.reactants import classify_reactant_smiles, iter_reactant_matches
from chemtools.router import detect_family_from_reaction


def print_header(title: str, char: str = "="):
    """Print formatted section header."""
    print(f"\n{char * 100}")
    print(title.upper().center(100))
    print(f"{char * 100}\n")


def print_divider(char: str = "-"):
    """Print a simple divider."""
    print(char * 100)


def test_reaction_and_reactant_classification():
    """Test both reaction type and reactant type classification for all samples."""
    
    print_header("COMPREHENSIVE REACTION AND REACTANT CLASSIFICATION TEST")
    
    print(f"Testing {len(SAMPLE_REACTIONS)} sample reactions...\n")
    
    # Statistics trackers
    reaction_stats = {
        'total': 0,
        'detected': 0,
        'unknown': 0,
        'families': Counter(),
        'confidence_high': 0,
        'confidence_med': 0,
        'confidence_low': 0,
    }
    
    reactant_stats = {
        'total_reactions': 0,
        'reactions_with_reactants': 0,
        'reactions_without_reactants': 0,
        'total_reactants_found': 0,
        'reactant_types': Counter(),
        'functional_groups': Counter(),
    }
    
    # Detailed results for analysis
    results: List[Dict[str, Any]] = []
    
    # Test each reaction
    for idx, sample in enumerate(SAMPLE_REACTIONS, 1):
        # Skip the first entry if it's the placeholder
        if isinstance(sample, str) and sample.startswith("Select a sample"):
            continue
            
        reaction_stats['total'] += 1
        reactant_stats['total_reactions'] += 1
        
        # Parse the reaction string
        # Format: "SMILES (Name)" or just "SMILES"
        if isinstance(sample, str):
            if '(' in sample and sample.endswith(')'):
                smiles, name_part = sample.rsplit('(', 1)
                smiles = smiles.strip()
                name = name_part.rstrip(')').strip()
            else:
                smiles = sample.strip()
                name = f"Reaction {idx}"
            expected = 'Unknown'  # Sample reactions don't have expected families
        else:
            # If it's a dict format (shouldn't happen with current data)
            name = sample.get('name', f'Reaction {idx}')
            smiles = sample.get('smiles', '')
            expected = sample.get('expected_family', 'Unknown')
        
        # 1. Detect reaction type
        detection_result = detect_family_from_reaction(smiles)
        detected_family = detection_result.get('family', 'Unknown')
        confidence = detection_result.get('confidence', 0.0)
        
        # Track reaction type stats
        if detected_family and detected_family != 'Unknown':
            reaction_stats['detected'] += 1
            reaction_stats['families'][detected_family] += 1
            
            if confidence >= 0.8:
                reaction_stats['confidence_high'] += 1
            elif confidence >= 0.5:
                reaction_stats['confidence_med'] += 1
            else:
                reaction_stats['confidence_low'] += 1
        else:
            reaction_stats['unknown'] += 1
            reaction_stats['families']['Unknown'] += 1
        
        # 2. Classify reactants WITHOUT context (baseline)
        # Extract individual reactant SMILES and classify each
        from chemtools.smiles import normalize_reaction
        normalized = normalize_reaction(smiles)
        baseline_reactants = []
        if normalized and 'reactants' in normalized:
            for r in normalized['reactants']:
                r_smiles = r.get('smiles_norm') or r.get('largest_smiles') or ''
                if r_smiles:
                    matches = list(iter_reactant_matches(r_smiles))
                    baseline_reactants.extend(matches)
        
        # 3. Classify reactants WITH context (Two-Pass Approach)
        # Use auto-detected reaction type
        context_result = classify_reactants_with_context(
            smiles,
            reaction_type=None,  # Auto-detect
            auto_detect=True
        )
        
        # Extract reactant info
        reactant_matches = context_result.reactants if context_result else []
        num_reactants = len(reactant_matches)
        
        if num_reactants > 0:
            reactant_stats['reactions_with_reactants'] += 1
            reactant_stats['total_reactants_found'] += num_reactants
            
            for match in reactant_matches:
                reactant_stats['reactant_types'][match.category] += 1
                # Track functional groups from the underlying match
                if hasattr(match.match, 'member_type'):
                    reactant_stats['functional_groups'][match.match.member_type] += 1
        else:
            reactant_stats['reactions_without_reactants'] += 1
        
        # Store detailed result
        result = {
            'index': idx,
            'name': name,
            'smiles': smiles,
            'expected_family': expected,
            'detected_family': detected_family,
            'confidence': confidence,
            'match': detected_family.lower() == expected.lower() if expected != 'Unknown' else None,
            'baseline_reactants': len(baseline_reactants),
            'context_reactants': num_reactants,
            'reactant_types': [m.category for m in reactant_matches] if reactant_matches else [],
            'detection_method': context_result.detection_method if context_result else 'unknown',
        }
        results.append(result)
        
        # Print progress every 50 reactions
        if idx % 50 == 0:
            print(f"Processed {idx}/{len(SAMPLE_REACTIONS)} reactions...")
    
    # Print summary statistics
    print_header("REACTION TYPE DETECTION SUMMARY")
    
    print(f"Total Reactions: {reaction_stats['total']}")
    print(f"Successfully Detected: {reaction_stats['detected']} ({reaction_stats['detected']/reaction_stats['total']*100:.1f}%)")
    print(f"Unknown/Failed: {reaction_stats['unknown']} ({reaction_stats['unknown']/reaction_stats['total']*100:.1f}%)")
    print()
    
    print("Confidence Distribution:")
    print(f"  High (>=0.8):   {reaction_stats['confidence_high']} ({reaction_stats['confidence_high']/reaction_stats['total']*100:.1f}%)")
    print(f"  Medium (>=0.5): {reaction_stats['confidence_med']} ({reaction_stats['confidence_med']/reaction_stats['total']*100:.1f}%)")
    print(f"  Low (<0.5):     {reaction_stats['confidence_low']} ({reaction_stats['confidence_low']/reaction_stats['total']*100:.1f}%)")
    print()
    
    print("Top 15 Detected Families:")
    for family, count in reaction_stats['families'].most_common(15):
        pct = count / reaction_stats['total'] * 100
        print(f"  {family:30s}: {count:3d} ({pct:5.1f}%)")
    
    # Print reactant classification summary
    print_header("REACTANT TYPE CLASSIFICATION SUMMARY")
    
    print(f"Total Reactions: {reactant_stats['total_reactions']}")
    print(f"Reactions with Reactants Found: {reactant_stats['reactions_with_reactants']} ({reactant_stats['reactions_with_reactants']/reactant_stats['total_reactions']*100:.1f}%)")
    print(f"Reactions with NO Reactants: {reactant_stats['reactions_without_reactants']} ({reactant_stats['reactions_without_reactants']/reactant_stats['total_reactions']*100:.1f}%)")
    print(f"Total Reactants Classified: {reactant_stats['total_reactants_found']}")
    
    if reactant_stats['reactions_with_reactants'] > 0:
        avg_reactants = reactant_stats['total_reactants_found'] / reactant_stats['reactions_with_reactants']
        print(f"Average Reactants per Reaction (when found): {avg_reactants:.2f}")
    print()
    
    if reactant_stats['reactant_types']:
        print("Top 15 Reactant Types Found:")
        for rtype, count in reactant_stats['reactant_types'].most_common(15):
            print(f"  {rtype:30s}: {count:3d}")
    else:
        print("WARNING: No reactant types were classified!")
        print("This indicates the SMARTS pattern matching issue is still present.")
    print()
    
    if reactant_stats['functional_groups']:
        print("Top 15 Functional Groups Detected:")
        for fg, count in reactant_stats['functional_groups'].most_common(15):
            print(f"  {fg:30s}: {count:3d}")
    
    # Print detailed examples
    print_header("DETAILED CLASSIFICATION EXAMPLES", "=")
    
    # Show 5 successful cases
    print("\n>>> SUCCESSFUL CLASSIFICATIONS (showing first 5)\n")
    successful = [r for r in results if r['context_reactants'] > 0]
    for result in successful[:5]:
        print_divider()
        print(f"Sample {result['index']}: {result['name']}")
        print(f"SMILES: {result['smiles'][:80]}...")
        print(f"Detected Family: {result['detected_family']} (confidence: {result['confidence']:.2f})")
        print(f"Reactants Found: {result['context_reactants']}")
        print(f"Reactant Types: {', '.join(result['reactant_types'])}")
        print(f"Detection Method: {result['detection_method']}")
    
    # Show 5 failed cases
    print("\n\n>>> FAILED CLASSIFICATIONS (showing first 5)\n")
    failed = [r for r in results if r['context_reactants'] == 0]
    for result in failed[:5]:
        print_divider()
        print(f"Sample {result['index']}: {result['name']}")
        print(f"SMILES: {result['smiles'][:80]}...")
        print(f"Detected Family: {result['detected_family']} (confidence: {result['confidence']:.2f})")
        print(f"Reactants Found: {result['context_reactants']} (FAILED)")
        print(f"Baseline Reactants: {result['baseline_reactants']}")
        print(f"Detection Method: {result['detection_method']}")
    
    # Accuracy analysis (when expected family is known)
    print_header("ACCURACY ANALYSIS", "=")
    
    with_expected = [r for r in results if r['expected_family'] != 'Unknown']
    if with_expected:
        correct = [r for r in with_expected if r['match'] is True]
        incorrect = [r for r in with_expected if r['match'] is False]
        
        print(f"Reactions with Expected Family: {len(with_expected)}")
        print(f"Correct Predictions: {len(correct)} ({len(correct)/len(with_expected)*100:.1f}%)")
        print(f"Incorrect Predictions: {len(incorrect)} ({len(incorrect)/len(with_expected)*100:.1f}%)")
        
        if incorrect:
            print("\nTop Misclassification Patterns:")
            misclass = Counter()
            for r in incorrect:
                misclass[(r['expected_family'], r['detected_family'])] += 1
            
            for (expected, detected), count in misclass.most_common(10):
                print(f"  Expected: {expected:25s} -> Detected: {detected:25s} ({count} cases)")
    
    # Two-Pass vs Baseline comparison
    print_header("TWO-PASS vs BASELINE COMPARISON", "=")
    
    baseline_total = sum(r['baseline_reactants'] for r in results)
    context_total = sum(r['context_reactants'] for r in results)
    
    print(f"Baseline Classification (no context):")
    print(f"  Total Reactants: {baseline_total}")
    print(f"  Reactions with Reactants: {len([r for r in results if r['baseline_reactants'] > 0])}")
    print()
    print(f"Two-Pass Classification (with context):")
    print(f"  Total Reactants: {context_total}")
    print(f"  Reactions with Reactants: {len([r for r in results if r['context_reactants'] > 0])}")
    print()
    
    if baseline_total == 0 and context_total == 0:
        print("WARNING: Both baseline and context-aware classification returned 0 reactants!")
        print("This confirms the SMARTS pattern matching issue affects the entire pipeline.")
    elif context_total > baseline_total:
        improvement = (context_total - baseline_total) / max(baseline_total, 1) * 100
        print(f"Two-Pass Approach found {improvement:.1f}% more reactants than baseline.")
    elif context_total < baseline_total:
        decrease = (baseline_total - context_total) / baseline_total * 100
        print(f"Two-Pass Approach found {decrease:.1f}% fewer reactants than baseline.")
    else:
        print("Both approaches found the same number of reactants.")
    
    print_header("TEST COMPLETED", "=")
    
    return results, reaction_stats, reactant_stats


def main():
    """Run all tests."""
    try:
        results, reaction_stats, reactant_stats = test_reaction_and_reactant_classification()
        
        # Save detailed results to file
        output_file = "full_classification_results.txt"
        print(f"\nSaving detailed results to {output_file}...")
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("DETAILED CLASSIFICATION RESULTS\n")
            f.write("=" * 100 + "\n\n")
            
            for result in results:
                f.write(f"Sample {result['index']}: {result['name']}\n")
                f.write(f"  SMILES: {result['smiles']}\n")
                f.write(f"  Expected: {result['expected_family']}\n")
                f.write(f"  Detected: {result['detected_family']} (confidence: {result['confidence']:.2f})\n")
                f.write(f"  Baseline Reactants: {result['baseline_reactants']}\n")
                f.write(f"  Context Reactants: {result['context_reactants']}\n")
                if result['reactant_types']:
                    f.write(f"  Reactant Types: {', '.join(result['reactant_types'])}\n")
                f.write(f"  Detection Method: {result['detection_method']}\n")
                f.write(f"  Match: {result['match']}\n")
                f.write("\n")
        
        print(f"OK - Results saved to {output_file}")
        
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
