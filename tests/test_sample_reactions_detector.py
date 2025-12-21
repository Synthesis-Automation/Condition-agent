#!/usr/bin/env python3
"""
Test sample reactions with the updated reaction type detector.
"""
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from sample_reactions import SAMPLE_REACTIONS
from chemtools import detect_reaction


def test_sample_reactions():
    """Test sample reactions with the updated detector."""
    # Parse the sample reactions (skip the 'Select...' placeholder)
    reactions = []
    for entry in SAMPLE_REACTIONS[1:100]:  # Test first ~100 reactions
        if '>>' in entry:
            # Extract SMILES (before the parenthesis with description)
            parts = entry.split(' (')
            smiles = parts[0].strip()
            desc = parts[1].rstrip(')') if len(parts) > 1 else ''
            reactions.append((smiles, desc))

    print(f'Testing {len(reactions)} reactions with updated detector...')
    print('=' * 80)

    results = {'correct': 0, 'incorrect': 0, 'unknown': 0}
    errors = []
    details = []

    for smiles, desc in reactions:
        try:
            r = detect_reaction(smiles, use_ml=False)  # Fast rule-based only
            family = r.get('family', 'Unknown')
            conf = r.get('confidence', 0)
            method = r.get('method', '')
            
            # Check if detection matches expected (from description)
            expected = None
            if 'Suzuki' in desc:
                expected = 'suzuki_miyaura'
            elif 'Sonogashira' in desc:
                expected = 'sonogashira'
            elif 'Heck' in desc:
                expected = 'heck'
            elif 'Negishi' in desc:
                expected = 'negishi'
            elif 'Kumada' in desc:
                expected = 'kumada'
            elif 'Stille' in desc:
                expected = 'stille'
            elif 'Buchwald' in desc or 'B-H' in desc:
                expected = 'c_n_cross_coupling'
            elif 'C-N' in desc:
                expected = 'c_n_cross_coupling'
            
            # Check match
            is_correct = False
            if expected:
                is_correct = (family == expected)
            
            status = '✓' if is_correct else ('?' if not expected else '✗')
            if is_correct:
                results['correct'] += 1
            elif expected:
                results['incorrect'] += 1
                errors.append((desc[:50], expected, family, conf))
            else:
                results['unknown'] += 1
            
            details.append({
                'desc': desc[:50],
                'expected': expected,
                'detected': family,
                'confidence': conf,
                'status': status
            })
                
        except Exception as e:
            print(f'ERROR: {desc[:40]} - {e}')

    # Print summary
    total_classified = results['correct'] + results['incorrect']
    accuracy = results['correct'] / total_classified * 100 if total_classified > 0 else 0
    
    print()
    print('=' * 80)
    print('SUMMARY')
    print('=' * 80)
    print(f"  Correct:   {results['correct']}")
    print(f"  Incorrect: {results['incorrect']}")
    print(f"  Unknown:   {results['unknown']}")
    print(f"  Accuracy:  {accuracy:.1f}%")
    print()

    if errors:
        print('MISCLASSIFIED REACTIONS:')
        print('-' * 80)
        for desc, exp, got, conf in errors[:15]:
            print(f'  {desc}')
            print(f'    Expected: {exp}, Got: {got} (conf: {conf:.2f})')
        print()

    # Print sample of correct detections
    print('SAMPLE CORRECT DETECTIONS:')
    print('-' * 80)
    correct_samples = [d for d in details if d['status'] == '✓'][:10]
    for d in correct_samples:
        print(f"  ✓ {d['desc']}")
        print(f"    Detected: {d['detected']} (conf: {d['confidence']:.2f})")
    
    return results


if __name__ == '__main__':
    test_sample_reactions()
