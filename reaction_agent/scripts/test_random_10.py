#!/usr/bin/env python
"""
Randomly pick 10 reactions from test_reactions_random_sampled.csv
and test them with the reaction analysis agent.
"""

import os
import sys
import csv
import random
from pathlib import Path
from datetime import datetime

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer, analyze_deterministic


def load_all_reactions(csv_path):
    """Load all reactions from CSV."""
    reactions = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            reactions.append(row)
    return reactions


def test_reaction(rxn_data, model, max_tokens=3000):
    """Test a single reaction."""

    rxn_smiles = rxn_data['reaction_smiles']
    rxn_id = rxn_data.get('reaction_id', 'unknown')
    rxn_type = rxn_data.get('reaction_type', 'unknown')

    print(f"\n{'─'*80}")
    print(f"Testing: {rxn_id}")
    print(f"Type: {rxn_type}")
    print(f"SMILES: {rxn_smiles[:80]}...")
    print(f"{'─'*80}")

    try:
        # Step 1: Deterministic analysis
        print("\n[1] Deterministic analysis...")
        det_result = analyze_deterministic(rxn_smiles, skip_mapping=False)

        mapping_qc = det_result.get('tool_facts', {}).get('mapping_qc', {})
        mapping_conf = mapping_qc.get('confidence', 0.0)
        mapping_ok = mapping_qc.get('ok', False)
        num_bond_changes = len(det_result.get('tool_facts', {}).get('bond_changes', []))

        print(f"    Mapping: {mapping_conf:.3f} {'✓' if mapping_ok else '✗'}")
        print(f"    Bond changes: {num_bond_changes}")

        # Step 2: LLM analysis
        print(f"\n[2] LLM analysis ({model})...")
        client = LLMClient(provider="openai", model=model)
        analyzer = ReactionSMILESAnalyzer(client, max_tokens=max_tokens)

        result = analyzer.analyze(rxn_smiles)

        interp = result.get('interpretation', {})
        overall_class = interp.get('overall_class', 'unknown')
        tags = interp.get('tags', [])
        confidence = interp.get('confidence', 0.0)
        num_events = len(interp.get('events', []))
        warnings = interp.get('warnings', [])

        print(f"    Classification: {overall_class}")
        print(f"    Tags: {', '.join(tags[:3])}")
        print(f"    Confidence: {confidence:.3f}")
        print(f"    Events: {num_events}")
        print(f"    Warnings: {len(warnings)}")

        # Summary
        status = "✓" if mapping_ok and confidence > 0.5 else "⚠"
        print(f"\n{status} Result: mapping={mapping_conf:.3f}, conf={confidence:.3f}")

        return {
            'reaction_id': rxn_id,
            'reaction_type': rxn_type,
            'success': True,
            'mapping_confidence': mapping_conf,
            'mapping_ok': mapping_ok,
            'num_bond_changes': num_bond_changes,
            'overall_class': overall_class,
            'tags': tags,
            'llm_confidence': confidence,
            'num_events': num_events,
            'num_warnings': len(warnings),
            'full_result': result
        }

    except Exception as e:
        print(f"\n✗ Error: {e}")
        return {
            'reaction_id': rxn_id,
            'reaction_type': rxn_type,
            'success': False,
            'error': str(e)
        }


def main():
    """Main testing function."""

    if not os.getenv("OPENAI_API_KEY"):
        print("Error: Set OPENAI_API_KEY to run test")
        return

    # Load all reactions
    csv_path = project_root / "examples" / "test_reactions_random_sampled.csv"

    if not csv_path.exists():
        print(f"Error: {csv_path} not found")
        return

    print("="*80)
    print("  Random 10 Reactions Test")
    print("="*80)

    all_reactions = load_all_reactions(csv_path)
    print(f"\nTotal reactions available: {len(all_reactions)}")

    # Randomly pick 10
    random.seed(42)  # For reproducibility
    selected = random.sample(all_reactions, min(10, len(all_reactions)))

    print(f"\nRandomly selected {len(selected)} reactions")
    print("\nSelected reactions:")
    for i, rxn in enumerate(selected, 1):
        print(f"  {i}. {rxn['reaction_id']} - {rxn['reaction_type']}")

    # Test each reaction
    model = "gpt-4o"  # Good balance of quality/cost
    max_tokens = 3000

    print(f"\n\nTesting with model: {model} ({max_tokens} tokens)")
    input("\nPress Enter to start testing...")

    results = []
    for i, rxn in enumerate(selected, 1):
        print(f"\n\n{'='*80}")
        print(f"  Reaction {i}/{len(selected)}")
        print(f"{'='*80}")

        result = test_reaction(rxn, model, max_tokens)
        results.append(result)

    # Summary
    print(f"\n\n{'='*80}")
    print("  SUMMARY")
    print(f"{'='*80}\n")

    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]

    print(f"Successful: {len(successful)}/{len(results)}")
    print(f"Failed: {len(failed)}/{len(results)}")

    if successful:
        avg_mapping = sum(r['mapping_confidence'] for r in successful) / len(successful)
        avg_llm_conf = sum(r['llm_confidence'] for r in successful) / len(successful)
        high_conf = sum(1 for r in successful if r['llm_confidence'] >= 0.7)

        print(f"\nAverage mapping confidence: {avg_mapping:.3f}")
        print(f"Average LLM confidence: {avg_llm_conf:.3f}")
        print(f"High confidence results (≥0.7): {high_conf}/{len(successful)}")

        print("\n\nDetailed Results:")
        print(f"{'ID':<25} {'Type':<30} {'Map':<6} {'LLM':<6} {'Events':<7} {'Status':<8}")
        print("─"*90)
        for r in successful:
            status = "✓" if r['mapping_ok'] and r['llm_confidence'] > 0.5 else "⚠"
            rxn_type = r['reaction_type'][:28] + ".." if len(r['reaction_type']) > 30 else r['reaction_type']
            print(f"{r['reaction_id']:<25} {rxn_type:<30} {r['mapping_confidence']:.3f}  {r['llm_confidence']:.3f}  {r['num_events']:<7} {status}")

    if failed:
        print(f"\n\nFailed Reactions:")
        for r in failed:
            print(f"  - {r['reaction_id']}: {r.get('error', 'unknown error')}")

    # Save results
    output_dir = project_root / "reaction_agent" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    import json
    output_file = output_dir / f"random_10_test_{timestamp}.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump({
            'timestamp': timestamp,
            'model': model,
            'max_tokens': max_tokens,
            'total_tested': len(results),
            'successful': len(successful),
            'failed': len(failed),
            'results': results
        }, f, indent=2, ensure_ascii=False)

    print(f"\n\n✓ Results saved to: {output_file}")
    print(f"\n{'='*80}\n")


if __name__ == "__main__":
    main()
