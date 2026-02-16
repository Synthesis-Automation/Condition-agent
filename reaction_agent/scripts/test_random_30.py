#!/usr/bin/env python
"""
Randomly pick 30 reactions from test_reactions_random_sampled.csv
and test them with the reaction analysis agent.

This is a more comprehensive test than the 10-reaction version.
"""

import os
import sys
import csv
import random
from pathlib import Path
from datetime import datetime
import json

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer, analyze_deterministic


class Colors:
    """ANSI color codes."""
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    END = '\033[0m'


def load_all_reactions(csv_path):
    """Load all reactions from CSV."""
    reactions = []
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            reactions.append(row)
    return reactions


def test_reaction(rxn_data, model, max_tokens=3000, show_details=False):
    """Test a single reaction."""

    rxn_smiles = rxn_data['reaction_smiles']
    rxn_id = rxn_data.get('reaction_id', 'unknown')
    rxn_type = rxn_data.get('reaction_type', 'unknown')

    if show_details:
        print(f"\n{'─'*80}")
        print(f"Testing: {rxn_id}")
        print(f"Type: {rxn_type}")
        print(f"SMILES: {rxn_smiles[:80]}...")
        print(f"{'─'*80}")

    try:
        # Step 1: Deterministic analysis
        det_result = analyze_deterministic(rxn_smiles, skip_mapping=False)

        mapping_qc = det_result.get('tool_facts', {}).get('mapping_qc', {})
        mapping_conf = mapping_qc.get('confidence', 0.0)
        mapping_ok = mapping_qc.get('ok', False)
        num_bond_changes = len(det_result.get('tool_facts', {}).get('bond_changes', []))

        if show_details:
            print(f"\n[1] Deterministic: mapping={mapping_conf:.3f} {'✓' if mapping_ok else '✗'}, bonds={num_bond_changes}")

        # Step 2: LLM analysis
        client = LLMClient(provider="openai", model=model)
        analyzer = ReactionSMILESAnalyzer(client, max_tokens=max_tokens)

        result = analyzer.analyze(rxn_smiles)

        interp = result.get('interpretation', {})
        overall_class = interp.get('overall_class', 'unknown')
        tags = interp.get('tags', [])
        confidence = interp.get('confidence', 0.0)
        num_events = len(interp.get('events', []))
        warnings = interp.get('warnings', [])

        if show_details:
            print(f"[2] LLM: class={overall_class}, conf={confidence:.3f}, events={num_events}")

        # Status
        status = "✓" if mapping_ok and confidence > 0.5 else "⚠" if confidence > 0.3 else "✗"

        print(f"{status} {rxn_id:<25} | map={mapping_conf:.3f} llm={confidence:.3f} | {rxn_type[:40]}")

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
            'status': status
        }

    except Exception as e:
        print(f"✗ {rxn_id:<25} | ERROR: {str(e)[:50]}")
        return {
            'reaction_id': rxn_id,
            'reaction_type': rxn_type,
            'success': False,
            'error': str(e),
            'status': '✗'
        }


def main():
    """Main testing function."""

    if not os.getenv("OPENAI_API_KEY"):
        print(f"{Colors.RED}Error: Set OPENAI_API_KEY to run test{Colors.END}")
        return

    # Load all reactions
    csv_path = project_root / "examples" / "test_reactions_random_sampled.csv"

    if not csv_path.exists():
        print(f"{Colors.RED}Error: {csv_path} not found{Colors.END}")
        return

    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}  Random 30 Reactions Test{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}\n")

    all_reactions = load_all_reactions(csv_path)
    print(f"Total reactions available: {len(all_reactions)}")

    # Randomly pick 30
    random.seed(42)  # For reproducibility
    selected = random.sample(all_reactions, min(30, len(all_reactions)))

    print(f"Randomly selected {len(selected)} reactions")

    # Show reaction type distribution
    type_counts = {}
    for rxn in selected:
        rtype = rxn['reaction_type']
        type_counts[rtype] = type_counts.get(rtype, 0) + 1

    print(f"\nReaction type distribution:")
    for rtype, count in sorted(type_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {count:2d}x {rtype}")
    if len(type_counts) > 10:
        print(f"  ... and {len(type_counts)-10} more types")

    # Test configuration
    model = "gpt-4o"  # Good balance
    max_tokens = 3000

    print(f"\n{Colors.BOLD}Testing with:{Colors.END} {model} ({max_tokens} tokens)")
    print(f"\n{Colors.CYAN}{'─'*80}{Colors.END}")
    print(f"{'Status':<3} {'ID':<25} | {'Confidence':<24} | {'Type'}")
    print(f"{Colors.CYAN}{'─'*80}{Colors.END}")

    results = []
    for i, rxn in enumerate(selected, 1):
        result = test_reaction(rxn, model, max_tokens, show_details=False)
        results.append(result)

        # Progress indicator every 10 reactions
        if i % 10 == 0:
            print(f"{Colors.CYAN}{'─'*80}{Colors.END}")
            print(f"{Colors.BOLD}Progress: {i}/{len(selected)} reactions tested{Colors.END}")
            print(f"{Colors.CYAN}{'─'*80}{Colors.END}")

    # Detailed Summary
    print(f"\n\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}  SUMMARY{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}\n")

    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]

    print(f"{Colors.BOLD}Overall:{Colors.END}")
    print(f"  Total tested:  {len(results)}")
    print(f"  Successful:    {Colors.GREEN}{len(successful)}{Colors.END}")
    print(f"  Failed:        {Colors.RED}{len(failed)}{Colors.END}")

    if successful:
        # Mapping confidence distribution
        high_map = sum(1 for r in successful if r['mapping_confidence'] >= 0.7)
        med_map = sum(1 for r in successful if 0.4 <= r['mapping_confidence'] < 0.7)
        low_map = sum(1 for r in successful if r['mapping_confidence'] < 0.4)

        print(f"\n{Colors.BOLD}Mapping Confidence:{Colors.END}")
        print(f"  High (≥0.7):   {Colors.GREEN}{high_map:2d}{Colors.END} ({high_map/len(successful)*100:.0f}%)")
        print(f"  Medium (0.4-0.7): {Colors.YELLOW}{med_map:2d}{Colors.END} ({med_map/len(successful)*100:.0f}%)")
        print(f"  Low (<0.4):    {Colors.RED}{low_map:2d}{Colors.END} ({low_map/len(successful)*100:.0f}%)")

        # LLM confidence distribution
        high_llm = sum(1 for r in successful if r['llm_confidence'] >= 0.7)
        med_llm = sum(1 for r in successful if 0.5 <= r['llm_confidence'] < 0.7)
        low_llm = sum(1 for r in successful if r['llm_confidence'] < 0.5)

        print(f"\n{Colors.BOLD}LLM Confidence:{Colors.END}")
        print(f"  High (≥0.7):   {Colors.GREEN}{high_llm:2d}{Colors.END} ({high_llm/len(successful)*100:.0f}%)")
        print(f"  Medium (0.5-0.7): {Colors.YELLOW}{med_llm:2d}{Colors.END} ({med_llm/len(successful)*100:.0f}%)")
        print(f"  Low (<0.5):    {Colors.RED}{low_llm:2d}{Colors.END} ({low_llm/len(successful)*100:.0f}%)")

        # Averages
        avg_mapping = sum(r['mapping_confidence'] for r in successful) / len(successful)
        avg_llm_conf = sum(r['llm_confidence'] for r in successful) / len(successful)
        avg_events = sum(r['num_events'] for r in successful) / len(successful)
        avg_bonds = sum(r['num_bond_changes'] for r in successful) / len(successful)

        print(f"\n{Colors.BOLD}Averages:{Colors.END}")
        print(f"  Mapping confidence:  {avg_mapping:.3f}")
        print(f"  LLM confidence:      {avg_llm_conf:.3f}")
        print(f"  Events per reaction: {avg_events:.1f}")
        print(f"  Bond changes:        {avg_bonds:.1f}")

        # Correlation
        print(f"\n{Colors.BOLD}Mapping vs LLM Confidence:{Colors.END}")

        # Group by mapping confidence
        high_map_rxns = [r for r in successful if r['mapping_confidence'] >= 0.7]
        low_map_rxns = [r for r in successful if r['mapping_confidence'] < 0.4]

        if high_map_rxns:
            avg_llm_high_map = sum(r['llm_confidence'] for r in high_map_rxns) / len(high_map_rxns)
            print(f"  High mapping (≥0.7) → LLM avg: {avg_llm_high_map:.3f}")

        if low_map_rxns:
            avg_llm_low_map = sum(r['llm_confidence'] for r in low_map_rxns) / len(low_map_rxns)
            print(f"  Low mapping (<0.4)  → LLM avg: {avg_llm_low_map:.3f}")

        # Problematic reactions (candidates for LLM-assisted mapping)
        problematic = [r for r in successful if r['mapping_confidence'] < 0.4]
        if problematic:
            print(f"\n{Colors.BOLD}Problematic Reactions:{Colors.END} (mapping <0.4, candidates for LLM-assisted analysis)")
            for r in problematic[:5]:  # Show first 5
                print(f"  {r['reaction_id']:<25} | map={r['mapping_confidence']:.3f} | {r['reaction_type'][:40]}")
            if len(problematic) > 5:
                print(f"  ... and {len(problematic)-5} more")

    if failed:
        print(f"\n{Colors.BOLD}{Colors.RED}Failed Reactions:{Colors.END}")
        for r in failed:
            print(f"  {r['reaction_id']}: {r.get('error', 'unknown error')[:60]}")

    # Save results
    output_dir = project_root / "reaction_agent" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    output_file = output_dir / f"random_30_test_{timestamp}.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump({
            'timestamp': timestamp,
            'model': model,
            'max_tokens': max_tokens,
            'total_tested': len(results),
            'successful': len(successful),
            'failed': len(failed),
            'summary': {
                'avg_mapping_confidence': avg_mapping if successful else 0,
                'avg_llm_confidence': avg_llm_conf if successful else 0,
                'high_mapping_count': high_map if successful else 0,
                'medium_mapping_count': med_map if successful else 0,
                'low_mapping_count': low_map if successful else 0,
                'high_llm_count': high_llm if successful else 0,
                'medium_llm_count': med_llm if successful else 0,
                'low_llm_count': low_llm if successful else 0
            },
            'results': results
        }, f, indent=2, ensure_ascii=False)

    print(f"\n\n{Colors.GREEN}✓ Results saved to:{Colors.END} {output_file}")
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{Colors.YELLOW}Interrupted by user{Colors.END}")
    except Exception as e:
        print(f"\n{Colors.RED}✗ Error: {e}{Colors.END}")
        import traceback
        traceback.print_exc()
