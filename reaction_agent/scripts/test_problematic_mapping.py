#!/usr/bin/env python
"""
Test LLM-assisted mapping on the 3 problematic reactions from random test
that had poor rxnmapper confidence.
"""

import os
import sys
from pathlib import Path
import json

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from reaction_agent.scripts.llm_assisted_mapping import hybrid_mapping_workflow


def test_problematic_reactions():
    """Test hybrid mapping on reactions where rxnmapper failed."""

    if not os.getenv("OPENAI_API_KEY"):
        print("Set OPENAI_API_KEY to run test")
        return

    # The 3 problematic reactions from random_10 test
    reactions = [
        {
            "name": "Weinreb amide",
            "id": "31-049-CAS-19587793",
            "smiles": "CCCCCCCCc1cc[c]([Mg][Br])cc1.CON(C)C(=O)CC[C@@H]1C[C@@H](O[Si](C)(C)C(C)(C)C)CN1C(=O)OC(C)(C)C>>CCCCCCCCc1ccc(C(=O)CC[C@@H]2C[C@@H](O)CN2C(=O)OC(C)(C)C)cc1",
            "expected_mapping": 0.341
        },
        {
            "name": "Bucherer-Bergs",
            "id": "31-173-CAS-18125986",
            "smiles": "N#[C][Na].N.O=C(O)O.CCOc1ccc(C(C)=O)cc1>>CCOc1ccc(C2(C)NC(=O)NC2=O)cc1",
            "expected_mapping": 0.143
        },
        {
            "name": "Staudinger reduction",
            "id": "31-367-CAS-8556291",
            "smiles": "C[C@@H](C(=O)OCc1ccc(OCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)cc1)[C@H](N=[N+]=[N-])c1ccccc1.O=P(c1ccccc1)(c1ccccc1)c1ccccc1>>C[C@@H](C(=O)OCc1ccc(OCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)cc1)[C@H](N)c1ccccc1",
            "expected_mapping": 0.004
        }
    ]

    print("="*80)
    print("  LLM-Assisted Mapping: Problematic Reactions Test")
    print("="*80)
    print("\nTesting 3 reactions where rxnmapper had poor confidence")
    print(f"Expected mapping scores: {[r['expected_mapping'] for r in reactions]}")
    print("\nUsing o3-mini for detailed mechanistic analysis...")

    all_results = []

    for i, rxn in enumerate(reactions, 1):
        print(f"\n\n{'#'*80}")
        print(f"  Reaction {i}/3: {rxn['name']}")
        print(f"  ID: {rxn['id']}")
        print(f"  Expected mapping: {rxn['expected_mapping']:.3f}")
        print(f"{'#'*80}")

        result = hybrid_mapping_workflow(
            rxn['smiles'],
            confidence_threshold=0.6
        )

        all_results.append({
            'reaction_name': rxn['name'],
            'reaction_id': rxn['id'],
            'expected_mapping': rxn['expected_mapping'],
            'result': result
        })

    # Summary
    print(f"\n\n{'='*80}")
    print("  SUMMARY: LLM Analysis of Failed Mappings")
    print(f"{'='*80}\n")

    for res in all_results:
        name = res['reaction_name']
        expected = res['expected_mapping']
        final_conf = res['result']['final_confidence']
        needs_review = res['result'].get('needs_manual_review', False)

        print(f"\n{name} ({res['reaction_id']}):")
        print(f"  rxnmapper confidence: {expected:.3f} {'✗ FAILED' if expected < 0.4 else '⚠ LOW'}")
        print(f"  Final confidence: {final_conf:.3f}")
        print(f"  Manual review needed: {needs_review}")

        if 'llm_analysis' in res['result']:
            analysis = res['result']['llm_analysis']['llm_analysis']
            rxn_analysis = analysis.get('reaction_analysis', {})
            mapping_assessment = analysis.get('mapping_assessment', {})

            print(f"\n  LLM Analysis:")
            print(f"    Reaction type: {rxn_analysis.get('type', 'N/A')}")
            print(f"    Complexity: {rxn_analysis.get('complexity', 'N/A')}")

            stages = rxn_analysis.get('stages', [])
            if stages:
                print(f"    Stages identified: {len(stages)}")
                for stage in stages:
                    print(f"      - {stage.get('name', 'Unknown')}")

            errors = mapping_assessment.get('major_errors', [])
            if errors:
                print(f"\n  Mapping Errors Found ({len(errors)}):")
                for j, error in enumerate(errors[:3], 1):  # Show first 3
                    print(f"    {j}. {error[:80]}{'...' if len(error) > 80 else ''}")

            corrections = analysis.get('suggested_corrections', [])
            if corrections:
                print(f"\n  Suggested Corrections ({len(corrections)}):")
                for j, corr in enumerate(corrections[:2], 1):  # Show first 2
                    print(f"    {j}. [{corr['priority']}] {corr['issue'][:80]}...")

    # Save results
    output_dir = Path("reaction_agent/results")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "problematic_reactions_llm_assisted.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump({
            'test_type': 'llm_assisted_mapping_on_failed_rxnmapper',
            'num_reactions': len(reactions),
            'results': all_results
        }, f, indent=2, ensure_ascii=False)

    print(f"\n\n{'='*80}")
    print(f"✓ Results saved to: {output_file}")
    print(f"{'='*80}\n")


if __name__ == "__main__":
    test_problematic_reactions()
