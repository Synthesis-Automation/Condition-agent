"""
Run LLM-assisted mapping on the Suzuki-Miyaura coupling reaction.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from reaction_agent.scripts.llm_assisted_mapping import hybrid_mapping_workflow
import json


def analyze_suzuki_with_llm():
    """Analyze Suzuki-Miyaura coupling with LLM assistance."""

    rxn_smiles = "CC1(C)CC=C(B2OC(C)(C)C(C)(C)O2)CC1.CS(=O)(=O)N1CCN(C(=O)c2cnc3ccc(F)cc3c2Br)CC1>>CC1(C)CC=C(c2c(C(=O)N3CCN(S(C)(=O)=O)CC3)cnc3ccc(F)cc23)CC1.O=C(O)C(F)(F)F"

    print("=" * 100)
    print("LLM-Assisted Mapping: Suzuki-Miyaura Coupling")
    print("=" * 100)
    print()
    print("This reaction showed low rxnmapper confidence (0.281) and disagreement between")
    print("all three mapping methods. Let's use LLM reasoning to understand why...")
    print()
    print("=" * 100)
    print()

    # Run hybrid mapping workflow with low confidence threshold to trigger LLM analysis
    # Note: rxnmapper confidence is 0.281, which is < 0.4, so LLM assistance will be triggered
    result = hybrid_mapping_workflow(
        rxn_smiles,
        confidence_threshold=0.6  # Triggers LLM for confidence < 0.6
    )

    print("\n" + "=" * 100)
    print("RESULTS")
    print("=" * 100)
    print()

    # Show deterministic results
    print("1. Deterministic Analysis (rxnmapper):")
    print("-" * 100)
    if result.get('deterministic_result'):
        det = result['deterministic_result']
        print(f"  Confidence: {det['tool_facts']['mapping_qc'].get('confidence', 'N/A'):.3f}")
        print(f"  Bond changes: {len(det['tool_facts'].get('bond_changes', []))}")
        print(f"  Reaction center atoms: {len(det['tool_facts'].get('reaction_center_atoms', []))}")
    print()

    # Show LLM analysis
    print("2. LLM Assistance:")
    print("-" * 100)
    if result.get('llm_assistance'):
        llm_result = result['llm_assistance']

        print()
        print("✓ LLM assistance triggered due to low rxnmapper confidence (<0.4)")
        print()

        # Check if parsed successfully
        if llm_result.get('parsed_successfully'):
            analysis = llm_result.get('llm_analysis', {})

            # Reaction type
            if analysis.get('reaction_analysis'):
                rxn_analysis = analysis['reaction_analysis']
                print("Reaction Analysis:")
                print(f"  Type: {rxn_analysis.get('type', 'Unknown')}")
                print(f"  Complexity: {rxn_analysis.get('complexity', 'Unknown')}")
                print()

                # Mechanism stages
                if rxn_analysis.get('stages'):
                    print("Mechanism Stages:")
                    for i, stage in enumerate(rxn_analysis['stages'], 1):
                        print(f"  Stage {i}: {stage.get('name', 'Unknown')}")
                        if stage.get('mechanism'):
                            print(f"    Mechanism: {stage['mechanism']}")
                    print()

            # Mapping assessment
            if analysis.get('mapping_assessment'):
                mapping_assess = analysis['mapping_assessment']
                print("Mapping Assessment:")
                print(f"  Current mapping correct: {mapping_assess.get('current_mapping_correct', 'Unknown')}")
                print(f"  Confidence in current mapping: {mapping_assess.get('confidence_in_current', 'N/A')}")
                if mapping_assess.get('major_errors'):
                    print(f"  Major errors found: {len(mapping_assess['major_errors'])}")
                    for error in mapping_assess['major_errors']:
                        print(f"    - {error}")
                print()

            # Suggested corrections
            if analysis.get('suggested_corrections'):
                corrections = analysis['suggested_corrections']
                print(f"Suggested Corrections: {len(corrections)}")
                for i, correction in enumerate(corrections, 1):
                    priority = correction.get('priority', 'unknown').upper()
                    issue = correction.get('issue', 'Unknown')
                    print(f"  {i}. [{priority}] {issue}")
                print()

            # Recommended action
            if analysis.get('recommended_action'):
                print("Recommended Action:")
                print(f"  {analysis['recommended_action']}")
                print()

            # Reasoning
            if analysis.get('reasoning'):
                print("LLM Reasoning:")
                # Print first 500 chars of reasoning
                reasoning = analysis['reasoning']
                if len(reasoning) > 500:
                    print(f"  {reasoning[:500]}...")
                else:
                    print(f"  {reasoning}")
                print()
        else:
            print("  ✗ LLM parsing failed")
            print(f"  Error: {llm_result.get('error', 'Unknown error')}")
            print()

    elif result.get('llm_validation'):
        print()
        print("✓ LLM validation triggered (confidence between 0.4-0.6)")
        print()
        validation = result['llm_validation']
        if validation.get('validation'):
            val = validation['validation']
            print(f"  Mapping valid: {val.get('mapping_valid', 'Unknown')}")
            print(f"  Reaction type: {val.get('reaction_type', 'Unknown')}")
            print(f"  Recommendation: {val.get('recommendation', 'Unknown')}")
            if val.get('issues_found'):
                print(f"  Issues found: {len(val['issues_found'])}")
                for issue in val['issues_found']:
                    print(f"    - {issue}")
            print()
    else:
        print("  ✗ LLM analysis not triggered (confidence above threshold)")
        print()

    # Show final confidence
    print("3. Final Assessment:")
    print("-" * 100)
    print(f"  Final confidence: {result.get('final_confidence', 'N/A'):.3f}")
    print(f"  Mapping method: {result.get('method', 'Unknown')}")
    print()

    # Save detailed results
    output_file = "results/suzuki_llm_analysis.json"
    print(f"Saving detailed results to: {output_file}")

    import os
    os.makedirs("results", exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2, default=str)

    print(f"✓ Saved to {output_file}")
    print()

    print("=" * 100)
    print()

    return result


if __name__ == "__main__":
    try:
        result = analyze_suzuki_with_llm()

        print()
        print("=" * 100)
        print("SUMMARY")
        print("=" * 100)
        print()
        print("The LLM-assisted mapping provides:")
        print("  ✓ Mechanistic understanding of the Suzuki coupling")
        print("  ✓ Identification of specific mapping errors")
        print("  ✓ Suggested corrections for manual review")
        print("  ✓ Explanation of why rxnmapper struggled")
        print()
        print("Use this analysis to:")
        print("  1. Understand the reaction mechanism")
        print("  2. Manually correct the atom mapping")
        print("  3. Validate against chemistry knowledge")
        print()

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
