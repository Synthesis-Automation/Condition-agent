"""
Test tandem Suzuki coupling + acetal hydrolysis reaction.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools._atom_mapping import analyze_bond_changes_hybrid, map_by_local_environment
from reaction_agent.scripts.llm_assisted_mapping import hybrid_mapping_workflow


def test_tandem_reaction():
    """Test Suzuki coupling with acetal hydrolysis."""

    rxn_smiles = "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O"

    print("=" * 100)
    print("Testing Tandem Reaction: Suzuki Coupling + Acetal Hydrolysis")
    print("=" * 100)
    print()
    print("Expected transformations:")
    print("  1. Suzuki coupling: Ar-Br + R-B(OR)2 -> Ar-R")
    print("  2. Acetal hydrolysis: -CH(OCH3)2 -> -CHO")
    print()
    print("SMILES:")
    print(rxn_smiles)
    print()

    # Test 1: Local environment mapping
    print("=" * 100)
    print("1. Local Environment Mapping (Standalone)")
    print("=" * 100)
    print()

    local_result = map_by_local_environment(rxn_smiles)

    if local_result['success']:
        print("✓ SUCCESS!")
        print()
        print(f"Confidence: {local_result['confidence']:.3f}")
        print(f"Interpretation: {local_result['interpretation']}")
        print()
        print("Bond Changes:")
        print(f"  Broken bonds: {len(local_result['broken_bonds'])}")
        print(f"  Formed bonds: {len(local_result['formed_bonds'])}")
        print()
        print("Mapping Statistics:")
        stats = local_result['mapping_stats']
        print(f"  Coverage: {stats['coverage']:.1%}")
        print(f"  Total reactant atoms: {stats['total_reactant_atoms']}")
        print(f"  Total product atoms: {stats['total_product_atoms']}")
        print(f"  Mapped atoms: {stats['mapped_atoms']}")
        print(f"  Unmapped reactants: {stats['unmapped_reactants']}")
        print(f"  Unmapped products: {stats['unmapped_products']}")
    else:
        print(f"✗ FAILED: {local_result.get('error', 'Unknown error')}")

    print()
    print()

    # Test 2: Hybrid approach
    print("=" * 100)
    print("2. Hybrid Approach (All Methods)")
    print("=" * 100)
    print()

    hybrid_result = analyze_bond_changes_hybrid(
        rxn_smiles,
        use_rxnmapper=True,
        use_local_env=True,
        use_mcs=True
    )

    if hybrid_result['success']:
        print("✓ SUCCESS!")
        print()
        print(f"Best Method: {hybrid_result['method']}")
        print(f"Combined Confidence: {hybrid_result['combined_confidence']:.3f}")
        print(f"Validation: {hybrid_result['validation']}")
        print()

        # Show which methods succeeded
        print("Methods Performance:")

        if hybrid_result['rxnmapper_result'] and hybrid_result['rxnmapper_result'].get('success'):
            rxn_conf = hybrid_result['rxnmapper_result'].get('mapping_confidence', 'N/A')
            rxn_broken = len(hybrid_result['rxnmapper_result'].get('broken_bonds', []))
            rxn_formed = len(hybrid_result['rxnmapper_result'].get('formed_bonds', []))
            print(f"  ✓ RXNMapper: conf={rxn_conf:.3f}, {rxn_broken} broken, {rxn_formed} formed")
        else:
            print(f"  ✗ RXNMapper: failed")

        if hybrid_result['local_env_result'] and hybrid_result['local_env_result'].get('success'):
            local_conf = hybrid_result['local_env_result'].get('confidence', 'N/A')
            local_broken = len(hybrid_result['local_env_result'].get('broken_bonds', []))
            local_formed = len(hybrid_result['local_env_result'].get('formed_bonds', []))
            print(f"  ✓ Local Env: conf={local_conf:.3f}, {local_broken} broken, {local_formed} formed")
        else:
            print(f"  ✗ Local Env: failed")

        if hybrid_result['mcs_result'] and hybrid_result['mcs_result'].get('success'):
            mcs_conf = hybrid_result['mcs_result'].get('confidence', 'N/A')
            mcs_broken = hybrid_result['mcs_result'].get('likely_broken_bonds', 0)
            mcs_formed = hybrid_result['mcs_result'].get('likely_formed_bonds', 0)
            print(f"  ✓ MCS: conf={mcs_conf:.3f}, {mcs_broken} broken, {mcs_formed} formed")
        else:
            print(f"  ✗ MCS: failed")

        print()

        # Show agreements
        print("Method Agreements:")
        agreements = hybrid_result['agreement']
        for key, value in agreements.items():
            if value is not None:
                status = "✓ Agree" if value else "✗ Disagree"
                print(f"  {key}: {status}")

        print()

        # Recommended result
        rec = hybrid_result['recommended_result']
        print("Recommended Result:")
        print(f"  Method: {rec.get('method', 'Unknown')}")
        print(f"  Broken bonds: {len(rec.get('broken_bonds', []))}")
        print(f"  Formed bonds: {len(rec.get('formed_bonds', []))}")
        print(f"  Changed atoms: {len(rec.get('changed_atoms', []))}")
    else:
        print(f"✗ FAILED: {hybrid_result.get('error', 'Unknown error')}")

    print()
    print()

    # Test 3: LLM assistance if needed
    print("=" * 100)
    print("3. LLM-Assisted Analysis (if needed)")
    print("=" * 100)
    print()

    # Check if we need LLM assistance
    rxn_conf = 1.0  # default high
    if hybrid_result.get('rxnmapper_result'):
        rxn_conf = hybrid_result['rxnmapper_result'].get('mapping_confidence', 1.0)

    all_disagree = (
        hybrid_result['agreement'].get('rxnmapper_vs_local_env') == False and
        hybrid_result['agreement'].get('rxnmapper_vs_mcs') == False and
        hybrid_result['agreement'].get('local_env_vs_mcs') == False
    )

    if rxn_conf < 0.6 or all_disagree:
        print(f"Triggering LLM analysis (rxnmapper conf: {rxn_conf:.3f}, all methods disagree: {all_disagree})")
        print()

        llm_result = hybrid_mapping_workflow(rxn_smiles, confidence_threshold=0.7)

        if llm_result.get('llm_assistance') or llm_result.get('llm_validation'):
            llm_data = llm_result.get('llm_assistance') or llm_result.get('llm_validation')

            if llm_data.get('llm_analysis'):
                analysis = llm_data['llm_analysis']

                print("LLM Analysis:")
                if analysis.get('reaction_analysis'):
                    rxn_analysis = analysis['reaction_analysis']
                    print(f"  Type: {rxn_analysis.get('type', 'Unknown')}")
                    print(f"  Complexity: {rxn_analysis.get('complexity', 'Unknown')}")

                    if rxn_analysis.get('stages'):
                        print(f"  Stages: {len(rxn_analysis['stages'])}")
                        for i, stage in enumerate(rxn_analysis['stages'], 1):
                            print(f"    {i}. {stage.get('name', 'Unknown')}")

                print()

                if analysis.get('mapping_assessment'):
                    assess = analysis['mapping_assessment']
                    print(f"  Current mapping correct: {assess.get('current_mapping_correct', 'Unknown')}")
                    if assess.get('major_errors'):
                        print(f"  Major errors: {len(assess['major_errors'])}")

                print()

                if analysis.get('recommended_action'):
                    print(f"  Recommendation: {analysis['recommended_action']}")

            elif llm_data.get('validation'):
                validation = llm_data['validation']
                print(f"  Mapping valid: {validation.get('mapping_valid', 'Unknown')}")
                print(f"  Recommendation: {validation.get('recommendation', 'Unknown')}")
        else:
            print("LLM analysis available, check detailed results.")
    else:
        print(f"✗ Not needed (rxnmapper conf: {rxn_conf:.3f}, methods agree)")

    print()
    print()

    # Summary
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print()

    print("Expected Chemistry:")
    print("  Step 1: Suzuki coupling (C-Br + C-B → C-C)")
    print("  Step 2: Acetal hydrolysis (CH(OMe)2 → CHO)")
    print()

    if hybrid_result['success']:
        rec_method = hybrid_result['method']
        rec_conf = hybrid_result['combined_confidence']

        print(f"Best Mapping: {rec_method} (confidence: {rec_conf:.3f})")

        # Check if it's a multi-stage reaction
        rec_broken = len(hybrid_result['recommended_result'].get('broken_bonds', []))
        rec_formed = len(hybrid_result['recommended_result'].get('formed_bonds', []))

        print(f"Bond changes: {rec_broken} broken, {rec_formed} formed")

        if rec_broken > 4 or rec_formed > 4:
            print()
            print("⚠ Complex transformation detected (>4 bond changes)")
            print("  This is expected for tandem reactions")
            print("  Consider validating with:")
            print("  - Manual inspection of key transformations")
            print("  - LLM-assisted mechanistic analysis")

        print()

        if all_disagree:
            print("⚠ All methods disagree - manual review recommended")
    else:
        print("✗ Mapping failed - needs manual attention")

    print()
    print("=" * 100)


if __name__ == "__main__":
    try:
        test_tandem_reaction()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
