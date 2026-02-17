"""
Test script for local environment atom mapping.

Demonstrates the new deterministic mapping method for functional group transformations.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools._atom_mapping import map_by_local_environment, analyze_bond_changes_hybrid


def test_simple_functional_group_transformations():
    """Test local environment mapping on simple functional group transformations."""

    print("=" * 80)
    print("Testing Local Environment Mapping on Functional Group Transformations")
    print("=" * 80)
    print()

    # Test cases: functional group transformations where rxnmapper might struggle
    test_reactions = [
        {
            "name": "Alcohol to Acetate (Protection)",
            "smiles": "CC(O)c1ccccc1>>CC(OC(C)=O)c1ccccc1",
            "description": "Simple esterification"
        },
        {
            "name": "Bromide to Amine (SN2)",
            "smiles": "CCCBr.N>>CCCN",
            "description": "Nucleophilic substitution"
        },
        {
            "name": "Aldehyde to Alcohol (Reduction)",
            "smiles": "CC(=O)c1ccccc1>>CC(O)c1ccccc1",
            "description": "Carbonyl reduction"
        },
        {
            "name": "Amine to Amide",
            "smiles": "CCN.CC(=O)O>>CCNC(=O)C",
            "description": "Amide bond formation"
        },
    ]

    for i, test in enumerate(test_reactions, 1):
        print(f"\n{'='*80}")
        print(f"Test {i}: {test['name']}")
        print(f"Description: {test['description']}")
        print(f"{'='*80}")
        print(f"Reaction: {test['smiles']}")
        print()

        # Test local environment mapping alone
        print("--- Local Environment Mapping (Standalone) ---")
        result = map_by_local_environment(test['smiles'])

        if result['success']:
            print(f"✓ Success!")
            print(f"  Confidence: {result['confidence']:.3f}")
            print(f"  Mapped SMILES: {result['mapped_smiles']}")
            print(f"  Broken bonds: {len(result['broken_bonds'])}")
            print(f"  Formed bonds: {len(result['formed_bonds'])}")
            print(f"  Changed atoms: {len(result['changed_atoms'])}")
            print(f"  Spectator atoms: {len(result['spectator_atoms'])}")
            print(f"  Interpretation: {result['interpretation']}")

            if result.get('mapping_stats'):
                stats = result['mapping_stats']
                print(f"  Mapping coverage: {stats['coverage']:.2%}")
                print(f"  Unmapped reactants: {stats['unmapped_reactants']}")
                print(f"  Unmapped products: {stats['unmapped_products']}")
        else:
            print(f"✗ Failed: {result.get('error', 'Unknown error')}")

        print()

        # Test hybrid approach (comparing with rxnmapper if available)
        print("--- Hybrid Approach (All Methods) ---")
        hybrid_result = analyze_bond_changes_hybrid(
            test['smiles'],
            use_rxnmapper=True,
            use_local_env=True,
            use_mcs=True
        )

        if hybrid_result['success']:
            print(f"✓ Success!")
            print(f"  Best method: {hybrid_result['method']}")
            print(f"  Combined confidence: {hybrid_result['combined_confidence']:.3f}")
            print(f"  Validation: {hybrid_result['validation']}")

            # Show agreement between methods
            print("\n  Method Agreements:")
            for key, value in hybrid_result['agreement'].items():
                if value is not None:
                    status = "✓ Agree" if value else "✗ Disagree"
                    print(f"    {key}: {status}")

            # Show which methods succeeded
            print("\n  Methods that succeeded:")
            if hybrid_result['rxnmapper_result'] and hybrid_result['rxnmapper_result'].get('success'):
                rxn_conf = hybrid_result['rxnmapper_result'].get('mapping_confidence', 'N/A')
                print(f"    - RXNMapper (confidence: {rxn_conf})")
            if hybrid_result['local_env_result'] and hybrid_result['local_env_result'].get('success'):
                local_conf = hybrid_result['local_env_result'].get('confidence', 'N/A')
                print(f"    - Local environment (confidence: {local_conf})")
            if hybrid_result['mcs_result'] and hybrid_result['mcs_result'].get('success'):
                mcs_conf = hybrid_result['mcs_result'].get('confidence', 'N/A')
                print(f"    - MCS (confidence: {mcs_conf})")
        else:
            print(f"✗ Failed: {hybrid_result.get('error', 'Unknown error')}")

        print()

    print("=" * 80)
    print("Testing Complete!")
    print("=" * 80)


def test_when_rxnmapper_fails():
    """Test cases where rxnmapper might fail but local environment should work."""

    print("\n" + "=" * 80)
    print("Testing Cases Where RXNMapper Struggles")
    print("=" * 80)
    print()

    # These are simpler reactions where deterministic methods should excel
    simple_reactions = [
        "CCO>>CC=O",  # Alcohol oxidation
        "c1ccccc1Br.N>>c1ccccc1N",  # Simple aromatic substitution
        "CC(=O)O.CCO>>CC(=O)OCC",  # Esterification
    ]

    for rxn in simple_reactions:
        print(f"\nReaction: {rxn}")
        print("-" * 80)

        result = analyze_bond_changes_hybrid(
            rxn,
            use_rxnmapper=True,
            use_local_env=True,
            use_mcs=False
        )

        if result['success']:
            method = result['method']
            conf = result['combined_confidence']

            print(f"Best method: {method} (confidence: {conf:.3f})")

            # Highlight when local_env succeeds where rxnmapper struggles
            rxn_ok = result['rxnmapper_result'] and result['rxnmapper_result'].get('success')
            local_ok = result['local_env_result'] and result['local_env_result'].get('success')

            if local_ok and not rxn_ok:
                print("⭐ Local environment succeeded where RXNMapper failed!")
            elif local_ok and rxn_ok:
                rxn_conf = result['rxnmapper_result'].get('mapping_confidence', 0)
                local_conf = result['local_env_result'].get('confidence', 0)
                if local_conf > rxn_conf:
                    print(f"⭐ Local environment more confident ({local_conf:.3f} vs {rxn_conf:.3f})")
        else:
            print(f"Failed: {result.get('error')}")


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("Local Environment Atom Mapping Test Suite")
    print("=" * 80)

    try:
        test_simple_functional_group_transformations()
        test_when_rxnmapper_fails()

        print("\n" + "=" * 80)
        print("Summary")
        print("=" * 80)
        print("""
The local environment mapper provides:

✓ Deterministic atom mapping (no ML required)
✓ High confidence (0.7-0.9) for functional group transformations
✓ Complementary to rxnmapper
✓ Fast and reliable
✓ Works well when:
  - Most of the molecule is unchanged (spectators)
  - Clear functional group transformation
  - Simple substitution or addition/elimination

Use in hybrid workflow for best results:
  1. Try rxnmapper first (high accuracy for all types)
  2. If fails or low confidence, local environment kicks in
  3. MCS as final fallback

Example usage:
    from chemtools._atom_mapping import analyze_bond_changes_hybrid

    result = analyze_bond_changes_hybrid(
        "CCO>>CC=O",  # Alcohol oxidation
        use_rxnmapper=True,
        use_local_env=True,
        use_mcs=True
    )

    print(f"Method: {result['method']}")
    print(f"Confidence: {result['combined_confidence']:.3f}")
        """)

    except ImportError as e:
        print(f"\n❌ Import error: {e}")
        print("Make sure RDKit is installed: pip install rdkit")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
