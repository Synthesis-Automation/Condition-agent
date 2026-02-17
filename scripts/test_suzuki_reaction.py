"""
Test the Suzuki-Miyaura coupling reaction with local environment mapping.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools._atom_mapping import analyze_bond_changes_hybrid, map_by_local_environment


def test_suzuki_coupling():
    """Test Suzuki-Miyaura coupling with all mapping methods."""

    rxn_smiles = "CC1(C)CC=C(B2OC(C)(C)C(C)(C)O2)CC1.CS(=O)(=O)N1CCN(C(=O)c2cnc3ccc(F)cc3c2Br)CC1>>CC1(C)CC=C(c2c(C(=O)N3CCN(S(C)(=O)=O)CC3)cnc3ccc(F)cc23)CC1.O=C(O)C(F)(F)F"

    print("=" * 100)
    print("Testing Suzuki-Miyaura Coupling Reaction")
    print("=" * 100)
    print()
    print("Reaction Type: Cross-coupling (C-C bond formation)")
    print("Key transformation: Ar-Br + Ar-B(OR)2 -> Ar-Ar + Br- + B(OR)2OH")
    print()
    print("SMILES:")
    print(rxn_smiles)
    print()

    # Test local environment mapping standalone
    print("=" * 100)
    print("1. Local Environment Mapping (Standalone)")
    print("=" * 100)
    print()

    local_result = map_by_local_environment(rxn_smiles)

    if local_result['success']:
        print("✓ SUCCESS!")
        print()
        print(f"Confidence: {local_result['confidence']:.3f}")
        print(f"Method: {local_result['method']}")
        print()
        print(f"Mapped SMILES:")
        print(f"  {local_result['mapped_smiles']}")
        print()
        print("Bond Changes:")
        print(f"  Broken bonds: {len(local_result['broken_bonds'])}")
        for bond in local_result['broken_bonds']:
            print(f"    - {bond}")
        print(f"  Formed bonds: {len(local_result['formed_bonds'])}")
        for bond in local_result['formed_bonds']:
            print(f"    - {bond}")
        print()
        print(f"Reaction Center:")
        print(f"  Changed atoms: {len(local_result['changed_atoms'])} atoms")
        print(f"  Spectator atoms: {len(local_result['spectator_atoms'])} atoms")
        print()
        print("Mapping Statistics:")
        stats = local_result['mapping_stats']
        print(f"  Total reactant atoms: {stats['total_reactant_atoms']}")
        print(f"  Total product atoms: {stats['total_product_atoms']}")
        print(f"  Mapped atoms: {stats['mapped_atoms']}")
        print(f"  Coverage: {stats['coverage']:.1%}")
        print(f"  Unmapped reactants: {stats['unmapped_reactants']}")
        print(f"  Unmapped products: {stats['unmapped_products']}")
        print()
        print(f"Interpretation: {local_result['interpretation']}")
    else:
        print(f"✗ FAILED: {local_result.get('error', 'Unknown error')}")

    print()
    print()

    # Test hybrid approach (all methods)
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
        print("Methods That Succeeded:")
        methods_used = []

        if hybrid_result['rxnmapper_result'] and hybrid_result['rxnmapper_result'].get('success'):
            rxn_conf = hybrid_result['rxnmapper_result'].get('mapping_confidence', 'N/A')
            print(f"  ✓ RXNMapper (confidence: {rxn_conf:.3f})")
            methods_used.append('rxnmapper')
        else:
            print(f"  ✗ RXNMapper failed")

        if hybrid_result['local_env_result'] and hybrid_result['local_env_result'].get('success'):
            local_conf = hybrid_result['local_env_result'].get('confidence', 'N/A')
            print(f"  ✓ Local Environment (confidence: {local_conf:.3f})")
            methods_used.append('local_env')
        else:
            print(f"  ✗ Local Environment failed")

        if hybrid_result['mcs_result'] and hybrid_result['mcs_result'].get('success'):
            mcs_conf = hybrid_result['mcs_result'].get('confidence', 'N/A')
            print(f"  ✓ MCS (confidence: {mcs_conf:.3f})")
            methods_used.append('mcs')
        else:
            print(f"  ✗ MCS failed")

        print()

        # Show agreements
        print("Method Agreements (within 2 bonds tolerance):")
        for key, value in hybrid_result['agreement'].items():
            if value is not None:
                status = "✓ Agree" if value else "✗ Disagree"
                print(f"  {key}: {status}")

        print()

        # Show recommended result details
        rec = hybrid_result['recommended_result']
        print("Recommended Result Details:")
        print(f"  Broken bonds: {len(rec.get('broken_bonds', []))}")
        print(f"  Formed bonds: {len(rec.get('formed_bonds', []))}")
        print(f"  Changed atoms: {len(rec.get('changed_atoms', []))}")

        # Show specific bond changes if available
        if rec.get('broken_bonds'):
            print()
            print("  Broken bond details:")
            for bond in rec['broken_bonds'][:5]:  # Show first 5
                print(f"    - {bond}")

        if rec.get('formed_bonds'):
            print()
            print("  Formed bond details:")
            for bond in rec['formed_bonds'][:5]:  # Show first 5
                print(f"    - {bond}")

    else:
        print(f"✗ FAILED: {hybrid_result.get('error', 'Unknown error')}")

    print()
    print()

    # Analysis
    print("=" * 100)
    print("Analysis")
    print("=" * 100)
    print()

    print("Expected for Suzuki-Miyaura Coupling:")
    print("  • Broken bonds: 2-3 (C-Br, B-O bonds)")
    print("  • Formed bonds: 1-2 (new C-C bond)")
    print("  • Key transformation: Aryl bromide + Boronic ester -> Biaryl")
    print()

    if local_result['success'] and hybrid_result['success']:
        local_broken = len(local_result['broken_bonds'])
        local_formed = len(local_result['formed_bonds'])
        hybrid_broken = len(hybrid_result['recommended_result'].get('broken_bonds', []))
        hybrid_formed = len(hybrid_result['recommended_result'].get('formed_bonds', []))

        print("Comparison:")
        print(f"  Local Environment: {local_broken} broken, {local_formed} formed")
        print(f"  Hybrid (Best):     {hybrid_broken} broken, {hybrid_formed} formed")
        print()

        if local_broken == hybrid_broken and local_formed == hybrid_formed:
            print("✓ Methods agree on bond counts!")
        else:
            print("⚠ Methods disagree on bond counts - review recommended")

        print()

        # Highlight local environment performance
        if hybrid_result['method'] == 'local_env_only':
            print("⭐ SPECIAL: Local environment was the ONLY successful method!")
        elif 'local_env' in str(hybrid_result.get('validation', '')).lower():
            if '✓' in str(hybrid_result.get('validation', '')):
                print("✓ Local environment validated the primary method")

    print()
    print("=" * 100)


if __name__ == "__main__":
    try:
        test_suzuki_coupling()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
