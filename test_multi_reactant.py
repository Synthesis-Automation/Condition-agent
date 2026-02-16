#!/usr/bin/env python
"""
Test the fixed mapping failure system with a complex multi-reactant reaction.
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

# Multi-reactant reaction (4 reactants -> 1 product)
rxn_smiles = "O=[N+]([O-])c1cc(Cl)ccc1Cl.Sc1ccc(Br)cc1.O=C(Cl)CCCl.CN1CCNCC1>>CN1CCN(CCC(=O)Nc2cc(Cl)ccc2Sc2ccc(Br)cc2)CC1"

print("=" * 80)
print("Testing Multi-Reactant Reaction with Mapping Failure Fix")
print("=" * 80)
print(f"\nReaction SMILES:")
print(f"  {rxn_smiles}")
print("\nReactants:")
print("  1. O=[N+]([O-])c1cc(Cl)ccc1Cl  (3,5-dichloronitrobenzene)")
print("  2. Sc1ccc(Br)cc1               (4-bromothiophenol)")
print("  3. O=C(Cl)CCCl                 (chloroacetyl chloride)")
print("  4. CN1CCNCC1                   (N-methylpiperazine)")
print("\n" + "-" * 80)

# Initialize client and analyzer
print("\nInitializing analyzer with GPT-4o (auto mode will switch if needed)...")
client = LLMClient(provider="openai", model="gpt-4o", timeout=300)
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

# Run analysis with auto mode
print("\nRunning analysis with mode='auto'...")
print("-" * 80)

try:
    result = analyzer.analyze(rxn_smiles, mode="auto")

    print("\n" + "=" * 80)
    print("ANALYSIS RESULTS")
    print("=" * 80)

    # Mapping quality
    mapping_qc = result['tool_facts']['mapping_qc']
    bond_changes = result['tool_facts']['bond_changes']

    print(f"\n📊 DETERMINISTIC ANALYSIS:")
    print(f"  Mapping confidence: {mapping_qc.get('confidence', 0):.3f}")
    print(f"  Mapping OK: {mapping_qc.get('ok', False)}")
    print(f"  Bond changes detected: {len(bond_changes)}")

    if len(bond_changes) > 0:
        print(f"\n  Bond change summary:")
        for bc in bond_changes[:5]:
            print(f"    - {bc['id']}: {bc['change']} bond between :{bc['a1']} and :{bc['a2']} ({bc['bond']})")
        if len(bond_changes) > 5:
            print(f"    ... and {len(bond_changes) - 5} more")

    # Model selection
    metadata = result['metadata']
    print(f"\n🤖 MODEL SELECTION:")
    print(f"  Model used: {metadata.get('model_selected', metadata.get('model'))}")
    print(f"  Mode: {metadata.get('mode', 'N/A')}")
    print(f"  Reasoning effort: {metadata.get('reasoning_effort', 'None')}")

    # Check if direct analysis was used
    interpretation = result.get('interpretation', {})
    warnings = interpretation.get('warnings', [])
    if 'mapping_failed_used_direct_analysis' in warnings:
        print(f"\n⚠️  DIRECT SMILES ANALYSIS MODE WAS USED")

    # Interpretation
    print(f"\n🔬 LLM INTERPRETATION:")

    if 'error' in interpretation:
        print(f"  ❌ Error: {interpretation['error']}")
        print(f"\n  Raw response (first 500 chars):")
        print(f"  {interpretation.get('raw_response', 'N/A')[:500]}")
    else:
        print(f"  Reaction class: {interpretation.get('overall_class', 'N/A')}")
        print(f"  Tags: {', '.join(interpretation.get('tags', []))}")
        print(f"  Confidence: {interpretation.get('confidence', 0):.2f}")

        # Roles
        roles = interpretation.get('roles', {})
        if roles:
            print(f"\n  Molecular roles:")
            for role, value in roles.items():
                print(f"    - {role}: {value}")

        # Mechanism summary
        mechanism = interpretation.get('mechanism_summary', [])
        if mechanism:
            print(f"\n  Mechanism summary ({len(mechanism)} steps):")
            for i, step in enumerate(mechanism, 1):
                print(f"    {i}. {step}")

        # Events
        events = interpretation.get('events', [])
        if events:
            print(f"\n  Mechanistic events ({len(events)}):")
            for event in events:
                print(f"    - {event['event_id']}: {event['event_type']} (confidence: {event.get('confidence', 0):.2f})")
                print(f"      {event.get('short_rationale', 'N/A')}")

        # Warnings
        if warnings:
            print(f"\n  ⚠️  Warnings:")
            for warning in warnings:
                print(f"    - {warning}")

    # Performance metrics
    print(f"\n⏱️  PERFORMANCE:")
    print(f"  Total tokens: {metadata['total_tokens']}")
    print(f"  Latency: {metadata['latency_ms']:.0f} ms ({metadata['latency_ms']/1000:.1f}s)")

    print("\n" + "=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)

except Exception as e:
    print(f"\n❌ Analysis failed with error:")
    print(f"  {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
