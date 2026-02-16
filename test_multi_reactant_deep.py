#!/usr/bin/env python
"""
Test multi-reactant reaction with DEEP mode (force GPT-5.2).
Compare with auto mode results.
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
print("Testing Multi-Reactant Reaction: AUTO vs DEEP Mode Comparison")
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
print("\nInitializing analyzer...")
client = LLMClient(provider="openai", model="gpt-4o", timeout=300)
analyzer = ReactionSMILESAnalyzer(client, max_tokens=8000)

# Test 1: AUTO mode (already know it uses gpt-4o)
print("\n" + "=" * 80)
print("TEST 1: AUTO MODE")
print("=" * 80)

result_auto = analyzer.analyze(rxn_smiles, mode="auto")

print(f"\n📊 Model used: {result_auto['metadata'].get('model_selected', 'N/A')}")
print(f"   Reasoning: {result_auto['metadata'].get('reasoning_effort', 'None')}")
print(f"   Mapping confidence: {result_auto['tool_facts']['mapping_qc'].get('confidence', 0):.3f}")

interp_auto = result_auto.get('interpretation', {})
print(f"\n🔬 Results:")
print(f"   Class: {interp_auto.get('overall_class', 'N/A')}")
print(f"   Confidence: {interp_auto.get('confidence', 0):.2f}")
print(f"   Mechanism steps: {len(interp_auto.get('mechanism_summary', []))}")

# Test 2: DEEP mode (force GPT-5.2)
print("\n" + "=" * 80)
print("TEST 2: DEEP MODE (Force GPT-5.2)")
print("=" * 80)

result_deep = analyzer.analyze(rxn_smiles, mode="deep")

print(f"\n📊 Model used: {result_deep['metadata'].get('model_selected', 'N/A')}")
print(f"   Reasoning: {result_deep['metadata'].get('reasoning_effort', 'None')}")
print(f"   Mapping confidence: {result_deep['tool_facts']['mapping_qc'].get('confidence', 0):.3f}")

interp_deep = result_deep.get('interpretation', {})

if 'error' in interp_deep:
    print(f"\n❌ ERROR: {interp_deep['error']}")
    print(f"\nRaw response (first 1000 chars):")
    print(interp_deep.get('raw_response', 'N/A')[:1000])
else:
    print(f"\n🔬 Results:")
    print(f"   Class: {interp_deep.get('overall_class', 'N/A')}")
    print(f"   Tags: {', '.join(interp_deep.get('tags', []))}")
    print(f"   Confidence: {interp_deep.get('confidence', 0):.2f}")

    # Roles
    roles = interp_deep.get('roles', {})
    if roles:
        print(f"\n   Molecular roles:")
        for role, value in roles.items():
            print(f"     - {role}: {value[:80]}...")

    # Mechanism summary
    mechanism = interp_deep.get('mechanism_summary', [])
    print(f"\n   Mechanism summary ({len(mechanism)} steps):")
    for i, step in enumerate(mechanism, 1):
        print(f"     {i}. {step}")

    # Events
    events = interp_deep.get('events', [])
    if events:
        print(f"\n   Mechanistic events ({len(events)}):")
        for event in events:
            print(f"     - {event['event_id']}: {event['event_type']} (conf: {event.get('confidence', 0):.2f})")
            rationale = event.get('short_rationale', 'N/A')
            print(f"       {rationale[:100]}...")

    # Reaction centers (if present in deep mode)
    if 'reaction_centers' in interp_deep:
        centers = interp_deep['reaction_centers']
        print(f"\n   Reaction centers ({len(centers)}):")
        for center in centers:
            print(f"     - Center {center.get('center_id', '?')}: {center.get('description', 'N/A')[:80]}...")

    # Warnings
    warnings = interp_deep.get('warnings', [])
    if warnings:
        print(f"\n   ⚠️  Warnings:")
        for warning in warnings:
            print(f"     - {warning}")

# Comparison
print("\n" + "=" * 80)
print("COMPARISON: AUTO vs DEEP")
print("=" * 80)

print(f"\n{'Metric':<30} {'AUTO (gpt-4o)':<25} {'DEEP (gpt-5.2)':<25}")
print("-" * 80)

auto_class = interp_auto.get('overall_class', 'N/A')
deep_class = interp_deep.get('overall_class', 'N/A')
auto_conf = f"{interp_auto.get('confidence', 0):.2f}"
deep_conf = f"{interp_deep.get('confidence', 0):.2f}"
auto_steps = len(interp_auto.get('mechanism_summary', []))
deep_steps = len(interp_deep.get('mechanism_summary', []))
auto_events = len(interp_auto.get('events', []))
deep_events = len(interp_deep.get('events', []))
auto_time = f"{result_auto['metadata']['latency_ms']/1000:.1f}s"
deep_time = f"{result_deep['metadata']['latency_ms']/1000:.1f}s"
auto_tokens = result_auto['metadata']['total_tokens']
deep_tokens = result_deep['metadata']['total_tokens']

print(f"{'Reaction Class':<30} {auto_class:<25} {deep_class:<25}")
print(f"{'Confidence':<30} {auto_conf:<25} {deep_conf:<25}")
print(f"{'Mechanism Steps':<30} {auto_steps:<25} {deep_steps:<25}")
print(f"{'Events Identified':<30} {auto_events:<25} {deep_events:<25}")
print(f"{'Time (seconds)':<30} {auto_time:<25} {deep_time:<25}")
print(f"{'Total Tokens':<30} {auto_tokens:<25} {deep_tokens:<25}")

# Cost estimate
cost_auto = result_auto['metadata']['total_tokens'] * 0.0000025  # gpt-4o rough estimate
cost_deep = result_deep['metadata']['total_tokens'] * 0.000015   # gpt-5.2 rough estimate
print(f"{'Est. Cost':<30} ${cost_auto:<24.4f} ${cost_deep:<24.4f}")

print("\n" + "=" * 80)
print("ANALYSIS")
print("=" * 80)

# Check if deep mode identified more details
auto_mech = set(interp_auto.get('mechanism_summary', []))
deep_mech = set(interp_deep.get('mechanism_summary', []))

print(f"\nDoes DEEP mode identify more mechanistic details?")
print(f"  AUTO identified: {len(interp_auto.get('mechanism_summary', []))} steps")
print(f"  DEEP identified: {len(interp_deep.get('mechanism_summary', []))} steps")

if len(interp_deep.get('mechanism_summary', [])) > len(interp_auto.get('mechanism_summary', [])):
    print(f"  ✅ YES - DEEP mode found {len(interp_deep.get('mechanism_summary', [])) - len(interp_auto.get('mechanism_summary', []))} additional steps")
else:
    print(f"  ⚠️  NO - Similar depth")

# Check reaction class differences
if interp_auto.get('overall_class') != interp_deep.get('overall_class'):
    print(f"\n⚠️  REACTION CLASS DIFFERS:")
    print(f"  AUTO: {interp_auto.get('overall_class')}")
    print(f"  DEEP: {interp_deep.get('overall_class')}")

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)
