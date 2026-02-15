#!/usr/bin/env python
"""
Demo script for Reaction SMILES Analysis Agent (POC)

Demonstrates the full pipeline:
1. Deterministic analysis (cleaning, mapping, bond changes)
2. LLM interpretation with structured output

Usage:
    python reaction_agent/scripts/demo.py

Requirements:
    - OPENAI_API_KEY environment variable set
    - rdkit installed
    - rxnmapper installed (optional, for mapping)
"""

import os
import json
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer


def print_section(title: str):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def print_result(result: dict):
    """Pretty-print analysis result."""

    # Input section
    print_section("INPUT")
    print(f"Raw:   {result['input']['rxn_smiles_raw']}")
    print(f"Clean: {result['input']['rxn_smiles_clean']}")
    if result['input']['spectators']:
        print(f"Spectators: {', '.join(result['input']['spectators'])}")
    if result['input']['parse_warnings']:
        print(f"⚠️  Warnings: {', '.join(result['input']['parse_warnings'])}")

    # Tool facts section
    print_section("TOOL FACTS (Deterministic)")
    tool_facts = result['tool_facts']

    if tool_facts['mapped_rxn_smiles']:
        print(f"Mapped: {tool_facts['mapped_rxn_smiles']}")

    mapping_qc = tool_facts['mapping_qc']
    qc_status = "✓ OK" if mapping_qc.get('ok', False) else "✗ Failed"
    print(f"Mapping QC: {qc_status}")
    if 'confidence' in mapping_qc:
        print(f"  Confidence: {mapping_qc['confidence']:.2f}")
    if mapping_qc.get('notes'):
        for note in mapping_qc['notes']:
            print(f"  Note: {note}")

    if tool_facts['bond_changes']:
        print(f"\nBond Changes ({len(tool_facts['bond_changes'])}):")
        for bc in tool_facts['bond_changes']:
            print(f"  {bc['id']}: {bc['change']} bond between :{bc['a1']} and :{bc['a2']} ({bc['bond']})")

    if tool_facts['reaction_center_atoms']:
        print(f"\nReaction Center Atoms: {tool_facts['reaction_center_atoms']}")

    # Interpretation section
    print_section("LLM INTERPRETATION")
    interp = result['interpretation']

    if 'error' in interp:
        print(f"❌ Error: {interp['error']}")
        if 'raw_response' in interp:
            print(f"\nRaw response:\n{interp['raw_response'][:500]}...")
        return

    print(f"Reaction Class: {interp.get('overall_class', 'N/A')}")
    if interp.get('tags'):
        print(f"Tags: {', '.join(interp['tags'])}")

    print(f"\nConfidence: {interp.get('confidence', 0.0):.2f}")

    if interp.get('roles'):
        print("\nRoles:")
        for role, value in interp['roles'].items():
            print(f"  {role}: {value}")

    if interp.get('events'):
        print(f"\nMechanistic Events ({len(interp['events'])}):")
        for event in interp['events']:
            print(f"\n  {event['event_id']}: {event['event_type']}")
            print(f"    Bond changes: {', '.join(event.get('bond_change_refs', []))}")
            print(f"    Rationale: {event.get('short_rationale', 'N/A')}")
            print(f"    Confidence: {event.get('confidence', 0.0):.2f}")

    if interp.get('mechanism_summary'):
        print("\nMechanism Summary:")
        for i, step in enumerate(interp['mechanism_summary'], 1):
            print(f"  {i}. {step}")

    if interp.get('warnings'):
        print("\n⚠️  Warnings:")
        for warning in interp['warnings']:
            print(f"  - {warning}")

    # Metadata section
    print_section("METADATA")
    meta = result['metadata']
    print(f"Model: {meta['model']} ({meta['provider']})")
    print(f"Tokens: {meta['total_tokens']} (prompt: {meta['prompt_tokens']}, completion: {meta['completion_tokens']})")
    print(f"Latency: {meta['latency_ms']:.0f} ms")


def main():
    """Run test cases."""

    # Check for API key
    if not os.getenv("OPENAI_API_KEY"):
        print("❌ Error: OPENAI_API_KEY environment variable not set")
        print("\nPlease set it:")
        print("  export OPENAI_API_KEY='your-key-here'")
        sys.exit(1)

    print_section("Reaction SMILES Analysis Agent - POC Demo")
    print("\nInitializing LLM client...")

    # Initialize client
    try:
        client = LLMClient(
            provider="openai",
            model="gpt-4o-mini",  # Use mini for faster/cheaper testing
            temperature=0.0,
            max_tokens=2000
        )
        print(f"✓ Client initialized: {client.model}")
    except Exception as e:
        print(f"❌ Failed to initialize client: {e}")
        sys.exit(1)

    # Initialize analyzer
    analyzer = ReactionSMILESAnalyzer(
        client=client,
        drop_spectators=True,
        temperature=0.0,
        max_tokens=2000
    )
    print("✓ Analyzer initialized")

    # Test cases
    test_reactions = [
        {
            "name": "SNAr (Nucleophilic Aromatic Substitution)",
            "smiles": "CN(C)C.Clc1cc(Cl)nc(Cl)c1>>CN(C)c1cc(Cl)nc(Cl)c1",
            "description": "Dimethylamine displacing chloride on trichloropyrimidine"
        },
        {
            "name": "Simple Alkylation",
            "smiles": "CCBr.N>>CCNC",
            "description": "Ethyl bromide + amine"
        },
        {
            "name": "Buchwald-Hartwig CN Coupling",
            "smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
            "description": "Bromobenzene + aniline coupling"
        },
    ]

    # Run test cases
    for i, test in enumerate(test_reactions, 1):
        print_section(f"TEST CASE {i}: {test['name']}")
        print(f"Description: {test['description']}")
        print(f"SMILES: {test['smiles']}")

        try:
            result = analyzer.analyze(test['smiles'])

            if 'error' in result:
                print(f"\n❌ Analysis failed: {result['error']}")
                continue

            print_result(result)

        except Exception as e:
            print(f"\n❌ Exception during analysis: {e}")
            import traceback
            traceback.print_exc()

    print_section("DEMO COMPLETE")
    print("\n✓ All tests completed")


if __name__ == "__main__":
    main()
