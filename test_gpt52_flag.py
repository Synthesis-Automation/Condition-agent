#!/usr/bin/env python
"""Test that --model gpt-5.2 flag works correctly."""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Test simple reaction
rxn = "CCBr.CCN>>CCNCC"

print("="*80)
print("Testing --model gpt-5.2 flag")
print("="*80)
print(f"\nReaction: {rxn}")
print("\nExpected behavior:")
print("  - Should detect 'gpt-5.2' as reasoning model")
print("  - Should automatically use 'deep' mode")
print("  - Should set reasoning_effort='low'")
print("\n" + "-"*80)

from llmtools.clients import LLMClient
from reaction_agent.agent import ReactionSMILESAnalyzer

# Simulate what the CLI does
model = "gpt-5.2"
is_gpt52_model = model.startswith("gpt-5.2") or model.startswith("gpt-o3") or model.startswith("o3-")
mode = "auto"  # Default

# Auto-detect if user wants to force GPT-5.2
effective_mode = mode
if is_gpt52_model and mode == "auto":
    effective_mode = "deep"
    print(f"✅ Detected reasoning model '{model}' - using deep reasoning mode")

# Set reasoning effort
reasoning_effort = None
if effective_mode == "deep":
    reasoning_effort = "low"

print(f"\nConfiguration:")
print(f"  Model: {model}")
print(f"  Mode: {mode} -> {effective_mode}")
print(f"  Reasoning effort: {reasoning_effort}")

# Initialize client and analyzer
try:
    client = LLMClient(provider="openai", model=model, timeout=300)
    analyzer = ReactionSMILESAnalyzer(
        client=client,
        drop_spectators=True,
        temperature=0.0,
        max_tokens=8000,
        reasoning_effort=reasoning_effort
    )

    print(f"\n✅ CLIENT INITIALIZED")
    print(f"   Client model: {client.model}")
    print(f"   Analyzer reasoning_effort: {analyzer.reasoning_effort}")

    # Run analysis
    print(f"\n⏳ Running analysis with mode='{effective_mode}'...")
    result = analyzer.analyze(rxn, mode=effective_mode)

    print(f"\n📊 RESULTS:")
    print(f"   Model used: {result['metadata'].get('model_selected', 'N/A')}")
    print(f"   Reasoning effort: {result['metadata'].get('reasoning_effort', 'N/A')}")
    print(f"   Mode: {result['metadata'].get('mode', 'N/A')}")
    print(f"   Latency: {result['metadata']['latency_ms']/1000:.1f}s")

    interp = result.get('interpretation', {})
    print(f"\n🔬 INTERPRETATION:")
    print(f"   Class: {interp.get('overall_class', 'N/A')}")
    print(f"   Confidence: {interp.get('confidence', 0):.2f}")

    # Verify it used gpt-5.2
    model_used = result['metadata'].get('model_selected', result['metadata'].get('model'))
    if model_used == "gpt-5.2":
        print(f"\n✅ SUCCESS: Used GPT-5.2 as expected!")
    else:
        print(f"\n❌ FAILED: Expected gpt-5.2, got {model_used}")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("TEST COMPLETE")
print("="*80)
