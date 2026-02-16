#!/usr/bin/env python
"""
Quick test of the first complex reaction with top 3 models.
Fast preview while comprehensive testing runs.
"""

import os
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer

# Test reaction: Rare C-N coupling
REACTION = "Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1.[F-][P+5]([F-])([F-])([F-])([F-])[F-]"
RXN_NAME = "Rare C-N Coupling"

def quick_test():
    """Quick test with 3 key models."""

    if not os.getenv("OPENAI_API_KEY"):
        print("Set API key: export OPENAI_API_KEY='sk-...'")
        return

    print("=" * 80)
    print(f"  Quick Test: {RXN_NAME}")
    print("=" * 80)
    print(f"\nReaction: {REACTION[:60]}...\n")

    configs = [
        {"model": "gpt-4o-mini", "tokens": 2500, "desc": "Fast"},
        {"model": "gpt-4o", "tokens": 3000, "desc": "Balanced"},
        {"model": "o3", "tokens": 3500, "desc": "Best Quality"},
    ]

    results = []

    for config in configs:
        print(f"\nTesting {config['desc']}: {config['model']} ({config['tokens']} tokens)")
        print("-" * 60)

        try:
            client = LLMClient(provider="openai", model=config['model'])
            analyzer = ReactionSMILESAnalyzer(client, max_tokens=config['tokens'])

            result = analyzer.analyze(REACTION)
            interp = result['interpretation']

            print(f"  ✓ Complete")
            print(f"    Class: {interp.get('overall_class', 'N/A')}")
            print(f"    Tags: {', '.join(interp.get('tags', []))}")
            print(f"    Confidence: {interp.get('confidence', 0.0):.3f}")
            print(f"    Events: {len(interp.get('events', []))}")
            print(f"    Tokens: {result['metadata']['total_tokens']}")
            print(f"    Time: {result['metadata']['latency_ms']/1000:.1f}s")

            results.append({
                'config': config,
                'confidence': interp.get('confidence', 0.0),
                'class': interp.get('overall_class', 'N/A'),
                'events': len(interp.get('events', [])),
            })

        except Exception as e:
            print(f"  ✗ Failed: {e}")

    # Quick summary
    if results:
        print("\n" + "=" * 80)
        print("  Quick Summary")
        print("=" * 80)

        best = max(results, key=lambda x: x['confidence'])
        print(f"\n✓ Best: {best['config']['desc']} - {best['config']['model']}")
        print(f"  Confidence: {best['confidence']:.3f}")
        print(f"  Class: {best['class']}")
        print(f"  Events: {best['events']}")

if __name__ == "__main__":
    quick_test()
