#!/usr/bin/env python
"""
Comprehensive model comparison for complex multi-stage reaction.

Testing C-N coupling + Sonogashira reaction:
c1ccc2c(c1)[nH]c1ccccc12.Fc1ccc(I)cc1.C#Cc1ccc(CC)cc1>>CCc1ccc(C#Cc2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1

Systematically compares different models and parameters.
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer


# The complex reaction to test
REACTION = "c1ccc2c(c1)[nH]c1ccccc12.Fc1ccc(I)cc1.C#Cc1ccc(CC)cc1>>CCc1ccc(C#Cc2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1"

REACTION_NAME = "Carbazole_CN_Sonogashira_tandem"


def print_separator(title=""):
    """Print a section separator."""
    if title:
        print(f"\n{'=' * 80}")
        print(f"  {title}")
        print('=' * 80)
    else:
        print('-' * 80)


def analyze_with_config(model: str, provider: str = "openai",
                       temperature: float = 0.0, max_tokens: int = 3000,
                       skip_mapping: bool = False) -> dict:
    """Analyze reaction with specific configuration."""

    config_name = f"{provider}/{model}"
    print(f"\n{config_name}")
    print(f"  Temperature: {temperature}, Max Tokens: {max_tokens}, Skip Mapping: {skip_mapping}")

    try:
        # Initialize client
        client = LLMClient(
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=max_tokens
        )

        # Create analyzer
        analyzer = ReactionSMILESAnalyzer(
            client=client,
            temperature=temperature,
            max_tokens=max_tokens
        )

        # Analyze
        result = analyzer.analyze(REACTION, skip_mapping=skip_mapping)

        # Extract key metrics
        interp = result.get('interpretation', {})
        metadata = result.get('metadata', {})

        print(f"  ✓ Complete")
        print(f"    Class: {interp.get('overall_class', 'N/A')}")
        print(f"    Tags: {', '.join(interp.get('tags', []))}")
        print(f"    Confidence: {interp.get('confidence', 0.0):.3f}")
        print(f"    Events: {len(interp.get('events', []))}")
        print(f"    Tokens: {metadata.get('total_tokens', 0)}")
        print(f"    Latency: {metadata.get('latency_ms', 0):.0f}ms")

        if interp.get('warnings'):
            print(f"    ⚠ Warnings: {len(interp['warnings'])}")

        return {
            "config": {
                "model": model,
                "provider": provider,
                "temperature": temperature,
                "max_tokens": max_tokens,
                "skip_mapping": skip_mapping
            },
            "result": result,
            "summary": {
                "overall_class": interp.get('overall_class', 'unknown'),
                "tags": interp.get('tags', []),
                "confidence": interp.get('confidence', 0.0),
                "num_events": len(interp.get('events', [])),
                "num_mechanism_steps": len(interp.get('mechanism_summary', [])),
                "warnings": interp.get('warnings', []),
                "tokens": metadata.get('total_tokens', 0),
                "latency_ms": metadata.get('latency_ms', 0)
            }
        }

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        return {
            "config": {
                "model": model,
                "provider": provider,
                "temperature": temperature,
                "max_tokens": max_tokens
            },
            "error": str(e)
        }


def run_comparison():
    """Run comprehensive model comparison."""

    print_separator("Complex Multi-Stage Reaction Analysis")
    print(f"\nReaction: {REACTION_NAME}")
    print(f"SMILES: {REACTION[:80]}...")
    print(f"\nDescription:")
    print("  - Reactant 1: Carbazole (carbazole)")
    print("  - Reactant 2: 4-fluoro-1-iodobenzene")
    print("  - Reactant 3: 4-ethylphenylacetylene")
    print("  - Expected: C-N coupling (Buchwald-Hartwig) + Sonogashira coupling")

    # Check API key
    if not os.getenv("OPENAI_API_KEY"):
        print("\n✗ Error: OPENAI_API_KEY not set")
        print("Set it: export OPENAI_API_KEY='sk-...'")
        return

    results = []

    # Test configurations
    print_separator("TEST CONFIGURATIONS")

    configs = [
        # Fast models - baseline
        {
            "name": "Baseline: Fast & Cheap",
            "model": "gpt-4o-mini",
            "temperature": 0.0,
            "max_tokens": 2000
        },

        # Same but more tokens for complex reaction
        {
            "name": "Fast with More Tokens",
            "model": "gpt-4o-mini",
            "temperature": 0.0,
            "max_tokens": 3500
        },

        # Better quality models
        {
            "name": "High Quality: GPT-4o",
            "model": "gpt-4o",
            "temperature": 0.0,
            "max_tokens": 3000
        },

        {
            "name": "High Quality: More Tokens",
            "model": "gpt-4o",
            "temperature": 0.0,
            "max_tokens": 4000
        },

        # Reasoning models
        {
            "name": "Reasoning: o3-mini",
            "model": "o3-mini",
            "temperature": 0.0,  # Ignored for o-series
            "max_tokens": 3000
        },

        {
            "name": "Reasoning: o3-mini (More Tokens)",
            "model": "o3-mini",
            "temperature": 0.0,
            "max_tokens": 4000
        },

        # Latest GPT-5 series
        {
            "name": "Latest: GPT-5.2",
            "model": "gpt-5.2",
            "temperature": 0.0,  # Ignored for GPT-5
            "max_tokens": 3000
        },

        {
            "name": "Latest: GPT-5.2 (More Tokens)",
            "model": "gpt-5.2",
            "temperature": 0.0,
            "max_tokens": 4000
        },

        # Advanced reasoning
        {
            "name": "Most Capable: o3",
            "model": "o3",
            "temperature": 0.0,
            "max_tokens": 3500
        },
    ]

    print(f"\nTesting {len(configs)} configurations...")

    for i, config in enumerate(configs, 1):
        print_separator(f"Test {i}/{len(configs)}: {config['name']}")

        result = analyze_with_config(
            model=config['model'],
            temperature=config['temperature'],
            max_tokens=config['max_tokens']
        )

        results.append({
            "test_name": config['name'],
            **result
        })

    # Summary comparison
    print_separator("RESULTS SUMMARY")

    successful = [r for r in results if 'error' not in r]

    if not successful:
        print("\n✗ All tests failed")
        return

    print(f"\nSuccessful: {len(successful)}/{len(results)}")
    print("\nComparison Table:")
    print("-" * 120)
    print(f"{'Model':<20} {'Class':<25} {'Events':<8} {'Conf':<6} {'Tokens':<8} {'Latency':<10} {'Tags'}")
    print("-" * 120)

    for r in successful:
        config = r['config']
        summary = r['summary']

        model_str = config['model'][:19]
        class_str = summary['overall_class'][:24]
        tags_str = ', '.join(summary['tags'][:3])[:40]

        print(f"{model_str:<20} {class_str:<25} {summary['num_events']:<8} "
              f"{summary['confidence']:<6.3f} {summary['tokens']:<8} "
              f"{summary['latency_ms']:<10.0f} {tags_str}")

    # Find best by different criteria
    print("\n" + "=" * 80)
    print("BEST CONFIGURATIONS")
    print("=" * 80)

    # Best by confidence
    best_conf = max(successful, key=lambda x: x['summary']['confidence'])
    print(f"\n✓ Highest Confidence: {best_conf['test_name']}")
    print(f"  Model: {best_conf['config']['model']}")
    print(f"  Confidence: {best_conf['summary']['confidence']:.3f}")
    print(f"  Class: {best_conf['summary']['overall_class']}")

    # Best by detail (most events)
    best_detail = max(successful, key=lambda x: x['summary']['num_events'])
    print(f"\n✓ Most Detailed: {best_detail['test_name']}")
    print(f"  Model: {best_detail['config']['model']}")
    print(f"  Events: {best_detail['summary']['num_events']}")
    print(f"  Mechanism Steps: {best_detail['summary']['num_mechanism_steps']}")

    # Best value (best confidence per dollar)
    best_value = min(successful, key=lambda x: x['summary']['tokens'] / max(x['summary']['confidence'], 0.1))
    print(f"\n✓ Best Value: {best_value['test_name']}")
    print(f"  Model: {best_value['config']['model']}")
    print(f"  Confidence: {best_value['summary']['confidence']:.3f}")
    print(f"  Tokens: {best_value['summary']['tokens']}")

    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = project_root / "reaction_agent" / "results" / f"comparison_{REACTION_NAME}_{timestamp}.json"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"\n✓ Detailed results saved to: {output_file}")

    # Recommendation
    print_separator("RECOMMENDATION")

    if best_conf['summary']['confidence'] >= 0.7:
        print(f"\n✅ RECOMMENDED: {best_conf['test_name']}")
        print(f"   Model: {best_conf['config']['model']}")
        print(f"   Confidence: {best_conf['summary']['confidence']:.3f}")
        print(f"   Tokens: {best_conf['summary']['tokens']} (~${best_conf['summary']['tokens'] * 0.000001:.4f})")
        print(f"\nCommand to reproduce:")
        print(f'   python reaction_agent/cli.py \\')
        print(f'     --reaction "{REACTION}" \\')
        print(f'     --model {best_conf["config"]["model"]} \\')
        print(f'     --max-tokens {best_conf["config"]["max_tokens"]}')
    else:
        print(f"\n⚠ All models show relatively low confidence (<0.7)")
        print(f"   This is expected for complex multi-stage reactions")
        print(f"\n   Best option: {best_conf['test_name']}")
        print(f"   Consider manual verification of results")

    return results


if __name__ == "__main__":
    try:
        results = run_comparison()
        print("\n✓ Comparison complete!")
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
