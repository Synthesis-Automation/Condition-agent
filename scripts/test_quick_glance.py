"""
Comprehensive test of quick LLM glance with different models and prompts.

Tests:
- Different models (gpt-4o-mini, gpt-4o)
- Different prompt styles (structured, concise, chemistry_expert)
- Different reaction types (simple, complex, tandem)
- Performance metrics (latency, accuracy, cost)
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools.quick_reaction_glance import quick_reaction_glance, format_quick_glance_report
from llmtools.clients import LLMClient
import time
import json


# Test reactions with known answers
TEST_REACTIONS = [
    {
        "name": "Tandem (Suzuki + Acetal Hydrolysis)",
        "smiles": "COC(OC)c1ccccc1Br.C=CC(CCCC)B1OC(C)(C)C(C)(C)O1>>C=CC(CCCC)c1ccccc1C=O",
        "expected": {
            "reaction_types": ["Suzuki coupling", "acetal hydrolysis"],
            "complexity": "tandem",
            "patterns": ["aryl bromide", "boronic ester", "acetal", "aldehyde"]
        }
    },
    {
        "name": "Simple SN2",
        "smiles": "CCCBr.N>>CCCN",
        "expected": {
            "reaction_types": ["nucleophilic substitution", "SN2"],
            "complexity": "simple",
            "patterns": ["alkyl bromide", "amine"]
        }
    },
    {
        "name": "Complex Suzuki Coupling",
        "smiles": "CC1(C)CC=C(B2OC(C)(C)C(C)(C)O2)CC1.CS(=O)(=O)N1CCN(C(=O)c2cnc3ccc(F)cc3c2Br)CC1>>CC1(C)CC=C(c2c(C(=O)N3CCN(S(C)(=O)=O)CC3)cnc3ccc(F)cc23)CC1",
        "expected": {
            "reaction_types": ["Suzuki coupling", "cross-coupling"],
            "complexity": "moderate",
            "patterns": ["boronic ester", "aryl bromide", "heterocycle"]
        }
    },
    {
        "name": "Esterification",
        "smiles": "CC(=O)O.CCO>>CC(=O)OCC",
        "expected": {
            "reaction_types": ["esterification", "condensation"],
            "complexity": "simple",
            "patterns": ["carboxylic acid", "alcohol", "ester"]
        }
    },
    {
        "name": "Grignard Addition",
        "smiles": "CC(=O)C.CCCMgBr>>CC(O)(CCC)C",
        "expected": {
            "reaction_types": ["Grignard addition", "nucleophilic addition"],
            "complexity": "simple",
            "patterns": ["ketone", "Grignard reagent", "alcohol"]
        }
    }
]


def score_result(result: dict, expected: dict) -> dict:
    """Score LLM result against expected answer."""

    score = {
        'total': 0.0,
        'reaction_types_score': 0.0,
        'complexity_score': 0.0,
        'patterns_score': 0.0,
        'details': {}
    }

    if not result.get('success'):
        return score

    # Score reaction types (40% weight)
    detected_types = [t.lower() for t in result.get('reaction_types', [])]
    expected_types = [t.lower() for t in expected.get('reaction_types', [])]

    matches = sum(1 for exp in expected_types if any(exp in det for det in detected_types))
    if expected_types:
        score['reaction_types_score'] = matches / len(expected_types)
    score['details']['reaction_types_match'] = f"{matches}/{len(expected_types)}"

    # Score complexity (20% weight)
    detected_complexity = result.get('complexity', '').lower()
    expected_complexity = expected.get('complexity', '').lower()
    if detected_complexity == expected_complexity:
        score['complexity_score'] = 1.0
    score['details']['complexity_match'] = f"{detected_complexity} vs {expected_complexity}"

    # Score patterns (40% weight)
    detected_patterns = [p.lower() for p in result.get('patterns', [])]
    expected_patterns = [p.lower() for p in expected.get('patterns', [])]

    matches = sum(1 for exp in expected_patterns if any(exp in det for det in detected_patterns))
    if expected_patterns:
        score['patterns_score'] = matches / len(expected_patterns)
    score['details']['patterns_match'] = f"{matches}/{len(expected_patterns)}"

    # Total score
    score['total'] = (
        0.4 * score['reaction_types_score'] +
        0.2 * score['complexity_score'] +
        0.4 * score['patterns_score']
    )

    return score


def test_configuration(
    model: str,
    prompt_style: str,
    test_reactions: list
) -> dict:
    """Test a specific model + prompt configuration."""

    print(f"\n{'='*100}")
    print(f"Testing: {model} + {prompt_style}")
    print(f"{'='*100}\n")

    client = LLMClient(provider="openai", model=model)

    results = []
    total_time = 0
    total_tokens = 0
    total_cost = 0.0

    for test in test_reactions:
        print(f"Testing: {test['name']}")
        print(f"SMILES: {test['smiles'][:80]}...")
        print()

        start = time.time()
        result = quick_reaction_glance(
            test['smiles'],
            client,
            prompt_style=prompt_style
        )
        elapsed = time.time() - start

        if result.get('success'):
            print(f"✓ Success in {elapsed:.2f}s")
            print(f"  Reaction types: {result.get('reaction_types', [])}")
            print(f"  Complexity: {result.get('complexity', 'N/A')}")
            print(f"  Patterns: {result.get('patterns', [])}")
            print(f"  Summary: {result.get('summary', 'N/A')}")
            print()

            # Score result
            score = score_result(result, test['expected'])
            print(f"  Score: {score['total']:.2f} (types: {score['reaction_types_score']:.2f}, "
                  f"complexity: {score['complexity_score']:.2f}, patterns: {score['patterns_score']:.2f})")
            print()

            # Accumulate metrics
            meta = result.get('metadata', {})
            total_time += elapsed
            total_tokens += meta.get('total_tokens', 0)

            # Estimate cost (approximate)
            if 'gpt-4o-mini' in model:
                cost = (meta.get('total_tokens', 0) / 1000) * 0.00025  # ~$0.25/1M tokens average
            elif 'gpt-4o' in model:
                cost = (meta.get('total_tokens', 0) / 1000) * 0.005  # ~$5/1M tokens average
            else:
                cost = (meta.get('total_tokens', 0) / 1000) * 0.001

            total_cost += cost

            results.append({
                'test': test['name'],
                'success': True,
                'score': score,
                'latency': elapsed,
                'tokens': meta.get('total_tokens', 0),
                'cost': cost
            })

        else:
            print(f"✗ Failed: {result.get('error', 'Unknown error')}")
            print()

            results.append({
                'test': test['name'],
                'success': False,
                'error': result.get('error')
            })

    # Summary
    successful = [r for r in results if r.get('success')]
    avg_score = sum(r['score']['total'] for r in successful) / len(successful) if successful else 0
    avg_latency = sum(r['latency'] for r in successful) / len(successful) if successful else 0

    print(f"\n{'='*100}")
    print(f"SUMMARY: {model} + {prompt_style}")
    print(f"{'='*100}")
    print(f"Success rate: {len(successful)}/{len(results)}")
    print(f"Average score: {avg_score:.3f}")
    print(f"Average latency: {avg_latency:.2f}s")
    print(f"Total time: {total_time:.2f}s")
    print(f"Total tokens: {total_tokens}")
    print(f"Total cost: ${total_cost:.4f}")
    print()

    return {
        'model': model,
        'prompt_style': prompt_style,
        'results': results,
        'summary': {
            'success_rate': len(successful) / len(results),
            'avg_score': avg_score,
            'avg_latency': avg_latency,
            'total_time': total_time,
            'total_tokens': total_tokens,
            'total_cost': total_cost
        }
    }


def main():
    """Run comprehensive tests."""

    print("="*100)
    print("COMPREHENSIVE QUICK GLANCE TESTING")
    print("="*100)
    print()
    print(f"Testing {len(TEST_REACTIONS)} reactions with different models and prompts")
    print()

    # Configurations to test
    configs = [
        ("gpt-4o-mini", "structured"),
        ("gpt-4o-mini", "concise"),
        ("gpt-4o-mini", "chemistry_expert"),
        ("gpt-4o", "structured"),
        ("gpt-4o", "chemistry_expert"),
    ]

    all_results = []

    for model, prompt_style in configs:
        try:
            result = test_configuration(model, prompt_style, TEST_REACTIONS)
            all_results.append(result)
        except Exception as e:
            print(f"\n✗ Configuration failed: {model} + {prompt_style}")
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            print()

    # Final comparison
    print("\n" + "="*100)
    print("FINAL COMPARISON")
    print("="*100)
    print()

    # Sort by score
    all_results.sort(key=lambda x: x['summary']['avg_score'], reverse=True)

    print(f"{'Model':<20} {'Prompt':<20} {'Score':<8} {'Latency':<10} {'Cost':<10}")
    print("-" * 100)

    for result in all_results:
        model = result['model']
        prompt = result['prompt_style']
        score = result['summary']['avg_score']
        latency = result['summary']['avg_latency']
        cost = result['summary']['total_cost']

        print(f"{model:<20} {prompt:<20} {score:.3f}    {latency:.2f}s     ${cost:.4f}")

    print()

    # Recommendation
    if all_results:
        best = all_results[0]
        print("RECOMMENDATION:")
        print(f"  Best configuration: {best['model']} + {best['prompt_style']}")
        print(f"  Average score: {best['summary']['avg_score']:.3f}")
        print(f"  Average latency: {best['summary']['avg_latency']:.2f}s")
        print(f"  Cost per reaction: ${best['summary']['total_cost'] / len(TEST_REACTIONS):.5f}")
        print()

        # Cost-effectiveness analysis
        print("COST-EFFECTIVENESS:")
        for result in all_results[:3]:  # Top 3
            score_per_dollar = result['summary']['avg_score'] / (result['summary']['total_cost'] + 0.0001)
            score_per_second = result['summary']['avg_score'] / result['summary']['avg_latency']
            print(f"  {result['model']} + {result['prompt_style']}")
            print(f"    Score/$ = {score_per_dollar:.1f}, Score/sec = {score_per_second:.3f}")

    # Save detailed results
    output_file = "results/quick_glance_test_results.json"
    import os
    os.makedirs("results", exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    print()
    print(f"✓ Detailed results saved to {output_file}")
    print()


if __name__ == "__main__":
    main()
