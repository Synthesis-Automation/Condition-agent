"""
Test GPT-5.2 with Different Reasoning Levels

Directly compares:
1. gpt-4o (fast, no reasoning)
2. gpt-5.2 with low reasoning
3. gpt-5.2 with medium reasoning
4. gpt-5.2 with high reasoning

Much simpler than building custom deep reasoning prompts!
"""

import sys
import json
import time
from pathlib import Path
from typing import Dict, Any, List

# Add reaction_agent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from llmtools.clients import LLMClient
from core import analyze_reaction_smiles as analyze_deterministic
from agent import ReactionSMILESAnalyzer


def test_reaction_with_reasoning_levels(
    rxn_smiles: str,
    rxn_name: str = None
) -> Dict:
    """
    Test a single reaction with multiple reasoning levels.

    Returns comparison of:
    - gpt-4o (baseline)
    - gpt-5.2 low reasoning
    - gpt-5.2 medium reasoning
    - gpt-5.2 high reasoning
    """
    print(f"\n{'='*80}")
    print(f"Testing: {rxn_name or rxn_smiles[:60]}")
    print(f"{'='*80}")

    # Step 1: Deterministic analysis (shared)
    print("\n[1/5] Running deterministic analysis...")
    det_result = analyze_deterministic(rxn_smiles, drop_spectators=True, skip_mapping=False)

    mapping_conf = det_result.get('tool_facts', {}).get('mapping_qc', {}).get('confidence', 0.0)
    print(f"      Mapping confidence: {mapping_conf:.3f}")

    results = {
        "reaction": rxn_name or rxn_smiles[:60],
        "smiles": rxn_smiles,
        "mapping_confidence": mapping_conf,
        "tests": {}
    }

    # Test configurations
    configs = [
        ("gpt-4o", "baseline", None),
        ("gpt-5.2", "low", "low"),
        ("gpt-5.2", "medium", "medium"),
        ("gpt-5.2", "high", "high"),
    ]

    for idx, (model, label, reasoning_level) in enumerate(configs, start=2):
        print(f"\n[{idx}/5] Testing {model} ({label})...")

        try:
            client = LLMClient(provider="openai", model=model)
            analyzer = ReactionSMILESAnalyzer(
                client,
                temperature=0.0 if model == "gpt-4o" else None,  # GPT-5 uses built-in temp=1
                max_tokens=8000 if model == "gpt-5.2" else 3000
            )

            start = time.time()

            # For GPT-5.2, we need to pass reasoning_effort through to the client
            # Since ReactionSMILESAnalyzer doesn't support it yet, we'll patch it
            if reasoning_level:
                # Temporarily override the client's chat method to include reasoning_effort
                original_chat = client.chat

                def chat_with_reasoning(*args, **kwargs):
                    kwargs['reasoning_effort'] = reasoning_level
                    return original_chat(*args, **kwargs)

                client.chat = chat_with_reasoning

            result = analyzer.analyze(rxn_smiles, skip_mapping=False)
            elapsed = time.time() - start

            interp = result.get('interpretation', {})
            confidence = interp.get('confidence', 0.0)
            reaction_class = interp.get('reaction_class', 'Unknown')
            mechanism = interp.get('mechanism_type', '')

            print(f"      Confidence: {confidence:.3f}")
            print(f"      Class: {reaction_class}")
            print(f"      Time: {elapsed:.1f}s")

            results["tests"][label] = {
                "model": model,
                "reasoning_level": reasoning_level,
                "confidence": confidence,
                "reaction_class": reaction_class,
                "mechanism_type": mechanism,
                "time": elapsed,
                "full_interpretation": interp
            }

        except Exception as e:
            print(f"      ❌ Error: {e}")
            results["tests"][label] = {
                "model": model,
                "reasoning_level": reasoning_level,
                "error": str(e)
            }

    # Comparison analysis
    print(f"\n{'='*80}")
    print("COMPARISON")
    print(f"{'='*80}")

    baseline_conf = results["tests"].get("baseline", {}).get("confidence", 0)

    for label in ["low", "medium", "high"]:
        if label in results["tests"] and "error" not in results["tests"][label]:
            test = results["tests"][label]
            conf_diff = test["confidence"] - baseline_conf
            time_ratio = test["time"] / results["tests"]["baseline"]["time"]

            print(f"\nGPT-5.2 {label.upper()}:")
            print(f"  Confidence: {test['confidence']:.3f} ({conf_diff:+.3f} vs baseline)")
            print(f"  Time: {test['time']:.1f}s ({time_ratio:.1f}x baseline)")
            print(f"  Class: {test['reaction_class']}")

    return results


def main():
    """Test on problematic reactions from our dataset."""
    print("="*80)
    print("GPT-5.2 REASONING LEVEL COMPARISON")
    print("="*80)
    print("\nTesting:")
    print("- Baseline: gpt-4o (no reasoning)")
    print("- GPT-5.2 with low/medium/high reasoning effort")
    print("\nOn reactions with poor atom mapping (<0.4)")
    print("="*80)

    # Test cases - reactions where rxnmapper struggled
    test_cases = [
        {
            "name": "Tosylhydrazone + Iodonium → Pyrazole + Sulfone (GPT-5.2's example)",
            "smiles": "Cc1ccc(S(=O)(=O)NN=CC=Cc2ccccc2)cc1.O=S(=O)([O-])C(F)(F)F.c1ccc([I+]c2ccccc2)cc1>>Cc1ccc(S(=O)(=O)c2ccccc2)cc1.c1ccc(-c2ccn(-c3ccccc3)n2)cc1",
            "mapping": 0.34
        },
        {
            "name": "Aziridination",
            "smiles": "[Cl-].CCCCCCCCCCCCCCCC[N+](C)(C)CCl.N/N=C/c1ccccc1>>c1ccc(C2CN2)cc1",
            "mapping": 0.087
        },
        {
            "name": "Buchwald-Hartwig Amination",
            "smiles": "[Cl-].[Na+].Cc1ccc(Br)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1",
            "mapping": 0.268
        },
    ]

    all_results = []

    for case in test_cases:
        try:
            result = test_reaction_with_reasoning_levels(
                rxn_smiles=case["smiles"],
                rxn_name=case["name"]
            )
            all_results.append(result)

        except Exception as e:
            print(f"\n❌ Error on {case['name']}: {e}")
            all_results.append({
                "reaction": case["name"],
                "error": str(e)
            })

    # Save results
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)

    output_file = output_dir / "gpt52_reasoning_levels.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\n✓ Results saved to: {output_file}")

    # Summary
    successful = [r for r in all_results if "error" not in r]
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    if successful:
        # Compare average improvements
        for reasoning_level in ["low", "medium", "high"]:
            improvements = []
            for result in successful:
                baseline = result["tests"].get("baseline", {}).get("confidence", 0)
                test = result["tests"].get(reasoning_level, {})
                if "confidence" in test:
                    improvements.append(test["confidence"] - baseline)

            if improvements:
                avg_improvement = sum(improvements) / len(improvements)
                better_count = sum(1 for x in improvements if x > 0.1)
                print(f"\nGPT-5.2 {reasoning_level.upper()}:")
                print(f"  Average confidence gain: {avg_improvement:+.3f}")
                print(f"  Substantially better (>0.1): {better_count}/{len(improvements)}")

        print("\n" + "="*80)

        # Find best reasoning level
        best_level = None
        best_improvement = -999
        for level in ["low", "medium", "high"]:
            improvements = []
            for result in successful:
                baseline = result["tests"].get("baseline", {}).get("confidence", 0)
                test = result["tests"].get(level, {})
                if "confidence" in test:
                    improvements.append(test["confidence"] - baseline)
            if improvements:
                avg = sum(improvements) / len(improvements)
                if avg > best_improvement:
                    best_improvement = avg
                    best_level = level

        if best_level:
            print(f"\n✓ RECOMMENDATION: Use GPT-5.2 with '{best_level}' reasoning")
            print(f"  Average improvement: {best_improvement:+.3f} confidence")
            print(f"  Use for reactions with mapping <0.4 (27% of dataset)")
        else:
            print("\n~ Inconclusive: Need more testing to determine optimal level")

    else:
        print("\n❌ No successful tests completed")


if __name__ == "__main__":
    main()
