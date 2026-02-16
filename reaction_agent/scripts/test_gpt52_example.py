"""
Test GPT-5.2 on Tosylhydrazone Example from NEW-1.md

This is the exact reaction GPT-5.2 analyzed via web interface.
We'll test different reasoning levels and parameters to match/exceed that quality.
"""

import sys
import json
import time
from pathlib import Path

# Add both reaction_agent and parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))  # reaction_agent
sys.path.insert(0, str(Path(__file__).parent.parent.parent))  # Condition-agent

from llmtools.clients import LLMClient
from reaction_agent.core import analyze_reaction_smiles as analyze_deterministic

# The tosylhydrazone + iodonium example from NEW-1.md (line 3)
TOSYLHYDRAZONE_REACTION = {
    "name": "Tosylhydrazone + Diaryliodonium → N-aryl pyrazole + diaryl sulfone",
    "smiles": "Cc1ccc(S(=O)(=O)NN=CC=Cc2ccccc2)cc1.O=S(=O)([O-])C(F)(F)F.c1ccc([I+]c2ccccc2)cc1>>Cc1ccc(S(=O)(=O)c2ccccc2)cc1.c1ccc(-c2ccn(-c3ccccc3)n2)cc1",
    "expected_findings": [
        "N-aryl pyrazole synthesis",
        "diaryl sulfone formation",
        "N-S bond cleavage (tosyl leaves)",
        "pyrazole ring construction + N-arylation",
        "S-Ar formation (sulfone)",
        "Cu-catalyzed oxidative cyclization",
        "α,β-unsaturated N-tosylhydrazone substrate",
        "diaryliodonium reagent"
    ]
}


def test_with_config(model: str, reasoning_level: str = None, max_tokens: int = 8000):
    """Test a single configuration."""
    rxn_smiles = TOSYLHYDRAZONE_REACTION["smiles"]

    print(f"\n{'='*80}")
    print(f"Testing: {model}" + (f" (reasoning={reasoning_level})" if reasoning_level else ""))
    print(f"Max tokens: {max_tokens}")
    print(f"{'='*80}")

    # Step 1: Deterministic analysis
    print("\n[1/2] Running deterministic analysis...")
    det_result = analyze_deterministic(rxn_smiles, drop_spectators=True, skip_mapping=False)

    tool_facts = det_result.get("tool_facts", {})
    mapping_conf = tool_facts.get("mapping_qc", {}).get("confidence", 0.0)

    # Handle bond_changes - could be dict or list
    bond_changes_data = tool_facts.get("bond_changes", {})
    if isinstance(bond_changes_data, dict):
        bonds_broken = bond_changes_data.get('bonds_broken', [])
        bonds_formed = bond_changes_data.get('bonds_formed', [])
        reaction_center = bond_changes_data.get('reaction_center_atoms', [])
    else:
        bonds_broken = []
        bonds_formed = []
        reaction_center = []

    print(f"      Mapping confidence: {mapping_conf:.3f}")
    print(f"      Bonds broken: {len(bonds_broken)}")
    print(f"      Bonds formed: {len(bonds_formed)}")

    # Step 2: LLM analysis
    print(f"\n[2/2] Running LLM analysis...")

    # Build custom prompt that encourages detailed analysis like GPT-5.2
    custom_prompt = f"""Analyze this chemical reaction in detail.

## Reaction SMILES
{rxn_smiles}

## Deterministic Analysis
- Atom mapping confidence: {mapping_conf:.3f}
- Bonds broken: {len(bonds_broken)}
- Bonds formed: {len(bonds_formed)}
- Reaction center atoms: {len(reaction_center)} atoms

## Required Analysis

Provide a comprehensive analysis including:

1. **Reaction Identification**
   - What is the reaction class/name?
   - What are the key substrates and reagents?
   - Is this a known transformation? (reference if known)

2. **Mechanism Breakdown**
   - What are the distinct reaction centers/stages?
   - Which bonds break and form?
   - What intermediates are likely formed?
   - Is this a cascade/tandem reaction?

3. **SMARTS Pattern Detection** (for automated classification)
   - Provide SMARTS patterns for key substrate motifs (e.g., tosylhydrazone, iodonium)
   - Provide SMARTS patterns for products (e.g., pyrazole, sulfone)
   - Suggest scoring logic: +N points for observation X

4. **Product Analysis**
   - What products are formed?
   - Are there multiple products from one starting material?
   - What is the mechanistic relationship between products?

5. **Evidence Scoring**
   - Break down your confidence by evidence type
   - Assign numeric scores to each observation
   - Total score and classification threshold

Return as JSON with this structure:
{{
    "reaction_identification": {{
        "reaction_class": "...",
        "key_substrates": ["...", "..."],
        "reagents": ["..."],
        "known_transformation": true/false,
        "literature_precedent": "description if known"
    }},
    "mechanism": {{
        "reaction_centers": [
            {{
                "center_id": 1,
                "description": "...",
                "bonds_broken": ["..."],
                "bonds_formed": ["..."]
            }}
        ],
        "is_cascade": true/false,
        "intermediates": ["..."]
    }},
    "smarts_patterns": {{
        "substrates": {{
            "pattern_name": "SMARTS string"
        }},
        "products": {{
            "pattern_name": "SMARTS string"
        }},
        "scoring_logic": {{
            "observation": score
        }}
    }},
    "products": {{
        "product_list": ["..."],
        "relationships": "description"
    }},
    "evidence_scores": [
        {{"observation": "...", "score": 3, "source": "bond_change"}}
    ],
    "total_evidence_score": 10,
    "confidence": 0.9,
    "warnings": ["..."]
}}
"""

    try:
        # Use longer timeout for reasoning models
        timeout = 300 if reasoning_level else 60  # 5 min for reasoning, 1 min for standard

        client = LLMClient(provider="openai", model=model, timeout=timeout)

        start = time.time()
        response = client.chat(
            prompt=custom_prompt,
            reasoning_effort=reasoning_level,
            max_tokens=max_tokens
        )
        elapsed = time.time() - start

        print(f"      ✓ Completed in {elapsed:.1f}s")
        print(f"      Tokens: {response.total_tokens} ({response.prompt_tokens} prompt + {response.completion_tokens} completion)")

        # Try to parse JSON
        content = response.content.strip()
        if content.startswith("```"):
            lines = content.split("\n")
            if lines[0].startswith("```"):
                lines = lines[1:]
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            content = "\n".join(lines)

        try:
            result = json.loads(content)

            # Check against expected findings
            print(f"\n{'='*80}")
            print("ANALYSIS QUALITY CHECK")
            print(f"{'='*80}")

            findings_found = 0
            for expected in TOSYLHYDRAZONE_REACTION["expected_findings"]:
                content_lower = json.dumps(result).lower()
                if any(term in content_lower for term in expected.lower().split()):
                    findings_found += 1
                    print(f"      ✓ Found: {expected}")
                else:
                    print(f"      ✗ Missing: {expected}")

            coverage = findings_found / len(TOSYLHYDRAZONE_REACTION["expected_findings"])
            print(f"\n      Coverage: {findings_found}/{len(TOSYLHYDRAZONE_REACTION['expected_findings'])} ({coverage:.1%})")

            return {
                "model": model,
                "reasoning_level": reasoning_level,
                "max_tokens": max_tokens,
                "time": elapsed,
                "tokens": response.total_tokens,
                "prompt_tokens": response.prompt_tokens,
                "completion_tokens": response.completion_tokens,
                "coverage": coverage,
                "findings_found": findings_found,
                "result": result,
                "raw_content": response.content[:500]  # First 500 chars for debugging
            }

        except json.JSONDecodeError as e:
            print(f"      ⚠ Failed to parse JSON: {e}")
            print(f"      Raw output preview: {content[:200]}...")
            return {
                "model": model,
                "reasoning_level": reasoning_level,
                "error": "JSON parse failed",
                "raw_content": content[:500]
            }

    except Exception as e:
        print(f"      ❌ Error: {e}")
        return {
            "model": model,
            "reasoning_level": reasoning_level,
            "error": str(e)
        }


def main():
    """Test multiple configurations to find optimal settings."""
    print("="*80)
    print("GPT-5.2 PARAMETER OPTIMIZATION TEST")
    print("="*80)
    print(f"\nReaction: {TOSYLHYDRAZONE_REACTION['name']}")
    print(f"Expected findings: {len(TOSYLHYDRAZONE_REACTION['expected_findings'])}")
    print("\nTesting configurations to match GPT-5.2 web interface quality...")
    print("="*80)

    # Test configurations (progressive)
    configs = [
        # Baseline
        ("gpt-4o", None, 3000),

        # GPT-5.2 with different reasoning levels
        ("gpt-5.2", "low", 8000),
        ("gpt-5.2", "medium", 10000),
        ("gpt-5.2", "high", 16000),

        # Try o3-mini as alternative
        ("o3-mini", "medium", 12000),
    ]

    results = []

    for model, reasoning, max_tokens in configs:
        try:
            result = test_with_config(model, reasoning, max_tokens)
            results.append(result)

            # Short pause between requests
            time.sleep(2)

        except Exception as e:
            print(f"\n❌ Configuration failed: {e}")
            results.append({
                "model": model,
                "reasoning_level": reasoning,
                "error": str(e)
            })

    # Save results
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)

    output_file = output_dir / "gpt52_parameter_optimization.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n✓ Results saved to: {output_file}")

    # Summary
    print("\n" + "="*80)
    print("OPTIMIZATION RESULTS")
    print("="*80)

    successful = [r for r in results if "error" not in r and "coverage" in r]

    if successful:
        # Sort by coverage
        successful.sort(key=lambda x: x["coverage"], reverse=True)

        print(f"\nBest Configuration:")
        best = successful[0]
        print(f"  Model: {best['model']}")
        print(f"  Reasoning: {best.get('reasoning_level', 'N/A')}")
        print(f"  Coverage: {best['coverage']:.1%} ({best['findings_found']}/{len(TOSYLHYDRAZONE_REACTION['expected_findings'])})")
        print(f"  Time: {best['time']:.1f}s")
        print(f"  Tokens: {best['tokens']}")
        print(f"  Cost estimate: ${best['tokens'] * 0.000015:.4f}")  # Rough estimate

        print(f"\nAll Configurations (by coverage):")
        print(f"{'Model':<15} {'Reasoning':<10} {'Coverage':<10} {'Time':<10} {'Tokens':<10}")
        print("-" * 70)
        for r in successful:
            model = r['model']
            reasoning = r.get('reasoning_level') or 'N/A'  # Handle None
            coverage = f"{r['coverage']:.1%}"
            time_str = f"{r['time']:.1f}s"
            tokens = str(r['tokens'])
            print(f"{model:<15} {reasoning:<10} {coverage:<10} {time_str:<10} {tokens:<10}")

        # Recommendation
        print("\n" + "="*80)
        print("RECOMMENDATIONS")
        print("="*80)

        if best['coverage'] >= 0.75:
            print(f"✓ Excellent: {best['model']} with reasoning={best.get('reasoning_level')} achieves {best['coverage']:.0%} coverage")
            print(f"  Use this configuration for reactions with mapping <0.4")
        elif best['coverage'] >= 0.5:
            print(f"~ Good: {best['model']} achieves {best['coverage']:.0%} coverage")
            print(f"  Consider increasing reasoning level or max_tokens")
        else:
            print(f"⚠ Limited: Best coverage is only {best['coverage']:.0%}")
            print(f"  May need prompt refinement or different model")

        # Check if GPT-5.2 is better than baseline
        baseline = next((r for r in results if r.get('model') == 'gpt-4o'), None)
        gpt52_best = next((r for r in successful if r.get('model') == 'gpt-5.2'), None)

        if baseline and gpt52_best and 'coverage' in baseline:
            improvement = gpt52_best['coverage'] - baseline['coverage']
            cost_ratio = gpt52_best['time'] / baseline['time'] if baseline['time'] > 0 else 0

            print(f"\nGPT-5.2 vs Baseline (gpt-4o):")
            print(f"  Coverage improvement: {improvement:+.1%}")
            print(f"  Time cost: {cost_ratio:.1f}x slower")

            if improvement > 0.1:
                print(f"  ✓ Significant improvement - worth the extra cost")
            elif improvement > 0:
                print(f"  ~ Modest improvement - may be worth it for critical reactions")
            else:
                print(f"  ✗ No improvement - stick with baseline")

    else:
        print("\n❌ No successful tests completed")
        print("Check API access and error messages above")


if __name__ == "__main__":
    main()
