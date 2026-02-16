"""
Compare Fast vs Deep Reasoning Modes

Tests the same reactions with:
1. Fast mode: gpt-4o, single-shot, T=0.0
2. Deep mode: o3-mini, extended reasoning, T=1.0, multi-turn

Hypothesis: Deep mode will better identify complex mechanisms,
especially for reactions with poor atom mapping (<0.4 confidence).
"""

import sys
import json
from pathlib import Path
from typing import Dict, Any, List

# Add reaction_agent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from llmtools.clients import LLMClient
from core import analyze_reaction_smiles as analyze_deterministic
from agent import ReactionSMILESAnalyzer


# Deep reasoning prompt template - encourages hypothesis exploration
DEEP_REASONING_TEMPLATE = """You are analyzing a chemical reaction. Use extended reasoning to explore multiple mechanistic hypotheses before concluding.

## Reaction SMILES
Original: {rxn_smiles_raw}
Cleaned: {rxn_smiles_clean}

## Deterministic Analysis (from RDKit + rxnmapper)

{tool_facts_section}

## Your Task

**Phase 1: Hypothesis Generation (explore multiple possibilities)**
1. What reaction classes could this be? List 3-5 possibilities based on:
   - Functional groups present
   - Reagent types
   - Bond changes observed
   - Known reaction precedents

**Phase 2: Evidence Evaluation**
For each hypothesis, score evidence (0-5 scale):
- Direct bond changes: +3
- Motif matches: +2
- Supporting observations: +1
- Inconsistencies: -2

**Phase 3: Mechanism Verification**
Consider:
- Are all bond changes explained?
- Are byproducts/spectators consistent?
- Does stoichiometry make sense?
- Any unexplained atoms or groups?

**Phase 4: Self-Correction**
- What are the weaknesses in your top hypothesis?
- What alternative explanations exist?
- What would you need to be more confident?

## Output Format

Return JSON with:
{{
    "hypotheses_explored": [
        {{
            "mechanism": "hypothesis description",
            "evidence_score": 8.5,
            "pros": ["..."],
            "cons": ["..."]
        }}
    ],
    "final_classification": {{
        "reaction_class": "most likely class",
        "mechanism_type": "type",
        "evidence_breakdown": [
            {{"observation": "...", "score": 3, "source": "bond_change"}}
        ],
        "total_evidence_score": 8.5
    }},
    "confidence": 0.85,
    "reasoning_summary": "Brief summary of why the final answer was chosen",
    "uncertainties": ["What remains uncertain"],
    "warnings": ["Potential issues or caveats"]
}}

**Important**: Show your iterative reasoning. Explore alternatives before committing.
"""


class DeepReasoningAnalyzer:
    """Extended reasoning mode for complex reactions."""

    def __init__(self, client: LLMClient, max_tokens: int = 12000):
        self.client = client
        self.max_tokens = max_tokens

    def analyze(self, rxn_smiles: str, deterministic_result: Dict) -> Dict[str, Any]:
        """
        Analyze reaction with extended reasoning.

        Uses higher temperature (1.0) and larger token budget
        to allow hypothesis exploration and self-correction.
        """
        # Extract tool facts
        tool_facts = deterministic_result.get("tool_facts", {})

        # Build comprehensive tool facts section
        facts_lines = []
        facts_lines.append("### Cleaning Report")
        if "clean" in tool_facts:
            clean = tool_facts["clean"]
            facts_lines.append(f"- Original: {clean.get('raw_smiles', '')}")
            facts_lines.append(f"- Cleaned: {clean.get('clean_smiles', '')}")
            if clean.get("spectators"):
                facts_lines.append(f"- Spectators removed: {', '.join(clean['spectators'])}")

        facts_lines.append("\n### Atom Mapping Quality")
        if "mapping_qc" in tool_facts:
            qc = tool_facts["mapping_qc"]
            facts_lines.append(f"- Confidence: {qc.get('confidence', 0):.3f}")
            facts_lines.append(f"- Status: {'✓ Reliable' if qc.get('ok') else '⚠ Unreliable - treat mechanism with caution'}")
            if not qc.get('ok'):
                facts_lines.append(f"- Issue: {qc.get('issue', 'unknown')}")

        facts_lines.append("\n### Bond Changes Detected")
        if "bond_changes" in tool_facts:
            bc = tool_facts["bond_changes"]
            if bc.get("bonds_broken"):
                facts_lines.append(f"- Bonds broken: {bc['bonds_broken']}")
            if bc.get("bonds_formed"):
                facts_lines.append(f"- Bonds formed: {bc['bonds_formed']}")
            if bc.get("bonds_order_changed"):
                facts_lines.append(f"- Bond order changed: {bc['bonds_order_changed']}")
            if bc.get("reaction_center_atoms"):
                facts_lines.append(f"- Reaction center: {len(bc['reaction_center_atoms'])} atoms involved")

        tool_facts_section = "\n".join(facts_lines)

        # Build prompt
        prompt = DEEP_REASONING_TEMPLATE.format(
            rxn_smiles_raw=deterministic_result.get("input_smiles", rxn_smiles),
            rxn_smiles_clean=tool_facts.get("clean", {}).get("clean_smiles", rxn_smiles),
            tool_facts_section=tool_facts_section
        )

        # Call LLM with extended reasoning settings
        response = self.client.chat(
            prompt=prompt,
            temperature=1.0,  # Allow exploration (vs 0.0 in fast mode)
            max_tokens=self.max_tokens
        )

        # Parse response
        content = response.content.strip()

        # Strip markdown fences if present
        if content.startswith("```"):
            lines = content.split("\n")
            # Remove first line if it's ```json or ```
            if lines[0].startswith("```"):
                lines = lines[1:]
            # Remove last line if it's ```
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            content = "\n".join(lines)

        try:
            interpretation = json.loads(content)
        except json.JSONDecodeError as e:
            interpretation = {
                "error": f"Failed to parse JSON: {e}",
                "raw_output": content[:500],
                "confidence": 0.0
            }

        # Combine with deterministic results
        result = deterministic_result.copy()
        result["interpretation"] = interpretation
        result["mode"] = "deep_reasoning"

        return result


def compare_modes(
    rxn_smiles: str,
    rxn_name: str = None,
    fast_model: str = "gpt-4o",
    deep_model: str = "o3-mini"
) -> Dict:
    """
    Compare fast vs deep reasoning on a single reaction.

    Returns:
    {
        "reaction": reaction_name,
        "fast_mode": {...},
        "deep_mode": {...},
        "comparison": {
            "confidence_diff": float,
            "evidence_score": {fast, deep},
            "time_ratio": float,
            "cost_ratio": float
        }
    }
    """
    print(f"\n{'='*80}")
    print(f"Testing: {rxn_name or rxn_smiles[:60]}")
    print(f"{'='*80}")

    # Step 1: Deterministic analysis (shared by both modes)
    print("\n[1/3] Running deterministic analysis...")
    det_result = analyze_deterministic(rxn_smiles, drop_spectators=True, skip_mapping=False)

    mapping_conf = det_result.get('tool_facts', {}).get('mapping_qc', {}).get('confidence', 0.0)
    print(f"      Mapping confidence: {mapping_conf:.3f}")

    # Step 2: Fast mode (current system)
    print(f"\n[2/3] Running FAST mode ({fast_model})...")
    fast_client = LLMClient(provider="openai", model=fast_model)
    fast_analyzer = ReactionSMILESAnalyzer(fast_client, temperature=0.0, max_tokens=3000)

    import time
    start = time.time()
    fast_result = fast_analyzer.analyze(rxn_smiles, skip_mapping=False)
    fast_time = time.time() - start

    fast_conf = fast_result.get('interpretation', {}).get('confidence', 0.0)
    fast_class = fast_result.get('interpretation', {}).get('reaction_class', 'Unknown')
    print(f"      Confidence: {fast_conf:.3f}")
    print(f"      Class: {fast_class}")
    print(f"      Time: {fast_time:.1f}s")

    # Step 3: Deep mode (extended reasoning)
    print(f"\n[3/3] Running DEEP mode ({deep_model})...")
    deep_client = LLMClient(provider="openai", model=deep_model)
    deep_analyzer = DeepReasoningAnalyzer(deep_client, max_tokens=12000)

    start = time.time()
    deep_result = deep_analyzer.analyze(rxn_smiles, det_result)
    deep_time = time.time() - start

    deep_interp = deep_result.get('interpretation', {})
    deep_conf = deep_interp.get('confidence', 0.0)
    deep_final = deep_interp.get('final_classification', {})
    deep_class = deep_final.get('reaction_class', 'Unknown')
    deep_evidence = deep_final.get('total_evidence_score', 0.0)
    num_hypotheses = len(deep_interp.get('hypotheses_explored', []))

    print(f"      Confidence: {deep_conf:.3f}")
    print(f"      Class: {deep_class}")
    print(f"      Evidence score: {deep_evidence}")
    print(f"      Hypotheses explored: {num_hypotheses}")
    print(f"      Time: {deep_time:.1f}s")

    # Comparison
    comparison = {
        "confidence_diff": deep_conf - fast_conf,
        "time_ratio": deep_time / fast_time if fast_time > 0 else 0,
        "hypotheses_explored": num_hypotheses,
        "deep_evidence_score": deep_evidence,
        "fast_more_confident": fast_conf > deep_conf,
        "agreement": fast_class.lower() == deep_class.lower()
    }

    print(f"\n{'='*80}")
    print("COMPARISON")
    print(f"{'='*80}")
    print(f"Confidence change: {comparison['confidence_diff']:+.3f}")
    print(f"Time ratio: {comparison['time_ratio']:.1f}x slower")
    print(f"Classes agree: {comparison['agreement']}")
    print(f"Fast more confident: {comparison['fast_more_confident']}")

    return {
        "reaction": rxn_name or rxn_smiles[:60],
        "smiles": rxn_smiles,
        "mapping_confidence": mapping_conf,
        "fast_mode": {
            "model": fast_model,
            "confidence": fast_conf,
            "class": fast_class,
            "time": fast_time
        },
        "deep_mode": {
            "model": deep_model,
            "confidence": deep_conf,
            "class": deep_class,
            "evidence_score": deep_evidence,
            "hypotheses_explored": num_hypotheses,
            "time": deep_time
        },
        "comparison": comparison
    }


def test_on_problematic_reactions():
    """
    Test on reactions where rxnmapper failed (mapping <0.4).
    These are the cases where deep reasoning should help most.
    """
    # From test_random_30.py results - reactions with poor mapping
    test_cases = [
        {
            "name": "Reaction 2 - Aziridination",
            "smiles": "[Cl-].CCCCCCCCCCCCCCCC[N+](C)(C)CCl.N/N=C/c1ccccc1>>c1ccc(C2CN2)cc1",
            "mapping": 0.087  # Very poor
        },
        {
            "name": "Reaction 3 - Alkene Bromination",
            "smiles": "BrBr.CC(C)=CCC/C(C)=C/CO>>CC(C)=CCC/C(C)=C/CBr",
            "mapping": 0.274
        },
        {
            "name": "Reaction 14 - Buchwald-Hartwig",
            "smiles": "[Cl-].[Na+].Cc1ccc(Br)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1",
            "mapping": 0.268
        },
        {
            "name": "Reaction 19 - Complex heterocycle",
            "smiles": "[Na+].[OH-].COC(=O)c1ccc(Cl)cc1Cl.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccc(Cl)cc1Cl",
            "mapping": 0.351
        },
        {
            "name": "Reaction 27 - Tosylhydrazone (from previous test)",
            "smiles": "Cc1ccc(S(=O)(=O)NN=CC=Cc2ccccc2)cc1.O=S(=O)([O-])C(F)(F)F.c1ccc([I+]c2ccccc2)cc1>>Cc1ccc(S(=O)(=O)c2ccccc2)cc1.c1ccc(-c2ccn(-c3ccccc3)n2)cc1",
            "mapping": 0.34  # GPT-5.2's example
        }
    ]

    results = []

    for case in test_cases:
        try:
            result = compare_modes(
                rxn_smiles=case["smiles"],
                rxn_name=case["name"],
                fast_model="gpt-4o",
                deep_model="o3-mini"
            )
            results.append(result)

        except Exception as e:
            print(f"\n❌ Error on {case['name']}: {e}")
            results.append({
                "reaction": case["name"],
                "error": str(e)
            })

    return results


def main():
    """Run comparison tests and save results."""
    print("="*80)
    print("REASONING MODE COMPARISON TEST")
    print("="*80)
    print("\nComparing:")
    print("- Fast mode: gpt-4o, T=0.0, 3000 tokens")
    print("- Deep mode: o3-mini, T=1.0, 12000 tokens")
    print("\nTesting on reactions with poor atom mapping (<0.4)")
    print("="*80)

    results = test_on_problematic_reactions()

    # Save results
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)

    output_file = output_dir / "reasoning_mode_comparison.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n✓ Results saved to: {output_file}")

    # Summary statistics
    successful = [r for r in results if "error" not in r]

    if successful:
        avg_conf_diff = sum(r["comparison"]["confidence_diff"] for r in successful) / len(successful)
        avg_time_ratio = sum(r["comparison"]["time_ratio"] for r in successful) / len(successful)
        deep_better = sum(1 for r in successful if r["comparison"]["confidence_diff"] > 0.1)
        agreements = sum(1 for r in successful if r["comparison"]["agreement"])

        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(f"Tests completed: {len(successful)}/{len(results)}")
        print(f"Average confidence improvement: {avg_conf_diff:+.3f}")
        print(f"Deep mode substantially better (>0.1): {deep_better}/{len(successful)}")
        print(f"Classifications agree: {agreements}/{len(successful)}")
        print(f"Average time cost: {avg_time_ratio:.1f}x slower")
        print("="*80)

        # Recommendation
        if avg_conf_diff > 0.1:
            print("\n✓ RECOMMENDATION: Deep mode provides significant improvement")
            print("  Use for reactions with mapping confidence <0.4")
        elif avg_conf_diff > 0:
            print("\n~ RECOMMENDATION: Deep mode provides modest improvement")
            print("  Consider using for high-value analyses")
        else:
            print("\n✗ RECOMMENDATION: Deep mode does not improve over fast mode")
            print("  Investigate if prompt needs refinement")


if __name__ == "__main__":
    main()
