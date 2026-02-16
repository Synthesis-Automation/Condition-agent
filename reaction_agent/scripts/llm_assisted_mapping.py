#!/usr/bin/env python
"""
LLM-Assisted Atom Mapping (Experimental)

Uses LLM reasoning to improve or validate atom mappings when rxnmapper
confidence is low or mapping fails.

Approaches:
1. Validator: LLM checks if rxnmapper mapping makes chemical sense
2. Corrector: LLM suggests improvements to low-confidence mappings
3. Fallback: LLM attempts mapping when rxnmapper fails completely
"""

import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple
import json

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import analyze_deterministic


def validate_mapping_with_llm(
    rxn_smiles: str,
    mapped_smiles: str,
    mapping_confidence: float,
    client: Optional[LLMClient] = None
) -> Dict:
    """
    Use LLM to validate rxnmapper results.

    Args:
        rxn_smiles: Original reaction SMILES
        mapped_smiles: Atom-mapped SMILES from rxnmapper
        mapping_confidence: rxnmapper confidence score
        client: LLM client (uses o3-mini if None)

    Returns:
        Dict with validation results and suggestions
    """

    if client is None:
        client = LLMClient(provider="openai", model="o3-mini")

    prompt = f"""You are a chemistry expert validating atom mappings in reactions.

**Reaction (unmapped)**: {rxn_smiles}
**Atom-mapped version**: {mapped_smiles}
**Tool confidence**: {mapping_confidence:.3f}

Your task: Validate if this atom mapping is chemically reasonable.

**Analysis Steps**:
1. Identify the reaction type and mechanism
2. Check if reactive centers are mapped correctly
3. Verify functional group transformations make sense
4. Look for obvious mapping errors (e.g., oxygen mapped to carbon)
5. Assess if stereochemistry is preserved appropriately

**Return JSON**:
{{
  "mapping_valid": true/false,
  "reaction_type": "cross_coupling/nucleophilic_substitution/etc",
  "issues_found": [
    "list of specific problems, or empty if none"
  ],
  "confidence_adjustment": float (-0.5 to +0.3),
  "recommendation": "trust/review/reject",
  "reasoning": "brief explanation"
}}

Be conservative - only flag clear chemical inconsistencies.
"""

    response = client.chat(prompt=prompt, temperature=0.0)

    # Parse response
    try:
        # Strip markdown fences if present
        content = response.content.strip()
        if content.startswith("```"):
            content = content.split("```")[1]
            if content.startswith("json"):
                content = content[4:]
        content = content.strip()

        validation = json.loads(content)

        # Adjust confidence
        adjusted_confidence = max(0.0, min(1.0,
            mapping_confidence + validation.get('confidence_adjustment', 0.0)))

        return {
            'original_confidence': mapping_confidence,
            'adjusted_confidence': adjusted_confidence,
            'validation': validation,
            'method': 'llm_validation'
        }

    except Exception as e:
        return {
            'error': f"Failed to parse LLM response: {e}",
            'original_confidence': mapping_confidence
        }


def llm_assisted_mapping(
    rxn_smiles: str,
    rxnmapper_result: Dict,
    client: Optional[LLMClient] = None
) -> Dict:
    """
    Attempt to improve mapping using LLM reasoning when rxnmapper fails.

    Args:
        rxn_smiles: Reaction SMILES
        rxnmapper_result: Result from rxnmapper (may have low confidence)
        client: LLM client (uses o3-mini if None)

    Returns:
        Improved mapping or fallback suggestion
    """

    if client is None:
        client = LLMClient(provider="openai", model="o3-mini", max_tokens=8000)

    mapping_conf = rxnmapper_result.get('mapping_qc', {}).get('confidence', 0.0)
    mapped_smiles = rxnmapper_result.get('mapped_rxn_smiles', rxn_smiles)

    prompt = f"""You are a chemistry expert helping with atom mapping in complex reactions.

**Reaction**: {rxn_smiles}
**Current mapping** (from automated tool): {mapped_smiles}
**Tool confidence**: {mapping_conf:.3f} (LOW - needs improvement)

**Your Task**: Use chemical reasoning to identify the most likely atom correspondences.

**Analysis Steps**:
1. Parse reactants and products
2. Identify reaction type(s) - could be multi-stage (e.g., coupling + cyclization)
3. For each stage:
   - Identify reactive centers
   - Trace bond formation/breaking
   - Determine atom correspondences
4. Propose corrected atom mapping if current one is wrong

**Additional Context**:
- Reaction may be tandem/sequential (multiple mechanistic steps)
- Focus on key transformations (ignore spectator ions)
- Use your chemistry knowledge about typical mechanisms

**Return JSON**:
{{
  "reaction_analysis": {{
    "type": "cross_coupling/nucleophilic_substitution/multi_stage/etc",
    "stages": [
      {{"name": "C-N coupling", "mechanism": "Buchwald-Hartwig"}},
      {{"name": "alkyne coupling", "mechanism": "Sonogashira"}}
    ],
    "complexity": "simple/moderate/complex"
  }},
  "mapping_assessment": {{
    "current_mapping_correct": true/false,
    "major_errors": ["list of errors found"],
    "confidence_in_current": 0.0-1.0
  }},
  "suggested_corrections": [
    {{"issue": "atom X mapped to Y but should be Z", "priority": "high/medium/low"}}
  ],
  "recommended_action": "accept/manual_review/expert_needed",
  "reasoning": "detailed explanation of your analysis"
}}

Think step-by-step about the chemistry before responding.
"""

    response = client.chat(prompt=prompt, temperature=0.0)

    try:
        content = response.content.strip()
        if content.startswith("```"):
            content = content.split("```")[1]
            if content.startswith("json"):
                content = content[4:]
        content = content.strip()

        analysis = json.loads(content)

        return {
            'original_confidence': mapping_conf,
            'llm_analysis': analysis,
            'method': 'llm_assisted_mapping',
            'recommendation': analysis.get('recommended_action', 'review')
        }

    except Exception as e:
        return {
            'error': f"Failed to parse LLM response: {e}",
            'original_confidence': mapping_conf
        }


def hybrid_mapping_workflow(rxn_smiles: str, confidence_threshold: float = 0.6) -> Dict:
    """
    Complete workflow: rxnmapper → LLM validation → LLM assistance if needed.

    Args:
        rxn_smiles: Reaction SMILES
        confidence_threshold: Trigger LLM assistance below this threshold

    Returns:
        Comprehensive mapping result with validation
    """

    print(f"\n{'='*80}")
    print("  Hybrid Mapping Workflow")
    print(f"{'='*80}\n")
    print(f"Reaction: {rxn_smiles[:60]}...")

    # Step 1: Standard deterministic analysis
    print("\n[1] Running rxnmapper...")
    det_result = analyze_deterministic(rxn_smiles, skip_mapping=False)

    mapping_qc = det_result.get('tool_facts', {}).get('mapping_qc', {})
    mapping_conf = mapping_qc.get('confidence', 0.0)
    mapping_ok = mapping_qc.get('ok', False)
    mapped_smiles = det_result.get('tool_facts', {}).get('mapped_rxn_smiles', rxn_smiles)

    print(f"    rxnmapper confidence: {mapping_conf:.3f} {'✓' if mapping_ok else '✗'}")

    result = {
        'rxn_smiles': rxn_smiles,
        'deterministic_result': det_result,
        'rxnmapper_confidence': mapping_conf,
        'steps_performed': []
    }

    # Step 2: LLM validation if confidence is borderline
    if 0.4 <= mapping_conf < confidence_threshold:
        print(f"\n[2] Confidence is borderline ({mapping_conf:.3f}) - validating with LLM...")

        validation = validate_mapping_with_llm(
            rxn_smiles,
            mapped_smiles,
            mapping_conf
        )

        result['llm_validation'] = validation
        result['steps_performed'].append('llm_validation')

        adj_conf = validation.get('adjusted_confidence', mapping_conf)
        print(f"    Adjusted confidence: {mapping_conf:.3f} → {adj_conf:.3f}")
        print(f"    Recommendation: {validation.get('validation', {}).get('recommendation', 'N/A')}")

        result['final_confidence'] = adj_conf

    # Step 3: LLM assistance if confidence is very low
    elif mapping_conf < 0.4:
        print(f"\n[2] Confidence is very low ({mapping_conf:.3f}) - requesting LLM assistance...")

        assistance = llm_assisted_mapping(rxn_smiles, det_result['tool_facts'])

        result['llm_assistance'] = assistance
        result['steps_performed'].append('llm_assistance')

        print(f"    LLM recommendation: {assistance.get('recommendation', 'N/A')}")

        # Don't automatically increase confidence, but flag for review
        result['final_confidence'] = mapping_conf
        result['needs_manual_review'] = True

    else:
        print(f"\n[2] Confidence is good ({mapping_conf:.3f}) - skipping LLM validation")
        result['final_confidence'] = mapping_conf
        result['needs_manual_review'] = False

    # Summary
    print(f"\n{'='*80}")
    print("  Summary")
    print(f"{'='*80}")
    print(f"  rxnmapper confidence: {mapping_conf:.3f}")
    print(f"  Final confidence:     {result['final_confidence']:.3f}")
    print(f"  Manual review needed: {result.get('needs_manual_review', False)}")
    print(f"  Steps performed:      {', '.join(result['steps_performed']) or 'none (high confidence)'}")
    print(f"{'='*80}\n")

    return result


def test_hybrid_mapping():
    """Test hybrid mapping on reactions with varying difficulty."""

    if not os.getenv("OPENAI_API_KEY"):
        print("Set OPENAI_API_KEY to run test")
        return

    reactions = [
        {
            "name": "Simple SNAr (expect high confidence)",
            "smiles": "Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1"
        },
        {
            "name": "Complex Tandem (expect low confidence)",
            "smiles": "c1ccc2c(c1)[nH]c1ccccc12.Fc1ccc(I)cc1.C#Cc1ccc(CC)cc1>>CCc1ccc(C#Cc2ccc(-n3c4ccccc4c4ccccc43)cc2)cc1"
        }
    ]

    for rxn in reactions:
        print(f"\n\n{'#'*80}")
        print(f"  Testing: {rxn['name']}")
        print(f"{'#'*80}")

        result = hybrid_mapping_workflow(rxn['smiles'], confidence_threshold=0.6)

        # Save result
        output_dir = Path("reaction_agent/results")
        output_dir.mkdir(parents=True, exist_ok=True)

        filename = rxn['name'].replace(' ', '_').lower() + '_hybrid_mapping.json'
        with open(output_dir / filename, 'w') as f:
            json.dump(result, f, indent=2)

        print(f"\n✓ Results saved to: {output_dir / filename}")


if __name__ == "__main__":
    test_hybrid_mapping()
