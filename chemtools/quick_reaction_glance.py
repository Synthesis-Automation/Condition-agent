"""
Quick LLM-based reaction recognition (Tier 2 interpretation).

Provides fast, LLM-powered pattern recognition as a middle ground between:
- Tier 1: String patterns (free, instant, limited coverage ~30-40%)
- Tier 2: Quick glance (cheap, 1-3s, good coverage ~80-90%)
- Tier 3: Deep analysis (expensive, 5-30s, complete understanding ~95%+)

This module uses cheap, fast models (gpt-4o-mini, claude-haiku) to quickly
identify reaction types and patterns from SMILES without detailed mechanism.
"""

from typing import Dict, Any, Optional
import json
import re
from llmtools.clients import LLMClient
import time


def quick_reaction_glance(
    rxn_smiles: str,
    client: LLMClient,
    prompt_style: str = "structured",
    thorough: bool = False
) -> Dict[str, Any]:
    """
    Fast LLM-based reaction pattern recognition.

    Uses a cheap model with concise prompt to quickly identify:
    - Main reaction type(s)
    - Key functional groups and patterns
    - Reaction complexity
    - One-line summary

    Args:
        rxn_smiles: Reaction SMILES (reactants>>products)
        client: LLM client (should use cheap model like gpt-4o-mini)
        prompt_style: "structured" or "concise" or "chemistry_expert" or "comprehensive"
        thorough: If True, use comprehensive analysis (more tokens, better for complex reactions)

    Returns:
        Dict with:
        - reaction_types: list of identified reaction types
        - all_transformations: list of ALL detected changes (if thorough=True)
        - protecting_groups: dict with added/removed (if thorough=True)
        - patterns: list of detected patterns
        - complexity: "simple" | "moderate" | "complex" | "tandem"
        - summary: one-line description
        - confidence: 0.0-1.0 estimate
        - metadata: timing, model info
    """

    start_time = time.time()

    # Parse reaction
    if '>>' not in rxn_smiles:
        return {
            'error': 'Invalid SMILES: missing >>',
            'success': False
        }

    reactants, products = rxn_smiles.split('>>')

    # Select prompt and token limit based on mode
    if thorough or prompt_style == "comprehensive":
        prompt = _get_comprehensive_prompt(reactants, products)
        max_tokens = 1000  # More detailed analysis
    elif prompt_style == "structured":
        prompt = _get_structured_prompt(reactants, products)
        max_tokens = 300
    elif prompt_style == "concise":
        prompt = _get_concise_prompt(reactants, products)
        max_tokens = 300
    elif prompt_style == "chemistry_expert":
        prompt = _get_chemistry_expert_prompt(reactants, products)
        max_tokens = 300
    else:
        raise ValueError(f"Unknown prompt_style: {prompt_style}")

    # Call LLM
    try:
        # Check if model is GPT-5 or o-series (needs reasoning_effort)
        is_gpt5_or_o_series = any(
            client.model.startswith(prefix) for prefix in ["gpt-5", "o3", "o4"]
        )

        if is_gpt5_or_o_series:
            # GPT-5/o-series: use reasoning_effort for better analysis
            response = client.chat(
                prompt=prompt,
                max_tokens=max_tokens,
                reasoning_effort="low"  # Low is optimal for quick analysis
            )
        else:
            # Standard models: use temperature
            response = client.chat(
                prompt=prompt,
                temperature=0.0,
                max_tokens=max_tokens
            )
    except Exception as e:
        return {
            'error': f'LLM call failed: {e}',
            'success': False
        }

    elapsed = time.time() - start_time

    # Parse response
    try:
        # Try to extract JSON from response
        content = response.content.strip()

        # Remove markdown fences if present
        content = re.sub(r'^```(?:json)?\s*\n', '', content, flags=re.MULTILINE)
        content = re.sub(r'\n```\s*$', '', content, flags=re.MULTILINE)
        content = content.strip()

        result = json.loads(content)

        # Validate required fields
        required_fields = ['reaction_types', 'patterns', 'complexity', 'summary']
        for field in required_fields:
            if field not in result:
                result[field] = "unknown" if field in ['complexity', 'summary'] else []

        # Add metadata
        result['success'] = True
        result['metadata'] = {
            'model': response.model,
            'provider': response.provider,
            'latency_ms': response.latency_ms,
            'total_time_seconds': elapsed,
            'total_tokens': response.total_tokens,
            'prompt_style': prompt_style
        }

        return result

    except json.JSONDecodeError as e:
        return {
            'error': f'Failed to parse JSON: {e}',
            'raw_response': response.content[:500],
            'success': False,
            'metadata': {
                'model': response.model,
                'latency_ms': response.latency_ms,
                'total_time_seconds': elapsed
            }
        }


def _get_structured_prompt(reactants: str, products: str) -> str:
    """Structured prompt with clear JSON schema."""
    return f"""Quickly identify this organic reaction. Be concise.

Reactants: {reactants}
Products: {products}

Analyze and return JSON with:
1. reaction_types: List main reaction type(s) (e.g., "Suzuki coupling", "SN2", "oxidation")
2. patterns: Key functional groups/patterns detected (e.g., "aryl bromide", "boronic ester", "acetal")
3. complexity: "simple", "moderate", "complex", or "tandem"
4. summary: One-line description (max 80 chars)
5. confidence: Your confidence 0.0-1.0 in this identification

Return ONLY valid JSON in this format:
{{
  "reaction_types": ["type1", "type2"],
  "patterns": ["pattern1", "pattern2", "pattern3"],
  "complexity": "simple",
  "summary": "Brief description of the transformation",
  "confidence": 0.9
}}"""


def _get_concise_prompt(reactants: str, products: str) -> str:
    """Ultra-concise prompt for fastest response."""
    return f"""Identify this reaction in JSON format:

Reactants: {reactants}
Products: {products}

{{
  "reaction_types": ["list reaction types"],
  "patterns": ["list key patterns"],
  "complexity": "simple|moderate|complex|tandem",
  "summary": "one line",
  "confidence": 0.0-1.0
}}

Return only JSON."""


def _get_chemistry_expert_prompt(reactants: str, products: str) -> str:
    """Chemistry expert persona for better pattern recognition."""
    return f"""You are an expert organic chemist. Quickly recognize this reaction's type and key features.

Reactants: {reactants}
Products: {products}

Identify:
- What reaction type(s) is this? (e.g., Suzuki coupling, SN2, oxidation, condensation)
- What key functional groups/patterns do you see? (e.g., aryl halide, carbonyl, acetal)
- Is this simple, moderate, complex, or a tandem/multi-step reaction?
- Summarize in one sentence (max 80 characters)

Return as JSON:
{{
  "reaction_types": ["primary type", "secondary type if tandem"],
  "patterns": ["pattern1", "pattern2", "pattern3"],
  "complexity": "simple|moderate|complex|tandem",
  "summary": "Concise summary",
  "confidence": 0.0-1.0
}}

JSON only, no explanation."""


def _get_comprehensive_prompt(reactants: str, products: str) -> str:
    """Comprehensive chemistry analysis prompt for thorough examination."""
    return f"""You are an expert organic chemist. Analyze this reaction THOROUGHLY and COMPREHENSIVELY.

Reactants: {reactants}
Products: {products}

Carefully examine and think step-by-step:

1. STRUCTURAL COMPARISON:
   - What structures/patterns appear in reactants but NOT in products? (removed/transformed)
   - What structures/patterns appear in products but NOT in reactants? (added/formed)

2. PROTECTING GROUPS (very important!):
   - Look for: THP (tetrahydropyran like C1CCCCO1), Boc, TBDMS, Bn, Ac, Cbz, acetals, etc.
   - Were any protecting groups ADDED? (protection)
   - Were any protecting groups REMOVED? (deprotection)

3. MAIN TRANSFORMATION:
   - Primary reaction type (Suzuki coupling, SN2, oxidation, etc.)
   - C-C bond formation? Bond breaking?

4. ADDITIONAL TRANSFORMATIONS:
   - Side reactions
   - Workup transformations (acid/base labile groups)
   - Is this tandem/sequential?

Return detailed JSON:
{{
  "main_transformation": "primary reaction type",
  "reaction_types": ["all reaction types"],
  "all_changes": ["list EVERY structural change you notice"],
  "protecting_groups": {{
    "removed": ["any deprotections with details"],
    "added": ["any protections with details"]
  }},
  "side_reactions": ["any additional transformations"],
  "complexity": "simple|moderate|complex|tandem",
  "summary": "Complete description including ALL transformations",
  "confidence": 0.0-1.0,
  "reasoning": "brief explanation of your analysis"
}}

CRITICAL: Respond with ONLY valid JSON. No markdown (```), no text before or after. Start with {{ and end with }}.
JSON only, nothing else."""


def should_run_quick_glance(
    string_patterns_result: Dict[str, Any],
    mapping_confidence: float = 1.0,
    mode: str = "auto"
) -> bool:
    """
    Decide whether to run quick LLM glance based on Tier 1 results.

    Args:
        string_patterns_result: Result from interpret_reaction_pattern()
        mapping_confidence: Confidence from rxnmapper
        mode: "always", "auto", "never", "if_uncertain"

    Returns:
        True if quick glance should be run
    """

    if mode == "always":
        return True

    if mode == "never":
        return False

    if mode == "auto":
        # Run if string patterns found nothing
        if len(string_patterns_result.get('likely_reaction_types', [])) == 0:
            return True

        # Run if only 1 pattern detected (low confidence)
        if len(string_patterns_result.get('patterns_detected', [])) <= 1:
            return True

        # Run if mapping confidence is borderline
        if 0.4 <= mapping_confidence <= 0.6:
            return True

        return False

    if mode == "if_uncertain":
        # Only run if string patterns failed completely
        if len(string_patterns_result.get('likely_reaction_types', [])) == 0:
            return True
        return False

    return False


def format_quick_glance_report(result: Dict[str, Any]) -> str:
    """
    Format quick glance result for display.

    Args:
        result: Result from quick_reaction_glance()

    Returns:
        Formatted text report
    """

    if not result.get('success'):
        return f"Quick glance failed: {result.get('error', 'Unknown error')}"

    lines = []
    lines.append("=" * 80)
    lines.append("QUICK LLM GLANCE (Tier 2)")
    lines.append("=" * 80)
    lines.append("")

    # Summary
    summary = result.get('summary', 'N/A')
    lines.append(f"Summary: {summary}")
    lines.append("")

    # Reaction types
    reaction_types = result.get('reaction_types', [])
    if reaction_types:
        lines.append("Reaction Type(s):")
        for rxn_type in reaction_types:
            lines.append(f"  • {rxn_type}")
        lines.append("")

    # Patterns
    patterns = result.get('patterns', [])
    if patterns:
        lines.append("Key Patterns:")
        for pattern in patterns:
            lines.append(f"  • {pattern}")
        lines.append("")

    # Complexity
    complexity = result.get('complexity', 'unknown')
    confidence = result.get('confidence', 0.0)
    lines.append(f"Complexity: {complexity.upper()}")
    lines.append(f"Confidence: {confidence:.2f}")
    lines.append("")

    # Metadata
    if 'metadata' in result:
        meta = result['metadata']
        lines.append(f"Model: {meta.get('model', 'N/A')}")
        lines.append(f"Latency: {meta.get('latency_ms', 0):.0f} ms")
        lines.append(f"Tokens: {meta.get('total_tokens', 0)}")

    lines.append("=" * 80)

    return '\n'.join(lines)
