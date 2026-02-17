"""
Reaction SMILES Analysis Agent

Combines deterministic cheminformatics (RDKit + rxnmapper) with LLM interpretation
to analyze reaction SMILES and provide mechanistic insights.

Design: docs/reaction_smiles_analysis_agent_simple_v1.md
"""

import json
import re
from typing import Dict, Any, Optional

from llmtools.clients import LLMClient
from reaction_agent.prompts import get_template, get_direct_smiles_template
from reaction_agent.core import analyze_reaction_smiles as analyze_deterministic
from chemtools.quick_reaction_glance import quick_reaction_glance, should_run_quick_glance


def _strip_markdown_fences(text: str) -> str:
    """Remove markdown code fences from JSON output."""
    text = re.sub(r'^```(?:json)?\s*\n', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n```\s*$', '', text, flags=re.MULTILINE)
    return text.strip()


def _format_bond_changes(bond_changes: list) -> str:
    """Format bond changes for prompt."""
    if not bond_changes:
        return "No bond changes detected"

    lines = []
    for bc in bond_changes:
        bc_id = bc.get("id", "?")
        change = bc.get("change", "unknown")
        a1 = bc.get("a1", "?")
        a2 = bc.get("a2", "?")
        bond = bc.get("bond", "unknown")
        lines.append(f"- {bc_id}: {change} bond between atoms :{a1} and :{a2} ({bond})")

    return "\n".join(lines)


def analyze_reaction_smiles(
    rxn_smiles: str,
    client: LLMClient,
    drop_spectators: bool = True,
    skip_mapping: bool = False,
    temperature: float = 0.0,
    max_tokens: int = 2000,
    reasoning_effort: Optional[str] = None
) -> Dict[str, Any]:
    """
    Analyze reaction SMILES using deterministic tools + LLM interpretation.

    This is the main function that orchestrates:
    1. Deterministic cleaning (reaction_agent.core)
    2. LLM interpretation with structured output

    Rxnmapper has been removed - Tier 2 (DeepSeek-v3.2) provides more accurate
    analysis directly from SMILES without atom mapping.

    Args:
        rxn_smiles: Raw reaction SMILES string (reactants>>products)
        client: LLM client instance
        drop_spectators: Remove simple ions/salts from mapping
        skip_mapping: Skip mapping and bond changes (always True, kept for backwards compat)
        temperature: LLM temperature (0.0 for deterministic)
        max_tokens: Max LLM output tokens
        reasoning_effort: For reasoning models (gpt-5.2, o3): "low", "medium", "high"

    Returns:
        Dict with structure:
        {
            "schema_version": "reaction_analysis.v1",
            "input": {...},  # raw/clean SMILES, warnings
            "interpretation": {...},  # LLM analysis
            "quick_glance": {...},  # Tier 2 DeepSeek analysis
            "auto_interpretation": {...},  # Tier 1 string patterns
            "metadata": {...}  # LLM provider, tokens, latency
        }
    """
    # Always skip mapping - Tier 2 is more accurate without it
    skip_mapping = True

    # Step 1: Run deterministic analysis
    try:
        deterministic_result = analyze_deterministic(
            rxn_smiles,
            drop_spectators=drop_spectators,
            skip_mapping=skip_mapping
        )
    except Exception as e:
        return {
            "error": f"Deterministic analysis failed: {e}",
            "schema_version": "reaction_analysis.v1",
            "input": {"rxn_smiles_raw": rxn_smiles}
        }

    input_data = deterministic_result["input"]
    auto_interpretation = deterministic_result.get("auto_interpretation")  # Get automatic interpretation

    # Step 1.5: Quick LLM glance (Tier 2) - optional fast analysis
    quick_glance_result = None

    if auto_interpretation and auto_interpretation.get('interpretation'):
        string_patterns = auto_interpretation['interpretation']

        # Use DeepSeek-v3.2 (Aliyun) for comprehensive chemistry analysis
        # Best accuracy: correctly identifies both Suzuki coupling AND THP deprotection
        # Run on ALL reactions for maximum coverage (prioritizing accuracy over cost)
        if should_run_quick_glance(string_patterns, 1.0, mode="always"):  # No mapping, use confidence=1.0
            try:
                # Create client for quick glance using DeepSeek-v3.2 with comprehensive mode
                quick_client = LLMClient(provider="aliyun", model="deepseek-v3.2")
                quick_glance_result = quick_reaction_glance(
                    input_data.get("rxn_smiles_clean", ""),
                    quick_client,
                    prompt_style="comprehensive",  # Thorough analysis
                    thorough=True  # Enable comprehensive mode
                )

                if quick_glance_result.get('success'):
                    summary = quick_glance_result.get('summary', 'N/A')
                    # Show protecting group changes if detected
                    pg_changes = quick_glance_result.get('protecting_groups', {})
                    if pg_changes.get('removed') or pg_changes.get('added'):
                        pg_note = []
                        if pg_changes.get('removed'):
                            pg_note.append(f"PG removed: {', '.join(pg_changes['removed'][:2])}")
                        if pg_changes.get('added'):
                            pg_note.append(f"PG added: {', '.join(pg_changes['added'][:2])}")
                        print(f"💡 Quick glance: {summary}")
                        print(f"   {' | '.join(pg_note)}")
                    else:
                        print(f"💡 Quick glance: {summary}")

            except Exception as e:
                print(f"⚠️  Quick glance failed: {e}")
                quick_glance_result = {"error": str(e), "success": False}

    # Always use direct SMILES analysis for Tier 3 (no mapping needed)

    # Step 2: Build LLM prompt (direct SMILES analysis, no mapping)
    template = get_direct_smiles_template()

    # Parse reactants and products
    rxn_clean = input_data.get("rxn_smiles_clean", "")
    parts = rxn_clean.split(">>")
    reactants_smiles = parts[0] if len(parts) > 0 else ""
    products_smiles = parts[1] if len(parts) > 1 else ""

    # Include Tier 2 context if available to ensure consistency
    tier2_context = ""
    if quick_glance_result and quick_glance_result.get('success'):
        tier2_rxn_types = quick_glance_result.get('reaction_types', [])
        if tier2_rxn_types:
            tier2_context = f"\n\nIMPORTANT CONTEXT: Quick analysis identified this as: {', '.join(tier2_rxn_types)}. Please verify and provide detailed mechanistic analysis consistent with this."

    prompt = template.format(
        reactants_smiles=reactants_smiles,
        products_smiles=products_smiles,
        rxn_smiles_clean=rxn_clean,
        mapping_confidence=1.0  # Not using mapping
    ) + tier2_context

    # Step 3: Call LLM
    try:
        response = client.chat(
            prompt=prompt,
            temperature=temperature,
            max_tokens=max_tokens,
            reasoning_effort=reasoning_effort
        )
    except Exception as e:
        return {
            "schema_version": "reaction_analysis.v1",
            "input": input_data,
            "error": f"LLM call failed: {e}"
        }

    # Step 4: Parse LLM response
    try:
        # Strip markdown fences if present
        cleaned_response = _strip_markdown_fences(response.content)
        interpretation = json.loads(cleaned_response)
    except json.JSONDecodeError as e:
        # If JSON parsing fails, return raw response with error
        interpretation = {
            "error": f"Failed to parse LLM JSON: {e}",
            "raw_response": response.content,
            "warnings": ["LLM response not valid JSON"],
            "confidence": 0.0
        }

    # Step 5: Assemble final result
    result = {
        "schema_version": "reaction_analysis.v1",
        "input": input_data,
        "interpretation": interpretation,
        "metadata": {
            "model": response.model,
            "provider": response.provider,
            "total_tokens": response.total_tokens,
            "prompt_tokens": response.prompt_tokens,
            "completion_tokens": response.completion_tokens,
            "latency_ms": response.latency_ms,
            "temperature": temperature,
        }
    }

    # Include automatic interpretation if available
    if auto_interpretation:
        result["auto_interpretation"] = auto_interpretation

    # Include quick glance if available
    if quick_glance_result:
        result["quick_glance"] = quick_glance_result

    return result


class ReactionSMILESAnalyzer:
    """
    Agent for analyzing reaction SMILES with LLM interpretation.

    Simple POC implementation combining deterministic cheminformatics with
    LLM reasoning for reaction mechanism and classification.

    Example:
        >>> from llmtools.clients import LLMClient
        >>> from reaction_agent import ReactionSMILESAnalyzer
        >>>
        >>> client = LLMClient(provider="openai", model="gpt-4o")
        >>> analyzer = ReactionSMILESAnalyzer(client)
        >>>
        >>> result = analyzer.analyze("CCBr>>CCN")
        >>> print(result["interpretation"]["overall_class"])
    """

    def __init__(
        self,
        client: LLMClient,
        drop_spectators: bool = True,
        temperature: float = 0.0,
        max_tokens: int = 2000,
        reasoning_effort: Optional[str] = None
    ):
        """
        Initialize analyzer.

        Args:
            client: LLM client for interpretation
            drop_spectators: Remove simple ions/salts from analysis
            temperature: LLM temperature (0.0 = deterministic)
            max_tokens: Max tokens for LLM response
            reasoning_effort: For reasoning models (gpt-5.2, o3): "low", "medium", "high"
        """
        self.client = client
        self.drop_spectators = drop_spectators
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.reasoning_effort = reasoning_effort

    def analyze(
        self,
        rxn_smiles: str,
        mode: str = "auto"
    ) -> Dict[str, Any]:
        """
        Analyze a reaction SMILES string.

        Args:
            rxn_smiles: Reaction SMILES (reactants>>products)
            mode: Analysis mode:
                - "auto": Smart switching based on mapping confidence (recommended)
                          Uses gpt-4o for good mapping (≥0.4), gpt-5.2 for poor mapping (<0.4)
                - "fast": Always use gpt-4o (faster, cheaper)
                - "deep": Always use gpt-5.2 with reasoning (slower, better quality)
                - "expert": gpt-5.2 with highest reasoning (very slow, maximum quality)

        Returns:
            Analysis result dict with input, tool_facts, interpretation, metadata
        """
        # Store original model and settings
        original_model = self.client.model
        original_max_tokens = self.max_tokens
        original_reasoning = self.reasoning_effort

        # Always skip mapping - Tier 2 is more accurate without it
        skip_mapping = True

        # Step 1: Run deterministic analysis (needed for auto mode)
        deterministic_result = analyze_deterministic(
            rxn_smiles,
            drop_spectators=self.drop_spectators,
            skip_mapping=skip_mapping
        )

        # Step 2: Model selection based on mode
        if mode == "auto":
            # Smart switching based on mapping quality
            mapping_conf = deterministic_result.get("tool_facts", {}).get(
                "mapping_qc", {}
            ).get("confidence", 1.0)

            if mapping_conf < 0.4:
                # Poor mapping - use GPT-5.2 with deep reasoning
                print(f"⚠ Poor mapping ({mapping_conf:.3f}) - switching to GPT-5.2 deep reasoning...")
                self.client.model = "gpt-5.2"
                self.max_tokens = 8000
                self.reasoning_effort = "low"  # Test showed "low" is optimal
            else:
                # Good mapping - use fast mode
                self.client.model = original_model if "gpt-5" not in original_model else "gpt-4o"
                self.max_tokens = 3000
                self.reasoning_effort = None

        elif mode == "fast":
            # Always fast - use gpt-4o
            self.client.model = "gpt-4o"
            self.max_tokens = 3000
            self.reasoning_effort = None

        elif mode == "deep":
            # Always deep - use gpt-5.2 with low reasoning
            self.client.model = "gpt-5.2"
            self.max_tokens = 8000
            self.reasoning_effort = "low"

        elif mode == "expert":
            # Expert mode - gpt-5.2 with medium reasoning
            self.client.model = "gpt-5.2"
            self.max_tokens = 10000
            self.reasoning_effort = "medium"

        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'auto', 'fast', 'deep', or 'expert'")

        # Step 3: Run analysis with selected model
        result = analyze_reaction_smiles(
            rxn_smiles=rxn_smiles,
            client=self.client,
            drop_spectators=self.drop_spectators,
            skip_mapping=skip_mapping,
            temperature=self.temperature if self.reasoning_effort is None else None,
            max_tokens=self.max_tokens,
            reasoning_effort=self.reasoning_effort
        )

        # Step 4: Add mode info to metadata
        if "metadata" not in result:
            result["metadata"] = {}

        result["metadata"]["mode"] = mode
        result["metadata"]["model_selected"] = self.client.model
        result["metadata"]["reasoning_effort"] = self.reasoning_effort

        # Restore original settings
        self.client.model = original_model
        self.max_tokens = original_max_tokens
        self.reasoning_effort = original_reasoning

        return result

    def analyze_batch(
        self,
        rxn_smiles_list: list
    ) -> list:
        """
        Analyze multiple reactions.

        Args:
            rxn_smiles_list: List of reaction SMILES

        Returns:
            List of analysis results
        """
        results = []
        for rxn in rxn_smiles_list:
            result = self.analyze(rxn)
            results.append(result)
        return results
