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
from reaction_agent.prompts import get_template
from reaction_agent.core import analyze_reaction_smiles as analyze_deterministic


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
    2. Atom mapping via rxnmapper
    3. Bond change extraction
    4. LLM interpretation with structured output

    Args:
        rxn_smiles: Raw reaction SMILES string (reactants>>products)
        client: LLM client instance
        drop_spectators: Remove simple ions/salts from mapping
        skip_mapping: Skip mapping and bond changes (faster, less detailed)
        temperature: LLM temperature (0.0 for deterministic)
        max_tokens: Max LLM output tokens
        reasoning_effort: For reasoning models (gpt-5.2, o3): "low", "medium", "high"

    Returns:
        Dict with structure:
        {
            "schema_version": "reaction_analysis.v1",
            "input": {...},  # raw/clean SMILES, warnings
            "tool_facts": {...},  # mapping, bond changes
            "interpretation": {...},  # LLM analysis
            "metadata": {...}  # LLM provider, tokens, latency
        }
    """
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
    tool_facts = deterministic_result["tool_facts"]

    # Step 2: Build LLM prompt
    template = get_template()

    # Format bond changes for prompt
    bond_changes_text = _format_bond_changes(tool_facts.get("bond_changes", []))

    # Format mapping QC
    mapping_qc = tool_facts.get("mapping_qc", {})
    mapping_qc_text = json.dumps(mapping_qc, indent=2)

    prompt = template.format(
        rxn_smiles_raw=input_data.get("rxn_smiles_raw", ""),
        rxn_smiles_clean=input_data.get("rxn_smiles_clean", ""),
        spectators=", ".join(input_data.get("spectators", [])) or "None",
        parse_warnings=", ".join(input_data.get("parse_warnings", [])) or "None",
        mapped_rxn_smiles=tool_facts.get("mapped_rxn_smiles") or "Not available",
        mapping_qc=mapping_qc_text,
        bond_changes_text=bond_changes_text,
        reaction_center_atoms=str(tool_facts.get("reaction_center_atoms", []))
    )

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
            "tool_facts": tool_facts,
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

    # Step 5: Apply QC gating
    if not mapping_qc.get("ok", False):
        # Degrade confidence if mapping failed
        current_confidence = interpretation.get("confidence", 1.0)
        interpretation["confidence"] = min(current_confidence, 0.3)

        if "warnings" not in interpretation:
            interpretation["warnings"] = []
        if "mapping_failed" not in interpretation["warnings"]:
            interpretation["warnings"].append("mapping_failed")

    # Step 6: Assemble final result
    return {
        "schema_version": "reaction_analysis.v1",
        "input": input_data,
        "tool_facts": tool_facts,
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
        skip_mapping: bool = False,
        mode: str = "auto"
    ) -> Dict[str, Any]:
        """
        Analyze a reaction SMILES string.

        Args:
            rxn_smiles: Reaction SMILES (reactants>>products)
            skip_mapping: Skip atom mapping (faster but less detailed)
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
        rxn_smiles_list: list,
        skip_mapping: bool = False
    ) -> list:
        """
        Analyze multiple reactions.

        Args:
            rxn_smiles_list: List of reaction SMILES
            skip_mapping: Skip atom mapping for all reactions

        Returns:
            List of analysis results
        """
        results = []
        for rxn in rxn_smiles_list:
            result = self.analyze(rxn, skip_mapping=skip_mapping)
            results.append(result)
        return results
