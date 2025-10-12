"""
LLM-powered enhancements for reaction condition recommendations.

This module provides LLM-based functions to enhance recommendation quality:
- Multi-source synthesis: Combine ML/Rule/Protocol recommendations
- Condition explanations: Generate human-readable rationales
- Re-ranking: Chemistry-aware sorting beyond similarity
- Optimization suggestions: Propose alternative conditions
- Troubleshooting: Diagnose failed reactions

Author: Condition-Agent Team
Date: October 12, 2025
"""

from typing import Dict, List, Any, Optional
import json
import re
from llmtools.clients import LLMClient
from llmtools.prompts import MULTI_SOURCE_SYNTHESIS


def _strip_markdown_fences(text: str) -> str:
    """
    Remove markdown code fences from LLM response.
    
    Args:
        text: Raw LLM response potentially with ```json ... ```
        
    Returns:
        Cleaned text without fences
    """
    # Remove ```json ... ``` or ``` ... ```
    text = re.sub(r'^```(?:json)?\s*\n', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n```\s*$', '', text, flags=re.MULTILINE)
    return text.strip()


def _format_conditions_for_llm(conditions: List[Dict[str, Any]], source_name: str) -> str:
    """
    Format conditions list into readable text for LLM prompt.
    
    Args:
        conditions: List of condition dicts
        source_name: Name of the source (for context)
        
    Returns:
        Formatted string
    """
    if not conditions:
        return f"No {source_name} recommendations available.\n"
    
    lines = []
    for i, cond in enumerate(conditions, 1):
        # Extract key fields
        catalyst = cond.get('catalyst', 'N/A')
        ligand = cond.get('ligand', cond.get('Ligand', 'N/A'))
        solvent = cond.get('solvent', cond.get('Solvent', 'N/A'))
        temp = cond.get('temperature', cond.get('Temperature', 'N/A'))
        base = cond.get('base', cond.get('Base', 'N/A'))
        
        # Get confidence/similarity if available
        confidence = cond.get('confidence', cond.get('similarity', None))
        conf_str = f" (confidence: {confidence:.2f})" if confidence else ""
        
        lines.append(
            f"{i}. Catalyst: {catalyst}, Ligand: {ligand}, Solvent: {solvent}, "
            f"Temp: {temp}, Base: {base}{conf_str}"
        )
    
    return "\n".join(lines)


def synthesize_recommendations_llm(
    reaction_smiles: str,
    ml_results: Optional[Dict[str, Any]] = None,
    rule_results: Optional[Dict[str, Any]] = None,
    protocol_results: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
    llm_client: Optional[LLMClient] = None,
    temperature: float = 0.2,
    max_tokens: int = 1500,
    prompt_version: str = "v2",
) -> Dict[str, Any]:
    """
    Synthesize recommendations from multiple sources using LLM analysis.
    
    This is the core function for multi-source synthesis. It takes results
    from ML-based, rule-based, and protocol-based recommendation engines,
    and uses an LLM to analyze all sources, identify consensus, explain
    discrepancies, and provide a single best recommendation with backups.
    
    Args:
        reaction_smiles: Reaction SMILES string
        ml_results: ML-based recommendations (DRFP similarity)
        rule_results: Rule-based recommendations (SCDB patterns)
        protocol_results: Protocol-based recommendations (literature)
        constraints: Optional user constraints (scale, cost, air_sensitivity, etc.)
        llm_client: Configured LLM client
        temperature: LLM temperature (0.0-1.0, lower = more deterministic)
        max_tokens: Maximum response tokens
        prompt_version: Which prompt to use - "v1" (original) or "v2" (optimized, default)
                       V2 has 30% lower latency with enhanced chemistry guidelines
        
    Returns:
        Dict with:
            - status: "success" or "error"
            - synthesis: LLM analysis and recommendations
            - sources_used: Count of conditions from each source
            - llm_metadata: Model, tokens, latency
            
    Example:
        >>> from llmtools.clients import LLMClient
        >>> llm = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        >>> result = synthesize_recommendations_llm(
        ...     reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        ...     ml_results=ml_recs,
        ...     rule_results=rule_recs,
        ...     protocol_results=protocol_recs,
        ...     llm_client=llm
        ... )
        >>> print(result["synthesis"]["confidence_level"])
        "high"
    """
    # Validate inputs
    if not llm_client:
        return {
            "status": "error",
            "error": "LLM client is required for synthesis"
        }
    
    if not any([ml_results, rule_results, protocol_results]):
        return {
            "status": "error",
            "error": "At least one recommendation source is required"
        }
    
    # Extract conditions from each source
    ml_conditions = []
    rule_conditions = []
    protocol_conditions = []
    
    if ml_results:
        ml_conditions = ml_results.get("recommended_conditions", [])[:3]
    
    if rule_results:
        # Rule results might have different structure
        rule_conditions = rule_results.get("recommended_conditions", [])[:3]
        if not rule_conditions:
            # Try alternative key
            rule_conditions = rule_results.get("conditions", [])[:3]
    
    if protocol_results:
        protocol_conditions = protocol_results.get("recommended_conditions", [])[:3]
    
    # Format for LLM prompt
    ml_text = _format_conditions_for_llm(ml_conditions, "ML-based")
    rule_text = _format_conditions_for_llm(rule_conditions, "Rule-based")
    protocol_text = _format_conditions_for_llm(protocol_conditions, "Protocol-based")
    
    # Format constraints
    if constraints:
        constraints_text = json.dumps(constraints, indent=2)
    else:
        constraints_text = "None specified"
    
    # Select prompt template based on version
    if prompt_version == "v2":
        from llmtools.prompts import MULTI_SOURCE_SYNTHESIS_V2
        prompt_template = MULTI_SOURCE_SYNTHESIS_V2
    else:
        prompt_template = MULTI_SOURCE_SYNTHESIS
    
    # Build prompt
    prompt = prompt_template.format(
        reaction_smiles=reaction_smiles,
        ml_conditions=ml_text,
        rule_conditions=rule_text,
        protocol_conditions=protocol_text,
        constraints=constraints_text
    )
    
    # Call LLM
    try:
        response = llm_client.chat(
            prompt=prompt,
            temperature=temperature,
            max_tokens=max_tokens
        )
    except Exception as exc:
        return {
            "status": "error",
            "error": f"LLM call failed: {exc}"
        }
    
    # Parse response
    raw_content = (response.content or "").strip()
    cleaned_content = _strip_markdown_fences(raw_content)
    
    try:
        synthesis = json.loads(cleaned_content)
    except json.JSONDecodeError as exc:
        return {
            "status": "error",
            "error": f"JSON parse failed: {exc}",
            "raw_response": raw_content
        }
    
    # Success
    return {
        "status": "success",
        "synthesis": synthesis,
        "sources_used": {
            "ml_precedents": len(ml_conditions),
            "rule_matches": len(rule_conditions),
            "protocol_procedures": len(protocol_conditions)
        },
        "llm_metadata": {
            "model": response.model,
            "tokens": response.total_tokens,
            "latency_ms": response.latency_ms
        }
    }


def explain_conditions_llm(
    reaction_smiles: str,
    recommended_conditions: Dict[str, Any],
    llm_client: LLMClient,
    temperature: float = 0.3,
    max_tokens: int = 250,
) -> str:
    """
    Generate human-readable explanation for recommended conditions.
    
    Args:
        reaction_smiles: Reaction SMILES
        recommended_conditions: Dict with catalyst, solvent, temp, etc.
        llm_client: LLM client
        temperature: LLM temperature
        max_tokens: Max response length
        
    Returns:
        2-3 sentence explanation of why these conditions work
        
    Example:
        >>> explain = explain_conditions_llm(
        ...     reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>...",
        ...     recommended_conditions={
        ...         "catalyst": "Pd(PPh3)4",
        ...         "solvent": "toluene",
        ...         "temperature": "80°C",
        ...         "base": "K3PO4"
        ...     },
        ...     llm_client=llm
        ... )
        >>> print(explain)
        "Pd(PPh3)4 is ideal for this Suzuki coupling because..."
    """
    prompt = f"""
    Explain why these reaction conditions are appropriate for this transformation.
    
    Reaction SMILES: {reaction_smiles}
    
    Recommended conditions:
    - Catalyst: {recommended_conditions.get('catalyst', 'N/A')}
    - Ligand: {recommended_conditions.get('ligand', 'N/A')}
    - Solvent: {recommended_conditions.get('solvent', 'N/A')}
    - Temperature: {recommended_conditions.get('temperature', 'N/A')}
    - Base: {recommended_conditions.get('base', 'N/A')}
    
    Provide a concise explanation (2-3 sentences) covering:
    1. Why this catalyst is appropriate for this specific reaction
    2. Why this solvent/temperature combination is chosen
    3. The role of key reagents (base, ligand, additives)
    
    Focus on chemistry rationale, not just repeating the data.
    """
    
    response = llm_client.chat(prompt, temperature=temperature, max_tokens=max_tokens)
    return response.content.strip()


__all__ = [
    "synthesize_recommendations_llm",
    "explain_conditions_llm",
]
