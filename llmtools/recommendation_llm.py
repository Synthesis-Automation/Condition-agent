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
from chemtools.recommend.utils import friendly_family_label


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
    
    Handles multiple condition formats:
    - ML/Protocol format: flat keys (catalyst, ligand, solvent, etc.)
    - Rule-based format: chemicals array with roles + nested conditions
    
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
        # Check if this is rule-based format (has 'chemicals' array)
        if 'chemicals' in cond:
            # Rule-based format: extract from chemicals array
            catalyst = 'N/A'
            ligand = 'N/A'
            solvent = 'N/A'
            base = 'N/A'
            
            for chem in cond.get('chemicals', []):
                role = chem.get('role', '').lower()
                name = chem.get('name') or chem.get('abbreviation', 'N/A')
                
                if role in ['metal_catalyst', 'catalyst', 'metal_source']:
                    catalyst = name
                elif role == 'ligand':
                    ligand = name
                elif role == 'solvent':
                    solvent = name
                elif role == 'base':
                    base = name
            
            # Extract temperature from nested conditions
            nested_cond = cond.get('conditions', {})
            temp = nested_cond.get('temperature', 'N/A')
            if isinstance(temp, list) and len(temp) == 2:
                temp = f"{temp[0]}-{temp[1]}°C"
            elif isinstance(temp, (int, float)):
                temp = f"{temp}°C"
            
        else:
            # Flat format: direct key access
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


def convert_llm_synthesis_to_standard_format(
    reaction_smiles: str,
    synthesis_result: Dict[str, Any],
    requested_type: Optional[str] = None,
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Convert LLM synthesis result to standard output format for robotic execution.
    
    Transforms the analysis-focused LLM output into the same structure as
    ML/Rule/Protocol outputs, making it compatible with robotic systems.
    
    Args:
        reaction_smiles: Reaction SMILES string
        synthesis_result: Result from synthesize_recommendations_llm()
        requested_type: User-requested reaction type
        processing_time_ms: Total processing time
        
    Returns:
        Standard format output with meta, input, detection, recommended_conditions
        
    Example output structure:
        {
            "meta": {"model": "LLM-synthesis", ...},
            "input": {"reaction_smiles": "...", ...},
            "detection": {"family": "Suzuki_Miyaura", ...},
            "recommended_conditions": [
                {
                    "rank": 1,
                    "chemicals": [...],
                    "conditions": {...},
                    "confidence": 0.95,
                    "reasoning": {...}
                }
            ]
        }
    """
    from chemtools.output_formatter import (
        format_meta,
        format_input,
        format_detection,
        normalize_chemical_entry,
        normalize_conditions_block,
        starting_material_entries,
    )
    
    if synthesis_result.get("status") != "success":
        # Return error format
        return {
            "meta": format_meta(
                model_type="LLM-synthesis",
                status="error",
                processing_time_ms=processing_time_ms,
            ),
            "input": format_input(reaction_smiles=reaction_smiles),
            "error": synthesis_result.get("error", "Unknown error"),
        }
    
    synthesis = synthesis_result.get("synthesis", {})
    llm_metadata = synthesis_result.get("llm_metadata", {})
    sources_used = synthesis_result.get("sources_used", {})
    
    # Detect reaction type from synthesis or fallback
    detected_type = requested_type or "unknown"
    
    # Extract confidence level
    confidence_level = synthesis.get("confidence_level", "medium")
    confidence_map = {"high": 0.95, "medium": 0.75, "low": 0.50}
    confidence_score = confidence_map.get(confidence_level, 0.75)
    
    # Build recommended conditions list (primary + backups)
    recommendations = []
    
    # Primary recommendation
    recommended_condition = synthesis.get("recommended_condition", {})
    if recommended_condition:
        primary_rec = _build_recommendation_from_llm(
            rank=1,
            condition_dict=recommended_condition,
            reaction_smiles=reaction_smiles,
            confidence=confidence_score,
            reasoning=synthesis.get("confidence_reasoning"),
            warnings=synthesis.get("warnings", []),
        )
        recommendations.append(primary_rec)
    
    # Backup conditions
    backup_conditions = synthesis.get("backup_conditions", [])
    for i, backup in enumerate(backup_conditions[:2], start=2):  # Max 2 backups
        backup_rec = _build_recommendation_from_llm(
            rank=i,
            condition_dict=backup,
            reaction_smiles=reaction_smiles,
            confidence=confidence_score - (i-1) * 0.1,
            reasoning=backup.get("when_to_use"),
            warnings=synthesis.get("warnings", []),
        )
        recommendations.append(backup_rec)
    
    # Build standard output
    output = {
        "meta": format_meta(
            model_type="LLM-synthesis",
            status="success",
            processing_time_ms=processing_time_ms or llm_metadata.get("processing_time_ms"),
            model_version="1.0.0",
        ),
        "input": format_input(
            reaction_smiles=reaction_smiles,
            requested_reaction_type=requested_type,
            detected_family=detected_type,
            detected_family_label=friendly_family_label(detected_type),
        ),
        "detection": format_detection(
            detected_type=detected_type,
            confidence=confidence_score,
            method="llm-multi-source",
            family_label=friendly_family_label(detected_type),
        ),
        "recommended_conditions": recommendations,
        "extras": {
            "llm_synthesis": synthesis,
            "sources_used": sources_used,
            "llm_metadata": llm_metadata,
        }
    }
    
    return output


def _build_recommendation_from_llm(
    rank: int,
    condition_dict: Dict[str, Any],
    reaction_smiles: str,
    confidence: float,
    reasoning: Optional[str] = None,
    warnings: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Build a single recommendation entry from LLM condition dictionary.
    
    Args:
        rank: Recommendation rank (1, 2, 3, ...)
        condition_dict: LLM recommended condition dict
        reaction_smiles: Reaction SMILES
        confidence: Confidence score
        reasoning: Reasoning text
        warnings: List of warnings
        
    Returns:
        Standard recommendation dictionary
    """
    from chemtools.output_formatter import (
        normalize_chemical_entry,
        starting_material_entries,
        normalize_rule_string_value,
    )
    
    # Start with reactants
    chemicals = starting_material_entries(reaction_smiles)
    
    # Extract catalyst system from LLM output
    catalyst_name = condition_dict.get("catalyst")
    ligand_name = condition_dict.get("ligand")
    solvent_name = condition_dict.get("solvent")
    base_name = condition_dict.get("base")
    additive_name = condition_dict.get("additive")
    temperature = condition_dict.get("temperature")
    
    # Add chemicals with enrichment
    if catalyst_name and catalyst_name not in ["None", "N/A", None]:
        chemicals.append(normalize_rule_string_value(catalyst_name, "metal_catalyst"))
    
    if ligand_name and ligand_name not in ["None", "N/A", None, "None (pre-formed catalyst)", "pre-complexed"]:
        chemicals.append(normalize_rule_string_value(ligand_name, "ligand"))
    
    if base_name and base_name not in ["None", "N/A", None]:
        chemicals.append(normalize_rule_string_value(base_name, "base", amount="200%"))
    
    if solvent_name and solvent_name not in ["None", "N/A", None]:
        chemicals.append(normalize_rule_string_value(solvent_name, "solvent"))
    
    if additive_name and additive_name not in ["None", "N/A", None]:
        chemicals.append(normalize_rule_string_value(additive_name, "additive"))
    
    # Parse temperature
    temp_value = None
    temp_range = None
    if temperature:
        # Handle formats like "80°C", "80-90°C", "80-100 °C"
        import re
        temp_match = re.search(r'(\d+)(?:\s*-\s*(\d+))?\s*°?C?', str(temperature))
        if temp_match:
            temp_low = float(temp_match.group(1))
            temp_high = float(temp_match.group(2)) if temp_match.group(2) else temp_low
            temp_value = (temp_low + temp_high) / 2
            temp_range = [temp_low, temp_high]
    
    # Build conditions block
    conditions = {}
    if temp_value:
        conditions["temperature"] = temp_range if temp_range else temp_value
    
    # Default time for synthesis
    conditions["time"] = [6.0, 24.0]  # Typical reaction time range
    conditions["atmosphere"] = None  # Will be inferred from chemistry
    
    # Build recommendation
    recommendation = {
        "rank": rank,
        "chemicals": chemicals,
        "conditions": conditions,
        "confidence": round(confidence, 4),
    }
    
    # Add reasoning section for robotic logs
    if reasoning or warnings:
        recommendation["summary"] = {
            "rationale": condition_dict.get("rationale") or reasoning or "",
            "warnings": warnings or [],
        }
    
    # Add source info
    recommendation["source"] = {
        "method": "llm-multi-source-synthesis",
        "confidence_level": "high" if confidence > 0.85 else "medium" if confidence > 0.65 else "low",
    }
    
    return recommendation


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
