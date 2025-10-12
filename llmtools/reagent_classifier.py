"""
LLM-Based Reagent Classifier
============================

Pure LLM workflow for reagent taxonomy classification.
Replaces deterministic token matching with intelligent chemistry-based classification.

Workflow:
1. resolve_identity() - Get chemical info from PubChem/CAS
2. classify_role() - LLM determines reagent role  
3. assign_fields() - LLM picks family and assigns role-specific fields
4. verify_entry() - LLM quality control check

All functions return structured dicts with confidence scores.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from llmtools.clients import LLMClient
from llmtools.prompts import (
    REAGENT_ENTRY_VERIFICATION,
    REAGENT_FIELD_ASSIGNMENT,
    REAGENT_ROLE_CLASSIFICATION,
)
from llmtools.reagent_review import _strip_markdown_fences


# Valid roles - must match ROLE_CONFIG in reagent_taxonomy_qt.py
VALID_ROLES = {
    "ligand",
    "metal_precursor",
    "preformed_metal_catalyst",
    "base",
    "acid",
    "condensation_agent",
    "oxidant",
    "reductant",
    "additive",
    "solvent",
    "organo_catalyst",
    "enzyme",
    "other_reagent",
}


def classify_role(
    identity: Dict[str, Any],
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 300,
) -> Dict[str, Any]:
    """
    Use LLM to classify reagent role from chemical identity.
    
    Args:
        identity: Chemical identity dict with name, cas, smiles, etc.
        llm_client: Configured LLM client
        temperature: LLM temperature (0.0 = deterministic)
        max_tokens: Max response tokens
        
    Returns:
        {
            "role": "base",
            "confidence": 0.95,
            "reasoning": "Tertiary amine with basicity",
            "status": "success" | "error",
            "error": "..." (if status=error)
        }
    """
    # Format synonyms list
    synonyms = identity.get("synonyms", [])
    synonyms_str = ", ".join(synonyms[:5]) if synonyms else "None"
    
    # Build prompt
    prompt = REAGENT_ROLE_CLASSIFICATION.format(
        name=identity.get("name", "Unknown"),
        cas=identity.get("cas", "Unknown"),
        smiles=identity.get("smiles", "Unknown"),
        molecular_formula=identity.get("molecular_formula", "Unknown"),
        synonyms=synonyms_str,
    )
    
    # Call LLM
    try:
        response = llm_client.chat(
            prompt=prompt,
            temperature=temperature,
            max_tokens=max_tokens,
        )
    except Exception as exc:
        return {
            "status": "error",
            "error": f"LLM call failed: {exc}",
            "role": None,
            "confidence": 0.0,
        }
    
    # Parse response
    raw_content = (response.content or "").strip()
    cleaned_content = _strip_markdown_fences(raw_content)
    
    try:
        parsed = json.loads(cleaned_content)
    except json.JSONDecodeError as exc:
        return {
            "status": "error",
            "error": f"JSON parse failed: {exc}",
            "raw_response": raw_content,
            "role": None,
            "confidence": 0.0,
        }
    
    # Validate role
    role = parsed.get("role", "").strip().lower()
    if role not in VALID_ROLES:
        return {
            "status": "error",
            "error": f"Invalid role '{role}'. Must be one of: {', '.join(sorted(VALID_ROLES))}",
            "raw_response": raw_content,
            "role": None,
            "confidence": 0.0,
        }
    
    # Success
    return {
        "status": "success",
        "role": role,
        "confidence": float(parsed.get("confidence", 0.0)),
        "reasoning": parsed.get("reasoning", ""),
        "model": response.model,
        "tokens": response.total_tokens,
        "latency_ms": response.latency_ms,
    }


def _format_families_description(families: List[Dict[str, Any]]) -> str:
    """Format family list for LLM prompt."""
    lines = []
    for fam in families:
        fam_id = fam.get("family", "")
        definition = fam.get("definition", "")
        keywords = fam.get("keywords", [])
        kw_str = f" (keywords: {', '.join(keywords[:5])})" if keywords else ""
        lines.append(f"- **{fam_id}**: {definition}{kw_str}")
    return "\n".join(lines) if lines else "No families available"


def _format_fields_schema(role: str, field_names: List[str], registry_dir: Path) -> str:
    """Format field schema with allowed values."""
    # This would ideally pull from actual schema, but for now we'll use hardcoded
    
    field_schemas = {
        "base": {
            "basicity": ["weak", "moderate", "strong", "superbase"],
            "nucleophilicity": ["weak", "moderate", "strong"],
            "sterics": ["unhindered", "moderate", "hindered"],
        },
        "metal_precursor": {
            "metal": "Element symbol (e.g., Pd, Ni, Cu, Fe, Ru, Ir, Rh, Pt, Au, Ag)",
            "oxidation_states": "List of integers (e.g., [0], [2], [0, 2])",
        },
        "preformed_metal_catalyst": {
            "metal": "Element symbol (e.g., Pd, Ni, Cu)",
            "oxidation_states": "List of integers",
            "ligand_type": "Optional string describing ligand coordination",
        },
        "solvent": {
            "proticity": ["aprotic", "protic"],
            "polarity": ["polar", "nonpolar", "intermediate"],
            "coordination": ["coordinating", "weakly_coordinating", "non_coordinating"],
        },
        "ligand": {
            "donors": "List of donor atoms (e.g., ['P'], ['N', 'N'], ['P', 'P'])",
            "denticity": "Integer (1, 2, 3, 4, etc.)",
        },
        "oxidant": {
            "strength_band": ["weak", "moderate", "strong"],
        },
        "reductant": {
            "strength_band": ["weak", "moderate", "strong"],
        },
        "condensation_agent": {
            "strength_band": ["weak", "moderate", "strong"],
        },
        "acid": {
            "acidity": ["weak", "moderate", "strong", "superacid"],
        },
        "organo_catalyst": {
            "activation_mode": "String (e.g., 'lewis_base', 'phase_transfer', 'hydrogen_bonding')",
            "chirality": ["achiral", "chiral"],
        },
        "enzyme": {
            "source": "String (e.g., 'bacterial', 'fungal', 'mammalian', 'plant')",
            "cofactor_requirement": "String or null",
        },
    }
    
    schema = field_schemas.get(role, {})
    if not schema:
        return "No required fields"
    
    lines = []
    for field_name, allowed_values in schema.items():
        if isinstance(allowed_values, list):
            values_str = " | ".join(f'"{v}"' for v in allowed_values)
            lines.append(f'- **{field_name}**: {values_str}')
        else:
            lines.append(f'- **{field_name}**: {allowed_values}')
    
    return "\n".join(lines)


def _get_example_entries(role: str, family: Optional[str], registry_dir: Path, limit: int = 2) -> str:
    """Get example entries from existing registry."""
    try:
        role_file = registry_dir / f"{role}.json"
        if not role_file.exists():
            return "No examples available"
        
        entries = json.loads(role_file.read_text(encoding="utf-8"))
        if not isinstance(entries, list):
            return "No examples available"
        
        # Filter by family if specified
        if family:
            entries = [e for e in entries if family in e.get("roles", {}).get(role, {}).get("families", [])]
        
        if not entries:
            return "No examples available"
        
        # Take first few examples
        examples = []
        for entry in entries[:limit]:
            role_payload = entry.get("roles", {}).get(role, {})
            ex = {
                "name": entry.get("name"),
                "cas": entry.get("cas"),
                "families": role_payload.get("families", []),
            }
            # Add role-specific fields
            for key, val in role_payload.items():
                if key != "families":
                    ex[key] = val
            examples.append(json.dumps(ex, indent=2))
        
        return "\n\n".join(examples)
        
    except Exception:
        return "No examples available"


def assign_fields(
    identity: Dict[str, Any],
    role: str,
    registry_dir: Path,
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 600,
) -> Dict[str, Any]:
    """
    Use LLM to assign family and populate role-specific fields.
    
    Args:
        identity: Chemical identity dict
        role: Reagent role (from classify_role)
        registry_dir: Path to registry directory (for families and examples)
        llm_client: Configured LLM client
        temperature: LLM temperature
        max_tokens: Max response tokens
        
    Returns:
        {
            "status": "success" | "error",
            "family": "tertiary_amines_aliphatic",
            "fields": {"basicity": "strong", ...},
            "abbreviations": ["TEA"],
            "additional_synonyms": [],
            "confidence": 0.92,
            "reasoning": "...",
            "error": "..." (if status=error)
        }
    """
    # Load families for this role
    families_path = registry_dir / "reagent_schema" / "families_registry.json"
    if not families_path.exists():
        return {
            "status": "error",
            "error": f"Families registry not found: {families_path}",
        }
    
    try:
        families_data = json.loads(families_path.read_text(encoding="utf-8"))
        all_families = families_data.get("entries", [])
        role_families = [f for f in all_families if f.get("role") == role]
    except Exception as exc:
        return {
            "status": "error",
            "error": f"Failed to load families: {exc}",
        }
    
    if not role_families:
        return {
            "status": "error",
            "error": f"No families found for role '{role}'",
        }
    
    # Format prompt components
    families_desc = _format_families_description(role_families)
    
    # Get expected fields for this role (from ROLE_PAYLOAD_FIELDS)
    # Hardcoded for now, ideally would import from reagent_taxonomy_qt.py
    field_names_map = {
        "base": ["basicity", "nucleophilicity", "sterics"],
        "metal_precursor": ["metal", "oxidation_states"],
        "preformed_metal_catalyst": ["metal", "oxidation_states", "ligand_type"],
        "solvent": ["proticity", "polarity", "coordination"],
        "ligand": ["donors", "denticity"],
        "oxidant": ["strength_band"],
        "reductant": ["strength_band"],
        "condensation_agent": ["strength_band"],
        "acid": ["acidity"],
        "organo_catalyst": ["activation_mode", "chirality"],
        "enzyme": ["source", "cofactor_requirement"],
    }
    field_names = field_names_map.get(role, [])
    fields_schema = _format_fields_schema(role, field_names, registry_dir)
    
    examples = _get_example_entries(role, None, registry_dir, limit=2)
    
    # Format synonyms
    synonyms = identity.get("synonyms", [])
    synonyms_str = ", ".join(synonyms[:5]) if synonyms else "None"
    
    # Build prompt
    prompt = REAGENT_FIELD_ASSIGNMENT.format(
        name=identity.get("name", "Unknown"),
        cas=identity.get("cas", "Unknown"),
        smiles=identity.get("smiles", "Unknown"),
        molecular_formula=identity.get("molecular_formula", "Unknown"),
        role=role,
        families_description=families_desc,
        fields_schema=fields_schema,
        examples=examples,
    )
    
    # Call LLM
    try:
        response = llm_client.chat(
            prompt=prompt,
            temperature=temperature,
            max_tokens=max_tokens,
        )
    except Exception as exc:
        return {
            "status": "error",
            "error": f"LLM call failed: {exc}",
        }
    
    # Parse response
    raw_content = (response.content or "").strip()
    cleaned_content = _strip_markdown_fences(raw_content)
    
    try:
        parsed = json.loads(cleaned_content)
    except json.JSONDecodeError as exc:
        return {
            "status": "error",
            "error": f"JSON parse failed: {exc}",
            "raw_response": raw_content,
        }
    
    # Validate family exists
    family = parsed.get("family", "").strip()
    if not any(f.get("family") == family for f in role_families):
        return {
            "status": "error",
            "error": f"Invalid family '{family}' for role '{role}'",
            "raw_response": raw_content,
        }
    
    # Success
    return {
        "status": "success",
        "family": family,
        "fields": parsed.get("fields", {}),
        "abbreviations": parsed.get("abbreviations", []),
        "additional_synonyms": parsed.get("additional_synonyms", []),
        "confidence": float(parsed.get("confidence", 0.0)),
        "reasoning": parsed.get("reasoning", ""),
        "model": response.model,
        "tokens": response.total_tokens,
        "latency_ms": response.latency_ms,
    }


def verify_entry(
    entry: Dict[str, Any],
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 500,
) -> Dict[str, Any]:
    """
    Use LLM to verify reagent entry for errors and inconsistencies.
    
    Args:
        entry: Complete reagent entry dict
        llm_client: Configured LLM client
        temperature: LLM temperature
        max_tokens: Max response tokens
        
    Returns:
        {
            "status": "success" | "error",
            "approved": true | false,
            "issues": [{"severity": "error|warning", "field": "...", "message": "..."}],
            "suggestions": ["...", "..."],
            "error": "..." (if status=error)
        }
    """
    # Format entry as JSON
    entry_json = json.dumps(entry, indent=2)
    
    # Build prompt
    prompt = REAGENT_ENTRY_VERIFICATION.format(entry_json=entry_json)
    
    # Call LLM
    try:
        response = llm_client.chat(
            prompt=prompt,
            temperature=temperature,
            max_tokens=max_tokens,
        )
    except Exception as exc:
        return {
            "status": "error",
            "error": f"LLM call failed: {exc}",
        }
    
    # Parse response
    raw_content = (response.content or "").strip()
    cleaned_content = _strip_markdown_fences(raw_content)
    
    try:
        parsed = json.loads(cleaned_content)
    except json.JSONDecodeError as exc:
        return {
            "status": "error",
            "error": f"JSON parse failed: {exc}",
            "raw_response": raw_content,
        }
    
    # Success
    return {
        "status": "success",
        "approved": bool(parsed.get("approved", False)),
        "issues": parsed.get("issues", []),
        "suggestions": parsed.get("suggestions", []),
        "model": response.model,
        "tokens": response.total_tokens,
        "latency_ms": response.latency_ms,
    }


__all__ = [
    "classify_role",
    "assign_fields",
    "verify_entry",
    "VALID_ROLES",
]
