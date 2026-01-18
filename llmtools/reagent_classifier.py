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

from chemtools.reagent import RegistryStore, load_families_registry_entries


# Valid roles - must match ROLE_CONFIG in reagent_taxonomy_qt.py
VALID_ROLES = {
    "ligand",
    "metal_catalyst",
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
    """Format family list for LLM prompt with metadata from unified taxonomy."""
    lines = []
    for fam in families:
        fam_id = fam.get("family", "")
        definition = fam.get("definition", "")
        keywords = fam.get("keywords", [])
        examples_pos = fam.get("examples_pos", [])
        notes = fam.get("notes", "")
        
        # Build rich description
        desc_parts = [f"- **{fam_id}**: {definition}"]
        
        # Add keywords if available
        if keywords:
            kw_str = ", ".join(keywords[:8])  # Increased from 5 to 8
            desc_parts.append(f"\n  Keywords: {kw_str}")
        
        # Add example CAS numbers if available
        if examples_pos:
            ex_str = ", ".join(examples_pos[:3])
            desc_parts.append(f"\n  Example CAS: {ex_str}")
        
        # Add notes if helpful
        if notes and len(notes) < 100:
            desc_parts.append(f"\n  Note: {notes}")
        
        lines.append("".join(desc_parts))
    
    return "\n\n".join(lines) if lines else "No families available"


def _load_schema_for_role(role: str, registry_dir: Path) -> Dict[str, Any]:
    """Return field definitions for a specific role."""
    default_schemas = {
        "base": {
            "basicity": ["weak", "moderate", "strong", "superbase"],
            "nucleophilicity": ["weak", "moderate", "strong"],
            "sterics": ["unhindered", "moderate", "hindered"]
        },
        "metal_precursor": {
            "metal": "Element symbol (Pd, Ni, Cu, Co, Fe, Ru, Rh, Ir, Au, Ag, Pt, Zn, etc.)",
            "oxidation_states": "List of integers [0], [1], [2], [0,2], etc."
        },
        "metal_catalyst": {
            "metal": "Element symbol (Pd, Ni, Cu, Co, Fe, Ru, Rh, Ir, Au, Ag, Pt, Zn)",
            "oxidation_states": "List of integers",
            "ligand_type": "Optional descriptor for pre-ligated complexes (phosphine, NHC, bipyridine, phenanthroline, pincer, carboxylate, other). Omit for simple salts."
        },
        "solvent": {
            "proticity": ["aprotic", "protic", "acidic", "basic"],
            "polarity": ["very_low", "low", "medium", "high", "very_high"],
            "coordination": ["non", "weak", "moderate", "strong"]
        },
        "ligand": {
            "donors": "List of donor atoms: ['P'], ['N','N'], ['P','P'], ['N','P'], etc.",
            "denticity": "Integer (1=monodentate, 2=bidentate, 3=tridentate, etc.)"
        },
        "oxidant": {
            "strength_band": ["weak", "moderate", "strong", "very_strong"]
        },
        "reductant": {
            "strength_band": ["weak", "moderate", "strong", "very_strong"]
        },
        "condensation_agent": {
            "strength_band": ["weak", "moderate", "strong"]
        },
        "acid": {
            "acidity": ["very_weak", "weak", "medium", "strong", "very_strong"]
        },
        "organo_catalyst": {
            "activation_mode": ["enamine", "iminium", "br酶nsted_acid", "hydrogen_bonding", "radical", "other"],
            "chirality": ["achiral", "racemic", "single_enantiomer"]
        },
        "enzyme": {
            "source": ["isolated", "whole_cell", "engineered", "other"],
            "cofactor_requirement": ["none", "metal", "organic", "multiple"]
        },
        "additive": {},  # No required fields beyond families
        "other_reagent": {}  # No required fields beyond families
    }

    return default_schemas.get(role, {})


def _format_fields_schema(role: str, registry_dir: Path) -> str:
    """Format field schema with allowed values from actual reagent schema."""
    schema = _load_schema_for_role(role, registry_dir)
    
    if not schema:
        return "No required fields beyond families"
    
    lines = []
    for field_name, allowed_values in schema.items():
        if isinstance(allowed_values, list):
            # Enumerated values
            values_str = " | ".join(f'"{v}"' for v in allowed_values)
            lines.append(f'- **{field_name}**: {values_str}')
        elif isinstance(allowed_values, str):
            # Description or free-form field
            if "|" in allowed_values:
                # Pipe-separated options in string format
                lines.append(f'- **{field_name}**: {allowed_values}')
            else:
                # Free-form description
                lines.append(f'- **{field_name}**: {allowed_values}')
        else:
            # Unknown format
            lines.append(f'- **{field_name}**: {str(allowed_values)}')
    
    return "\n".join(lines) if lines else "No required fields beyond families"


def _get_example_entries(role: str, family: Optional[str], registry_dir: Path, limit: int = 2) -> str:
    """Get example entries from existing registry."""
    try:
        store = RegistryStore(registry_dir)
        entries = store.role_entries.get(role, [])
        if family:
            entries = [
                entry
                for entry in entries
                if family in entry.get("roles", {}).get(role, {}).get("families", [])
            ]
        if not entries:
            return "No examples available"

        examples = []
        for entry in entries[:limit]:
            role_payload = entry.get("roles", {}).get(role, {})
            ex = {
                "name": entry.get("name"),
                "cas": entry.get("cas"),
                "families": role_payload.get("families", []),
            }
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
    # Load families for this role from unified taxonomy
    families_entries: List[Dict[str, Any]] = []
    try:
        families_entries = load_families_registry_entries()
    except FileNotFoundError as exc:
        return {
            "status": "error",
            "error": f"Families registry unavailable: {exc}",
        }

    role_families = [f for f in families_entries if f.get("role") == role]
    
    if not role_families:
        return {
            "status": "error",
            "error": f"No families found for role '{role}'",
        }
    
    # Format prompt components
    families_desc = _format_families_description(role_families)
    
    # Get field schema from defaults
    fields_schema = _format_fields_schema(role, registry_dir)
    
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


