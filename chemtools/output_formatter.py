"""
Enhanced output formatter for reaction condition recommendations.

Provides comprehensive structured output following the improved format specification.
"""

from __future__ import annotations

import copy
import datetime
import re
from typing import Any, Dict, List, Optional, Tuple
from . import reagent


def format_meta(
    model_type: str = "ML-precedent-knn",
    status: str = "success",
    processing_time_ms: Optional[float] = None,
    request_id: Optional[str] = None,
    schema_version: str = "2.0",
    model_version: str = "1.0.0",
) -> Dict[str, Any]:
    """Format metadata section."""
    meta = {
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "model": model_type,
        "model_version": model_version,
        "status": status,
        "version": schema_version,
    }
    
    if request_id:
        meta["request_id"] = request_id
    
    if processing_time_ms is not None:
        meta["processing_time_ms"] = round(processing_time_ms, 2)
    
    return meta


def format_input(
    reaction_smiles: str,
    requested_reaction_type: Optional[str] = None,
    detected_family: Optional[str] = None,
    constraints: Optional[Dict[str, Any]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Format input section."""
    input_data: Dict[str, Any] = {
        "reaction_smiles": reaction_smiles,
    }
    
    if detected_family:
        input_data["family"] = detected_family
    
    if requested_reaction_type:
        input_data["requested_reaction_type"] = requested_reaction_type
        # Provide an alias aligned with schema naming conventions.
        input_data["requested_family"] = requested_reaction_type
    
    if constraints:
        input_data["constraints"] = constraints
    
    if options:
        input_data["options"] = options
    
    return input_data


def format_detection(
    detected_type: str,
    confidence: Optional[float],
    method: str = "auto",
    alternatives: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Format detection section."""
    detection = {
        "family": detected_type,
        "detected_reaction_type": detected_type,
        "method": method,
    }
    
    if confidence is not None:
        detection["confidence"] = round(float(confidence), 4)
    
    if alternatives:
        detection["alternative_types"] = alternatives
    
    # Remove keys with None values to keep payload concise.
    detection = {key: value for key, value in detection.items() if value is not None}
    
    return detection


_REQUIRED_CHEM_FIELDS: Tuple[str, ...] = (
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "equivalents",
    "role",
)


def _normalize_chemical_entry(entry: Dict[str, Any], role_hint: Optional[str] = None) -> Dict[str, Any]:
    """Ensure chemical entries provide the required schema keys."""
    data: Dict[str, Any] = copy.deepcopy(entry or {})
    
    if role_hint and not data.get("role"):
        data["role"] = role_hint
    
    for key in _REQUIRED_CHEM_FIELDS:
        data.setdefault(key, None)
    
    return data


def _normalize_condition_value(value: Any, default_unit: Optional[str] = None) -> Any:
    """Normalize temperature/time values to a consistent structure."""
    if value is None:
        return None
    
    if isinstance(value, dict):
        normalized = copy.deepcopy(value)
        if "value" not in normalized:
            if "amount" in normalized:
                normalized["value"] = normalized.pop("amount")
            elif "quantity" in normalized:
                normalized["value"] = normalized.pop("quantity")
        if default_unit and "unit" not in normalized and "units" not in normalized:
            normalized.setdefault("unit", default_unit)
        return normalized
    
    if isinstance(value, (int, float)):
        result: Dict[str, Any] = {"value": value}
        if default_unit:
            result["unit"] = default_unit
        return result
    
    if isinstance(value, str):
        match = re.match(r"\s*([+-]?\d+(?:\.\d+)?)", value)
        if match:
            numeric = float(match.group(1))
            remainder = value[match.end():].strip()
            result = {"value": numeric}
            unit = remainder or default_unit
            if unit:
                result["unit"] = unit
            result["text"] = value
            return result
        return value
    
    return value


def _normalize_conditions_block(conditions: Any) -> Dict[str, Any]:
    """Normalize conditions dictionaries into the standard schema."""
    if not isinstance(conditions, dict):
        return {
            "temperature": None,
            "time": None,
            "atmosphere": None,
            "note": conditions,
        }
    
    normalized: Dict[str, Any] = {}
    
    # Temperature
    temp_value = (
        conditions.get("temperature")
        or conditions.get("temperature_C")
        or conditions.get("temp_C")
        or conditions.get("temp")
    )
    normalized["temperature"] = _normalize_condition_value(temp_value, "C")
    
    # Time
    time_value = (
        conditions.get("time")
        or conditions.get("time_h")
        or conditions.get("time_hr")
        or conditions.get("time_hours")
    )
    normalized["time"] = _normalize_condition_value(time_value, "h")
    
    # Atmosphere
    normalized["atmosphere"] = conditions.get("atmosphere") or conditions.get("environment")
    
    # Preserve any additional condition metadata
    extras = {
        key: copy.deepcopy(value)
        for key, value in conditions.items()
        if key not in {"temperature", "temperature_C", "temp_C", "temp", "time", "time_h", "time_hr", "time_hours", "atmosphere", "environment"}
    }
    if extras:
        normalized["extras"] = extras
    
    for key in ("temperature", "time", "atmosphere"):
        normalized.setdefault(key, None)
    
    return normalized


def _normalize_recommendation_entry(rec: Dict[str, Any], position: int) -> Dict[str, Any]:
    """Normalize a single recommendation entry to the standard schema."""
    source = rec or {}
    normalized: Dict[str, Any] = {}
    
    normalized["rank"] = source.get("rank") or position
    
    chemicals_source = source.get("chemicals")
    if not isinstance(chemicals_source, list) or not chemicals_source:
        chemicals_source = source.get("reagents") or []
    
    normalized_chemicals: List[Dict[str, Any]] = []
    for chem in chemicals_source:
        if isinstance(chem, dict):
            normalized_chemicals.append(_normalize_chemical_entry(chem))
        else:
            normalized_chemicals.append(_normalize_chemical_entry({"name": str(chem)}, None))
    normalized["chemicals"] = normalized_chemicals
    
    normalized["conditions"] = _normalize_conditions_block(source.get("conditions", {}))
    
    for key in ("confidence", "precedent_count", "summary", "component_scores", "reasons", "source", "reaction"):
        if key in source:
            normalized[key] = copy.deepcopy(source[key])
    
    extras_keys = set(source.keys()) - {
        "rank",
        "chemicals",
        "reagents",
        "conditions",
        "confidence",
        "precedent_count",
        "summary",
        "component_scores",
        "reasons",
        "source",
        "reaction",
    }
    if extras_keys:
        normalized["extras"] = {key: copy.deepcopy(source[key]) for key in extras_keys}
    
    return normalized


def _normalize_recommendations(recommendations_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Normalize a list of recommendation entries."""
    if not recommendations_data:
        return []
    return [
        _normalize_recommendation_entry(rec, position=index)
        for index, rec in enumerate(recommendations_data, start=1)
    ]


def _starting_material_entries(reaction_smiles: str) -> List[Dict[str, Any]]:
    """Build chemical entries for starting materials based on reaction SMILES."""
    if ">>" not in reaction_smiles:
        return []
    
    reactants = [frag for frag in reaction_smiles.split(">>", 1)[0].split(".") if frag]
    roles = ["electrophile", "nucleophile"]
    
    entries = []
    for idx, smiles in enumerate(reactants):
        role = roles[idx] if idx < len(roles) else "starting_material"
        entries.append(_normalize_chemical_entry({"smiles": smiles, "role": role}, role))
    
    return entries


def _parse_amount_to_equivalents(amount: str, role: str) -> Optional[float]:
    """Parse amount string to equivalents.
    
    Handles formats like:
    - "200%" or "200–250%" -> 2.0 or 2.0 (midpoint 2.25)
    - "0.5–1.5%" -> 0.005 or 0.005 (midpoint 0.01)
    - "2.0 eq" -> 2.0
    - "10–20%" for ligands -> 0.1 or 0.1 (midpoint 0.15)
    
    Args:
        amount: Amount string from rule file
        role: Reagent role (metal_source, ligand, base, solvent, etc.)
    
    Returns:
        Equivalents as float, or None if cannot parse
    """
    if not amount or not isinstance(amount, str):
        return None
    
    amount = amount.strip()
    
    # Check for explicit "eq" notation
    eq_match = re.search(r'([\d.]+)\s*eq', amount, re.IGNORECASE)
    if eq_match:
        try:
            return float(eq_match.group(1))
        except ValueError:
            return None
    
    # Check for percentage notation
    # Match patterns like "200%", "0.5–1.5%", "10-20%"
    pct_match = re.search(r'([\d.]+)(?:\s*[–\-]\s*([\d.]+))?\s*%', amount)
    if pct_match:
        try:
            low = float(pct_match.group(1))
            high = float(pct_match.group(2)) if pct_match.group(2) else low
            
            # For catalysts/ligands, % usually means mol%
            # For bases, % usually means % equivalents (e.g., 200% = 2.0 eq)
            if role in ['metal_precursor', 'metal_source', 'catalyst', 'ligand']:
                # Convert mol% to equivalents (e.g., 1.5% = 0.015 eq)
                midpoint = (low + high) / 2.0
                return midpoint / 100.0
            else:
                # For bases and other reagents, % is equivalents (e.g., 200% = 2.0 eq)
                midpoint = (low + high) / 2.0
                return midpoint / 100.0
        except (ValueError, TypeError):
            return None
    
    # Check for "M" notation (molarity - not equivalents)
    if 'M' in amount or 'm' in amount.lower():
        return None
    
    return None


def _normalize_rule_string_value(value: str, role: str, amount: Optional[str] = None) -> Dict[str, Any]:
    """Convert string-based rule entries into structured chemical dictionaries with database enrichment.
    
    Args:
        value: Reagent name/description from rule file
        role: Chemical role (metal_precursor, ligand, base, etc.)
        amount: Optional amount string to parse for equivalents
    """
    # Remove parenthetical notes to get clean name for lookup
    cleaned = re.sub(r"\s*\(.*?\)\s*$", "", value or "").strip()
    
    # Map output roles to database reagent types
    role_to_type_map = {
        "metal_precursor": "metal_precursor",
        "catalyst": "metal_precursor",
        "ligand": "ligand",
        "base": "base",
        "solvent": "solvent",
        "additive": "additive",
        "oxidant": "oxidant",
        "reductant": "reductant",
        "acid": "acid",
    }
    
    reagent_type = role_to_type_map.get(role, role)
    
    # Try to enrich with database information
    enriched_info = reagent.enrich_reagent_info(cleaned, reagent_type)
    
    # Build chemical entry
    chemical: Dict[str, Any] = {
        "name": enriched_info.get("name") if enriched_info.get("found") else cleaned or value or None,
        "role": role,
        "notes": value,
        "abbreviation": enriched_info.get("abbreviation") if enriched_info.get("found") else None,
        "cas": enriched_info.get("cas") if enriched_info.get("found") else None,
        "smiles": enriched_info.get("smiles") if enriched_info.get("found") else None,
    }
    
    # Parse equivalents from amount parameter (new standardized format)
    equivalents = None
    if amount:
        equivalents = _parse_amount_to_equivalents(amount, role)
    
    # Fallback: extract equivalents from the value string itself (legacy format)
    if equivalents is None:
        eq_match = re.search(r"([+-]?\d+(?:\.\d+)?)\s*eq", value or "", re.IGNORECASE)
        if eq_match:
            try:
                equivalents = float(eq_match.group(1))
            except ValueError:
                pass
    
    # Set equivalents if found
    if equivalents is not None:
        chemical["equivalents"] = equivalents
    else:
        chemical["equivalents"] = None
    
    # Extract mol% from the original value string for additional info
    mol_match = re.search(r"([+-]?\d+(?:\.\d+)?)\s*mol%?", value or "", re.IGNORECASE)
    if mol_match:
        try:
            chemical.setdefault("extras", {})["mol_percent"] = float(mol_match.group(1))
        except ValueError:
            pass
    
    return _normalize_chemical_entry(chemical, role)


def _convert_rule_match_to_recommendations(
    reaction_smiles: str,
    match_result: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """Convert SCDB match results into standardized recommendation entries.
    
    Supports both new standardized format (with 'reagents' array) and legacy format
    (with direct keys like 'pd_source', 'ligand', 'base', 'solvent').
    """
    conditions = match_result.get("conditions")
    if isinstance(conditions, dict):
        condition_entries = [conditions]
    elif isinstance(conditions, list):
        condition_entries = [cond for cond in conditions if isinstance(cond, dict)]
    else:
        condition_entries = []
    
    if not condition_entries:
        return []
    
    # Check if entry contains reagents data (new format is at match_result level via raw entry)
    entry_raw = match_result.get("trace", {}).get("selected", {})
    
    starting_materials = _starting_material_entries(reaction_smiles)
    recommendations: List[Dict[str, Any]] = []
    
    for idx, cond in enumerate(condition_entries, start=1):
        chemicals: List[Dict[str, Any]] = copy.deepcopy(starting_materials)
        
        # Role mapping for new standardized format
        role_mapping = {
            "metal_source": "metal_precursor",
            "metal_precursor": "metal_precursor",
            "catalyst": "metal_precursor",
            "ligand": "ligand",
            "base": "base",
            "solvent": "solvent",
            "additive": "additive",
            "partner": "partner",
            "boron_partner": "partner",
        }
        
        # Try new format first: check for reagents in the entry's raw data
        reagents_added = False
        if "reagents" in match_result:
            # Reagents at match result level (passed through from entry)
            reagents = match_result.get("reagents", [])
            if isinstance(reagents, list):
                for reagent in reagents:
                    if isinstance(reagent, dict):
                        name = reagent.get("name")
                        role = reagent.get("role", "reagent")
                        amount = reagent.get("amount")
                        
                        # Map role to output role
                        output_role = role_mapping.get(role, role)
                        
                        if name:
                            chemicals.append(_normalize_rule_string_value(str(name), output_role, amount=amount))
                            reagents_added = True
        
        # Fallback to legacy format if reagents not found in new format
        if not reagents_added:
            def add_entry(value: Any, role: str) -> None:
                if not value:
                    return
                if isinstance(value, dict):
                    chemicals.append(_normalize_chemical_entry(value, role))
                else:
                    chemicals.append(_normalize_rule_string_value(str(value), role))
            
            # Legacy format extraction
            add_entry(cond.get("pd_source") or cond.get("catalyst") or cond.get("metal_source"), "metal_precursor")
            add_entry(cond.get("ligand"), "ligand")
            add_entry(cond.get("base"), "base")
            add_entry(cond.get("solvent"), "solvent")
            add_entry(cond.get("additive"), "additive")
            add_entry(cond.get("boron_partner") or cond.get("partner"), "partner")
            add_entry(cond.get("condition_core"), "condition_core")
        
        recommendation = {
            "rank": idx,
            "chemicals": chemicals,
            "conditions": _normalize_conditions_block(cond),
        }
        
        loadings = cond.get("loadings")
        if isinstance(loadings, dict) and loadings:
            recommendation.setdefault("conditions", {}).setdefault("extras", {})["loadings"] = loadings
        
        source_meta = {
            "match_type": match_result.get("match_type"),
            "entry_id": match_result.get("entry_id"),
            "entry_name": match_result.get("entry_name"),
            "priority": match_result.get("priority"),
        }
        if any(value is not None for value in source_meta.values()):
            recommendation["source"] = source_meta
        
        recommendations.append(recommendation)
    
    return recommendations


def _build_standard_output(
    *,
    model_type: str,
    reaction_smiles: str,
    requested_type: Optional[str],
    detected_type: Optional[str],
    detection_confidence: Optional[float],
    detection_method: str,
    recommendations_data: List[Dict[str, Any]],
    processing_time_ms: Optional[float] = None,
    status: str = "success",
    extras: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Construct the canonical output structure shared across recommendation modes."""
    normalized_recommendations = _normalize_recommendations(recommendations_data)
    
    output: Dict[str, Any] = {
        "meta": format_meta(
            model_type=model_type,
            status=status,
            processing_time_ms=processing_time_ms,
        ),
        "input": format_input(
            reaction_smiles=reaction_smiles,
            requested_reaction_type=requested_type,
            detected_family=detected_type,
        ),
        "detection": format_detection(
            detected_type=detected_type or (requested_type or "unknown"),
            confidence=detection_confidence,
            method=detection_method,
        ),
        "recommended_conditions": normalized_recommendations,
    }
    
    if extras:
        output["extras"] = extras
    
    return output


def ensure_standard_output(
    response: Dict[str, Any],
    *,
    default_model: str,
    fallback_reaction_smiles: Optional[str] = None,
    fallback_requested_type: Optional[str] = None,
    extras: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Ensure an arbitrary response conforms to the canonical schema.
    
    If the response already satisfies the schema, it is returned unchanged.
    Otherwise we attempt to coerce it into the new structure.
    """
    required_keys = {"meta", "input", "detection", "recommended_conditions"}
    if isinstance(response, dict) and required_keys.issubset(response.keys()):
        return response
    
    response = response or {}
    input_section = response.get("input", {}) if isinstance(response, dict) else {}
    detection_section = response.get("detection", {}) if isinstance(response, dict) else {}
    
    reaction_smiles = input_section.get("reaction_smiles") or fallback_reaction_smiles or ""
    requested_type = (
        input_section.get("requested_family")
        or input_section.get("requested_reaction_type")
        or fallback_requested_type
    )
    detected_type = (
        detection_section.get("family")
        or detection_section.get("detected_reaction_type")
        or requested_type
    )
    detection_confidence = detection_section.get("confidence")
    detection_method = detection_section.get("method", "auto")
    processing_time_ms = (
        response.get("meta", {}).get("processing_time_ms")
        if isinstance(response.get("meta"), dict)
        else None
    )
    
    recomms = []
    if isinstance(response, dict):
        if isinstance(response.get("recommended_conditions"), list):
            recomms = response["recommended_conditions"]
        elif isinstance(response.get("recommendations"), list):
            recomms = response["recommendations"]
    
    combined_extras: Dict[str, Any] = {}
    if isinstance(response.get("extras"), dict):
        combined_extras.update(response["extras"])
    if extras:
        combined_extras.update(extras)
    combined_extras.setdefault("raw_response", copy.deepcopy(response))
    
    return _build_standard_output(
        model_type=default_model,
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type,
        detection_confidence=detection_confidence,
        detection_method=detection_method,
        recommendations_data=recomms,
        processing_time_ms=processing_time_ms,
        extras=combined_extras,
    )


def format_rule_match_result(
    reaction_smiles: str,
    match_result: Dict[str, Any],
    requested_type: Optional[str] = None,
    database_name: Optional[str] = None,
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """Format raw SCDB match results into the canonical output schema."""
    recommendations_data = _convert_rule_match_to_recommendations(reaction_smiles, match_result)
    detected_type = (
        match_result.get("reaction_type")
        or match_result.get("family")
        or requested_type
    )
    if not detected_type and database_name:
        detected_type = database_name
    extras = {
        "match": match_result,
    }
    model_label = f"Rule-based-{database_name}" if database_name else "Rule-based"
    return _build_standard_output(
        model_type=model_label,
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type,
        detection_confidence=1.0 if recommendations_data else None,
        detection_method="pattern-match",
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
        extras=extras,
    )


def format_fusion_output(
    reaction_smiles: str,
    requested_type: Optional[str],
    fusion_result: Dict[str, Any],
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """Format fusion recommendation results into the canonical schema."""
    formatted_section = dict(fusion_result.get("formatted") or {})
    detection_section = dict(formatted_section.get("detection") or fusion_result.get("detection") or {})
    detected_type = (
        detection_section.get("family")
        or detection_section.get("detected_reaction_type")
        or requested_type
    )
    confidence = detection_section.get("confidence")
    if not isinstance(confidence, (int, float)):
        confidence = None
    method = detection_section.get("method", "fusion")
    
    recommendations_data = (
        formatted_section.get("recommended_conditions")
        or fusion_result.get("recommended_conditions")
        or []
    )
    
    extras = {
        "fusion_meta": fusion_result.get("fusion_meta"),
        "raw_response": fusion_result,
    }
    
    meta_section = fusion_result.get("meta")
    meta_processing = None
    if isinstance(meta_section, dict):
        meta_processing = meta_section.get("processing_time_ms")
    if processing_time_ms is None:
        processing_time_ms = meta_processing
    
    return _build_standard_output(
        model_type="Fusion-hybrid",
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type,
        detection_confidence=confidence,
        detection_method=method,
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
        extras={key: value for key, value in extras.items() if value is not None},
    )



def enrich_reagent(
    name: str,
    reagent_type: str,
    role: str,
    equivalents: Optional[float] = None,
    equiv_range: Optional[Tuple[float, float]] = None,
    smiles: Optional[str] = None,
    reagent_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Enrich a single reagent with database information.
    
    Args:
        name: Reagent name or abbreviation
        reagent_type: Type for database lookup (ligand, base, solvent, etc.)
        role: Reagent role (catalyst, ligand, base, solvent, etc.)
        equivalents: Stoichiometric equivalents
        equiv_range: Range of equivalents [min, max]
        smiles: SMILES string (for starting materials)
        reagent_id: Optional unique identifier
    
    Returns:
        Enriched reagent dictionary
    """
    # Try to enrich from database
    enriched = {
        "id": reagent_id or f"{role.upper()[:3]}_1",
        "role": role,
    }
    
    # For starting materials, use SMILES
    if role in ["electrophile", "nucleophile", "starting_material"] and smiles:
        enriched.update({
            "name": None,
            "abbreviation": None,
            "cas": None,
            "smiles": smiles,
            "inchi_key": None,
        })
    else:
        # Look up in database
        db_info = reagent.enrich_reagent_info(name, reagent_type)
        
        if db_info.get("found"):
            enriched.update({
                "name": db_info.get("name"),
                "abbreviation": db_info.get("abbreviation"),
                "cas": db_info.get("cas"),
                "smiles": db_info.get("smiles"),
                "inchi_key": db_info.get("inchi_key"),
            })
            
            # Add category for catalysts
            if role == "catalyst":
                enriched["category"] = reagent_type
            
            # Add molecular weight if available
            # TODO: Calculate from SMILES or look up
            
            # Add properties from database
            if db_info.get("roles") and reagent_type in db_info["roles"]:
                role_props = db_info["roles"][reagent_type]
                enriched["properties"] = role_props
        else:
            # Not found in database - use provided info
            enriched.update({
                "name": name,
                "abbreviation": None,
                "cas": None,
                "smiles": smiles,
                "inchi_key": None,
            })
    
    # Add equivalents
    if equivalents is not None:
        if equiv_range:
            enriched["equivalents"] = {
                "value": equivalents,
                "range": list(equiv_range),
                "unit": "eq"
            }
        else:
            enriched["equivalents"] = {
                "value": equivalents,
                "range": [equivalents, equivalents],
                "unit": "eq"
            }
    
    # Add loading for catalysts and ligands
    if role in ["catalyst", "ligand"] and equivalents is not None:
        loading_pct = equivalents * 100
        enriched["loading"] = {
            "value": round(loading_pct, 1),
            "range": [round(equiv_range[0] * 100, 1), round(equiv_range[1] * 100, 1)] if equiv_range else [loading_pct, loading_pct],
            "unit": "mol%"
        }
    
    return enriched


def format_conditions(
    temperature: Optional[float] = None,
    temp_range: Optional[Tuple[float, float]] = None,
    time_hours: Optional[float] = None,
    time_range: Optional[Tuple[float, float]] = None,
    atmosphere: Optional[str] = None,
    pressure: Optional[float] = None,
) -> Dict[str, Any]:
    """Format reaction conditions."""
    conditions = {}
    
    if temperature is not None:
        conditions["temperature"] = {
            "value": temperature,
            "range": list(temp_range) if temp_range else [temperature, temperature],
            "unit": "°C"
        }
    
    if time_hours is not None:
        conditions["time"] = {
            "value": time_hours,
            "range": list(time_range) if time_range else [time_hours, time_hours],
            "unit": "hours"
        }
    
    if atmosphere:
        conditions["atmosphere"] = {
            "type": "inert",
            "gas": atmosphere,
            "requirement": "anhydrous" if "N" in atmosphere or "Ar" in atmosphere else "ambient"
        }
    
    if pressure:
        conditions["pressure"] = {
            "value": pressure,
            "unit": "atm"
        }
    
    return conditions


def format_recommendation(
    rank: int,
    confidence: float,
    support: int,
    reaction_smiles: str,
    reagents: List[Dict[str, Any]],
    conditions: Dict[str, Any],
    similarity_score: Optional[float] = None,
    expected_yield: Optional[float] = None,
    yield_range: Optional[Tuple[float, float]] = None,
    precedents: Optional[List[Dict[str, Any]]] = None,
    reasoning: Optional[Dict[str, Any]] = None,
    notes: Optional[str] = None,
    tags: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Format a single recommendation with enhanced metadata."""
    recommendation = {
        "rank": rank,
        "confidence": round(confidence, 4),
        "support": support,
        "reaction": {
            "smiles": reaction_smiles,
        },
        "reagents": reagents,
        "conditions": conditions,
    }
    
    if similarity_score is not None:
        recommendation["similarity_score"] = round(similarity_score, 4)
    
    if expected_yield is not None:
        recommendation["expected_outcome"] = {
            "yield": {
                "typical": round(expected_yield, 1),
                "range": [round(yield_range[0], 1), round(yield_range[1], 1)] if yield_range else [expected_yield, expected_yield],
                "unit": "%"
            }
        }
    
    if precedents:
        recommendation["precedents"] = {
            "count": support,
            "top_similar": precedents[:5]  # Top 5 most similar
        }
    
    # Add reasoning section
    if reasoning:
        recommendation["reasoning"] = reasoning
    elif support > 0:
        # Auto-generate basic reasoning if not provided
        recommendation["reasoning"] = {
            "method": "precedent-based",
            "basis": f"Based on {support} similar precedent reactions",
            "confidence_factors": [
                f"High structural similarity (score: {similarity_score:.3f})" if similarity_score and similarity_score > 0.8 else None,
                f"{support} successful precedent reactions" if support > 10 else None,
                "Standard conditions for this reaction type",
            ],
        }
        # Remove None values
        recommendation["reasoning"]["confidence_factors"] = [
            f for f in recommendation["reasoning"]["confidence_factors"] if f is not None
        ]
    
    # Add optional notes
    if notes:
        recommendation["notes"] = notes
    
    # Add optional tags
    if tags:
        recommendation["tags"] = tags
    
    return recommendation


def format_ml_output(
    reaction_smiles: str,
    requested_type: Optional[str],
    detected_type: str,
    detection_confidence: float,
    recommendations_data: List[Dict[str, Any]],
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Format ML recommendation output in the improved format.
    
    Args:
        reaction_smiles: Input reaction SMILES
        requested_type: User-requested reaction type
        detected_type: Auto-detected reaction type
        detection_confidence: Detection confidence score
        recommendations_data: List of recommendation data dicts
        processing_time_ms: Processing time in milliseconds
    
    Returns:
        Complete formatted output dictionary
    """
    return _build_standard_output(
        model_type="ML-precedent-knn",
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type,
        detection_confidence=detection_confidence,
        detection_method="rxn-insight-ml",
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
    )


def format_rule_output(
    reaction_smiles: str,
    requested_type: Optional[str],
    detected_type: Optional[str],
    recommendations_data: Any,
    database_name: str,
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Format rule-based recommendation output in the improved format.
    
    Args:
        reaction_smiles: Input reaction SMILES
        requested_type: User-requested reaction type
        detected_type: Auto-detected reaction type
        recommendations_data: List of recommendation data dicts
        database_name: Name of rule database used
        processing_time_ms: Processing time in milliseconds
    
    Returns:
        Complete formatted output dictionary
    """
    if isinstance(recommendations_data, dict):
        return format_rule_match_result(
            reaction_smiles=reaction_smiles,
            match_result=recommendations_data,
            requested_type=requested_type,
            database_name=database_name,
            processing_time_ms=processing_time_ms,
        )
    
    return _build_standard_output(
        model_type=f"Rule-based-{database_name}",
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type or requested_type,
        detection_confidence=1.0 if recommendations_data else None,
        detection_method="pattern-match",
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
    )


def parse_condition_string(condition_str: str) -> Tuple[str, Optional[float]]:
    """
    Parse condition strings like 'K2CO3 2.0 eq' into name and equivalents.
    
    Args:
        condition_str: Condition string
    
    Returns:
        Tuple of (name, equivalents)
    """
    import re
    
    # Try to extract equivalents
    match = re.search(r'(\d+\.?\d*)\s*eq', condition_str, re.IGNORECASE)
    
    if match:
        equiv = float(match.group(1))
        # Remove equiv part to get name
        name = re.sub(r'\s*\d+\.?\d*\s*eq.*$', '', condition_str, flags=re.IGNORECASE).strip()
        return name, equiv
    
    return condition_str.strip(), None


def expand_rule_conditions_to_recommendations(
    reaction_smiles: str,
    conditions_dict: Dict[str, Any],
    num_recommendations: int = 3,
) -> List[Dict[str, Any]]:
    """
    Expand rule-based conditions (which have multiple options) into separate recommendations.
    
    Args:
        reaction_smiles: Reaction SMILES
        conditions_dict: Conditions dictionary from rule matching
        num_recommendations: Number of recommendations to generate
    
    Returns:
        List of recommendation dictionaries
    """
    recommendations = []
    
    # Extract reactants from SMILES
    parts = reaction_smiles.split(">>")
    if len(parts) == 2:
        reactants = parts[0].split(".")
    else:
        reactants = []
    
    # Generate variations by picking different options
    # Strategy: First recommendation uses first option, second uses second option, etc.
    
    for i in range(num_recommendations):
        reagents = []
        reagent_counter = {"cat": 1, "lig": 1, "bas": 1, "sol": 1, "oxi": 1, "red": 1}
        
        # Add starting materials
        for j, reactant_smiles in enumerate(reactants[:2], 1):
            role = "electrophile" if j == 1 else "nucleophile"
            reagents.append({
                "id": f"SM{j}",
                "name": None,
                "abbreviation": None,
                "cas": None,
                "smiles": reactant_smiles,
                "inchi_key": None,
                "role": role,
                "equivalents": {
                    "value": 1.0 if j == 1 else 1.2,
                    "range": [1.0, 1.0] if j == 1 else [1.1, 1.5],
                    "unit": "eq"
                }
            })
        
        # Map condition keys to reagent info
        condition_map = {
            "pd_source": ("catalyst", "metal_precursor"),
            "cu_source": ("catalyst", "metal_precursor"),
            "ni_source": ("catalyst", "metal_precursor"),
            "metal_source": ("catalyst", "metal_precursor"),
            "catalyst": ("catalyst", "metal_precursor"),
            "ligand": ("ligand", "ligand"),
            "ligands": ("ligand", "ligand"),
            "base": ("base", "base"),
            "solvent": ("solvent", "solvent"),
            "solvents": ("solvent", "solvent"),
            "oxidant": ("oxidant", "oxidant"),
            "reductant": ("reductant", "reductant"),
        }
        
        # Process conditions
        for key, value in conditions_dict.items():
            if key in ["temperature_C", "time_h", "loadings"]:
                continue  # Handle separately
            
            if key not in condition_map:
                continue
            
            role, reagent_type = condition_map[key]
            counter_key = role[:3].lower()
            
            # Handle list of options - pick based on recommendation index
            if isinstance(value, list) and len(value) > 0:
                # Cycle through options
                option_idx = i % len(value)
                selected_value = value[option_idx]
                
                # Parse name and equivalents
                name, equiv = parse_condition_string(str(selected_value))
                
                # Enrich reagent
                enriched = enrich_reagent(
                    name=name,
                    reagent_type=reagent_type,
                    role=role,
                    equivalents=equiv if equiv else (0.1 if role in ["catalyst", "ligand"] else 2.0),
                    equiv_range=(0.05, 0.2) if role in ["catalyst", "ligand"] else (1.5, 2.5) if equiv else None,
                    reagent_id=f"{role.upper()[:3]}{reagent_counter[counter_key]}"
                )
                
                reagent_counter[counter_key] += 1
                reagents.append(enriched)
            
            elif isinstance(value, str):
                name, equiv = parse_condition_string(value)
                
                enriched = enrich_reagent(
                    name=name,
                    reagent_type=reagent_type,
                    role=role,
                    equivalents=equiv if equiv else (0.1 if role in ["catalyst", "ligand"] else 2.0),
                    equiv_range=(0.05, 0.2) if role in ["catalyst", "ligand"] else (1.5, 2.5) if equiv else None,
                    reagent_id=f"{role.upper()[:3]}{reagent_counter[counter_key]}"
                )
                
                reagent_counter[counter_key] += 1
                reagents.append(enriched)
        
        # Format conditions
        temp_ranges = conditions_dict.get("temperature_C", [])
        time_ranges = conditions_dict.get("time_h", [])
        
        if isinstance(temp_ranges, list) and len(temp_ranges) >= 2:
            temp_value = (temp_ranges[0] + temp_ranges[1]) / 2
            temp_range = (temp_ranges[0], temp_ranges[1])
        else:
            temp_value = 110
            temp_range = (80, 140)
        
        if isinstance(time_ranges, list) and len(time_ranges) >= 2:
            time_value = (time_ranges[0] + time_ranges[1]) / 2
            time_range = (time_ranges[0], time_ranges[1])
        else:
            time_value = 12
            time_range = (6, 24)
        
        formatted_conditions = format_conditions(
            temperature=temp_value,
            temp_range=temp_range,
            time_hours=time_value,
            time_range=time_range,
            atmosphere="N₂ or Ar",
            pressure=1.0,
        )
        
        # Create recommendation
        recommendation = {
            "rank": i + 1,
            "confidence": 0.95 - (i * 0.05),  # Slightly lower confidence for variants
            "support": 1,  # Rule-based support is 1 (exact match)
            "reaction": {
                "smiles": reaction_smiles,
            },
            "reagents": reagents,
            "conditions": formatted_conditions,
        }
        
        recommendations.append(recommendation)
    
    return recommendations


def create_standard_output(
    reaction_smiles: str,
    detected_type: str,
    recommendations: List[Dict[str, Any]],
    requested_type: Optional[str] = None,
    model_type: str = "ML-precedent-knn",
    detection_method: str = "auto",
    detection_confidence: float = 0.95,
    processing_time_ms: Optional[float] = None,
    request_id: Optional[str] = None,
    alternatives: Optional[List[Dict[str, Any]]] = None,
    constraints: Optional[Dict[str, Any]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create standardized JSON output for both ML and rule-based recommendations.
    
    This is the main function to use for generating standard outputs that match
    the format_of_recommended_reactions.json specification.
    
    Args:
        reaction_smiles: Input reaction SMILES
        detected_type: Detected reaction type
        recommendations: List of formatted recommendation dictionaries
        requested_type: User-requested reaction type (optional)
        model_type: Model type (ML-precedent-knn, Rule-based-SCDB, etc.)
        detection_method: Detection method (auto, rule-based, ml, pattern-match)
        detection_confidence: Detection confidence (0-1)
        processing_time_ms: Processing time in milliseconds
        request_id: Optional request ID for tracking
        alternatives: Alternative reaction type detections
        constraints: Applied constraints (inventory, blacklist, etc.)
        options: Request options (k, max_variants, etc.)
    
    Returns:
        Complete standardized output dictionary
    
    Example:
        >>> output = create_standard_output(
        ...     reaction_smiles="Brc1ccc(OC)cc1.Nc1ccccc1>>COc1ccc(Nc2ccccc2)cc1",
        ...     detected_type="C_N_Coupling_Cu",
        ...     recommendations=[rec1, rec2, rec3],
        ...     model_type="ML-precedent-knn",
        ...     detection_confidence=0.95,
        ... )
    """
    
    # Build metadata
    meta = {
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "model": model_type,
        "model_version": "1.0.0",
        "status": "success",
    }
    
    if request_id:
        meta["request_id"] = request_id
    
    if processing_time_ms is not None:
        meta["processing_time_ms"] = round(processing_time_ms, 2)
    
    # Build input section
    input_data = {
        "reaction_smiles": reaction_smiles,
    }
    
    if requested_type:
        input_data["requested_reaction_type"] = requested_type
    
    if constraints:
        input_data["constraints"] = constraints
    
    if options:
        input_data["options"] = options
    
    # Build detection section
    detection = {
        "detected_reaction_type": detected_type,
        "confidence": round(detection_confidence, 4),
        "method": detection_method,
    }
    
    if alternatives:
        detection["alternative_types"] = alternatives
    
    # Build complete output
    output = {
        "meta": meta,
        "input": input_data,
        "detection": detection,
        "recommended_conditions": recommendations,
    }
    
    return output


def convert_raw_recommendation_to_standard(
    raw_rec: Dict[str, Any],
    rank: int,
    reaction_smiles: str,
    add_reasoning: bool = True,
) -> Dict[str, Any]:
    """
    Convert a raw recommendation dictionary (from recommend_from_reaction)
    to the standard format.
    
    Args:
        raw_rec: Raw recommendation from chemtools.recommend
        rank: Rank number (1, 2, 3, ...)
        reaction_smiles: Original reaction SMILES
        add_reasoning: Whether to add auto-generated reasoning
    
    Returns:
        Standardized recommendation dictionary
    """
    
    # Extract reagents from combo and chemicals
    combo = raw_rec.get('combo', {})
    chemicals = raw_rec.get('chemicals', [])
    
    # Parse reactants from SMILES
    reactants = reaction_smiles.split(">>")[0].split(".")
    
    # Build reagents list in standard format
    reagents = []
    
    # Add starting materials
    for i, reactant_smiles in enumerate(reactants[:2], 1):
        role = "electrophile" if i == 1 else "nucleophile"
        reagents.append({
            "name": None,
            "cas": None,
            "smiles": reactant_smiles,
            "equivalents": 1.0 if i == 1 else 1.2,
            "role": role,
        })
    
    # Add catalyst system from chemicals
    for chem in chemicals:
        role = chem.get('role', '')
        cas = chem.get('cas')
        name = chem.get('name')
        
        if role in ['metal_precursor', 'ligand', 'base', 'solvent', 'additive']:
            reagent = {
                "name": name,
                "cas": cas,
                "smiles": None,
                "equivalents": None,
                "role": role,
            }
            
            # Try to get equivalents from combo
            if role == 'base':
                # Bases typically use 2 eq
                reagent['equivalents'] = 2.0
            elif role in ['metal_precursor', 'ligand']:
                # Catalysts typically use 0.05-0.1 eq (5-10 mol%)
                reagent['equivalents'] = 0.1
            
            reagents.append(reagent)
    
    # Extract conditions
    conditions_data = raw_rec.get('conditions', {})
    temp_str = conditions_data.get('temperature', '')
    time_str = conditions_data.get('time', '')
    
    # Parse temperature (e.g., "100 °C" -> 100)
    import re
    temp_match = re.search(r'(\d+)', temp_str)
    temperature = float(temp_match.group(1)) if temp_match else None
    
    # Parse time (e.g., "12 h" -> 12)
    time_match = re.search(r'(\d+)', time_str)
    time_hours = float(time_match.group(1)) if time_match else None
    
    conditions = {}
    
    if temperature:
        conditions["temperature"] = f"{int(temperature)}°C"
        if temperature >= 100:
            conditions["temperature_range"] = f"{int(temperature-20)}-{int(temperature+20)}°C"
    
    if time_hours:
        conditions["time"] = f"{int(time_hours)} hours"
        conditions["time_range"] = f"{max(1, int(time_hours/2))}-{int(time_hours*2)} hours"
    
    # Add atmosphere for metal-catalyzed reactions
    conditions["atmosphere"] = "Inert (N₂ or Ar)"
    
    # Build recommendation
    recommendation = {
        "rank": rank,
        "reaction": {
            "smiles": reaction_smiles,
        },
        "reagents": reagents,
        "conditions": conditions,
    }
    
    # Add reasoning if requested
    if add_reasoning:
        precedent_count = raw_rec.get('precedent_support', 0)
        similarity = raw_rec.get('similarity_score', 0.0)
        
        reasoning = {
            "method": "precedent-based",
            "basis": f"Based on analysis of {precedent_count} similar precedent reactions from literature database",
            "confidence_factors": [],
        }
        
        if precedent_count > 100:
            reasoning["confidence_factors"].append(f"Large precedent dataset ({precedent_count} reactions) provides high confidence")
        elif precedent_count > 10:
            reasoning["confidence_factors"].append(f"Good precedent support ({precedent_count} reactions)")
        
        if similarity > 0.8:
            reasoning["confidence_factors"].append(f"High structural similarity to known reactions (score: {similarity:.2f})")
        
        reasoning["confidence_factors"].append("Conditions optimized for this reaction class")
        
        recommendation["reasoning"] = reasoning
    
    return recommendation
