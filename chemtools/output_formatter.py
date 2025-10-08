"""
Enhanced output formatter for reaction condition recommendations.

Provides comprehensive structured output following the improved format specification.
"""

from __future__ import annotations

import datetime
from typing import Any, Dict, List, Optional, Tuple
from . import reagent_lookup


def format_meta(
    model_type: str = "ML-precedent-knn",
    status: str = "success",
    processing_time_ms: Optional[float] = None,
    request_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Format metadata section."""
    meta = {
        "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        "model": model_type,
        "model_version": "1.0.0",
        "status": status,
    }
    
    if request_id:
        meta["request_id"] = request_id
    
    if processing_time_ms is not None:
        meta["processing_time_ms"] = round(processing_time_ms, 2)
    
    return meta


def format_input(
    reaction_smiles: str,
    requested_reaction_type: Optional[str] = None,
    constraints: Optional[Dict[str, Any]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Format input section."""
    input_data = {
        "reaction_smiles": reaction_smiles,
        "requested_reaction_type": requested_reaction_type,
    }
    
    if constraints:
        input_data["constraints"] = constraints
    
    if options:
        input_data["options"] = options
    
    return input_data


def format_detection(
    detected_type: str,
    confidence: float,
    method: str = "auto",
    alternatives: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Format detection section."""
    detection = {
        "detected_reaction_type": detected_type,
        "confidence": round(confidence, 4),
        "method": method,
    }
    
    if alternatives:
        detection["alternative_types"] = alternatives
    
    return detection


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
        db_info = reagent_lookup.enrich_reagent_info(name, reagent_type)
        
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
    output = {
        "meta": format_meta(
            model_type="ML-precedent-knn",
            status="success",
            processing_time_ms=processing_time_ms,
        ),
        "input": format_input(
            reaction_smiles=reaction_smiles,
            requested_reaction_type=requested_type,
        ),
        "detection": format_detection(
            detected_type=detected_type,
            confidence=detection_confidence,
            method="rxn-insight-ml",
        ),
        "recommendations": [],
    }
    
    # Format each recommendation
    for i, rec_data in enumerate(recommendations_data[:3], 1):  # Limit to 3
        output["recommendations"].append(rec_data)
    
    return output


def format_rule_output(
    reaction_smiles: str,
    requested_type: Optional[str],
    detected_type: str,
    recommendations_data: List[Dict[str, Any]],
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
    output = {
        "meta": format_meta(
            model_type=f"Rule-based-{database_name}",
            status="success",
            processing_time_ms=processing_time_ms,
        ),
        "input": format_input(
            reaction_smiles=reaction_smiles,
            requested_reaction_type=requested_type,
        ),
        "detection": format_detection(
            detected_type=detected_type,
            confidence=1.0,  # Rule-based is deterministic
            method="pattern-match",
        ),
        "recommendations": [],
    }
    
    # Format each recommendation
    for i, rec_data in enumerate(recommendations_data[:3], 1):  # Limit to 3
        output["recommendations"].append(rec_data)
    
    return output


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
