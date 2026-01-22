"""
ML-based output formatting and standard output builders.

This module handles formatting of ML recommendation outputs and provides
the core builder function for creating standardized recommendation outputs.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List, Optional
from ..recommend.utils import friendly_family_label
from .base import format_meta, format_input, format_detection
from .normalization import normalize_recommendations


def build_standard_output(
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
    """
    Construct the canonical output structure shared across recommendation modes.
    
    This is the core builder function that assembles all sections into the
    standard output format.
    
    Args:
        model_type: Model identifier (ML-precedent-knn, Fusion-hybrid, etc.)
        reaction_smiles: Input reaction SMILES
        requested_type: User-requested reaction type
        detected_type: Auto-detected reaction type
        detection_confidence: Detection confidence score
        detection_method: Detection method used
        recommendations_data: List of recommendation dictionaries
        processing_time_ms: Processing time in milliseconds
        status: Status (success, error)
        extras: Additional metadata
        
    Returns:
        Complete standardized output dictionary
    """
    normalized_recommendations = normalize_recommendations(recommendations_data)
    
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
            detected_family_label=friendly_family_label(detected_type),
        ),
        "detection": format_detection(
            detected_type=detected_type or requested_type,
            confidence=detection_confidence,
            method=detection_method,
            family_label=friendly_family_label(detected_type or requested_type),
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
    
    Args:
        response: Response dictionary to validate/convert
        default_model: Default model type if not present
        fallback_reaction_smiles: Fallback reaction SMILES
        fallback_requested_type: Fallback requested type
        extras: Additional metadata to include
        
    Returns:
        Standardized output dictionary
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
    
    return build_standard_output(
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
    return build_standard_output(
        model_type="ML-precedent-knn",
        reaction_smiles=reaction_smiles,
        requested_type=requested_type,
        detected_type=detected_type,
        detection_confidence=detection_confidence,
        detection_method="rxn-insight-ml",
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
    )


def format_fusion_output(
    reaction_smiles: str,
    requested_type: Optional[str],
    fusion_result: Dict[str, Any],
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Format fusion recommendation results into the canonical schema.
    
    Args:
        reaction_smiles: Input reaction SMILES
        requested_type: User-requested reaction type
        fusion_result: Fusion result dictionary
        processing_time_ms: Processing time in milliseconds
        
    Returns:
        Formatted fusion output dictionary
    """
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
    
    return build_standard_output(
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


