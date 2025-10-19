"""
Base formatting functions for metadata, input, and detection sections.

This module provides core formatting functions that build the standard
output structure sections used across all recommendation types.
"""

from __future__ import annotations

import datetime
from typing import Any, Dict, List, Optional


def format_meta(
    model_type: str = "ML-precedent-knn",
    status: str = "success",
    processing_time_ms: Optional[float] = None,
    request_id: Optional[str] = None,
    schema_version: str = "2.0",
    model_version: str = "1.0.0",
) -> Dict[str, Any]:
    """
    Format metadata section.
    
    Args:
        model_type: Type of model used (ML-precedent-knn, Rule-based, etc.)
        status: Execution status (success, error)
        processing_time_ms: Processing time in milliseconds
        request_id: Optional unique request identifier
        schema_version: Output schema version
        model_version: Model version string
        
    Returns:
        Formatted metadata dictionary
    """
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
    """
    Format input section.
    
    Args:
        reaction_smiles: Input reaction SMILES string
        requested_reaction_type: User-requested reaction type/family
        detected_family: Automatically detected family
        constraints: Optional constraint rules
        options: Optional processing options
        
    Returns:
        Formatted input dictionary
    """
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
    """
    Format detection section.
    
    Args:
        detected_type: Detected reaction type/family
        confidence: Detection confidence score (0.0 to 1.0)
        method: Detection method (auto, manual, pattern-match, etc.)
        alternatives: Optional list of alternative type detections
        
    Returns:
        Formatted detection dictionary
    """
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
