"""
Recommendation Service - ML-based and rule-based condition recommendations.

This service handles:
- ML-based condition recommendations
- Rule-based recommendations with reranking
- Plate design from reactions
- Recommendation formatting and post-processing
"""

from typing import Dict, Any, Optional
import time
from chemtools import chem, output_formatter
from chemtools.contracts import (
    RecommendConditionsRequest,
    RecommendFromReactionRequest,
    PlateDesignRequest,
)
from chemtools.exceptions import ValidationError


def recommend_conditions(req: RecommendConditionsRequest) -> Dict[str, Any]:
    """
    Get ML-based condition recommendations.
    
    Returns clean, compact output matching the ui_simple.py format.
    This is robot-friendly and excludes large feature vectors.
    
    Args:
        req: RecommendConditionsRequest with reaction and parameters
        
    Returns:
        Formatted recommendation output with:
        - detection: Reaction type detection info
        - recommendations: List of recommended conditions
        - meta: Processing metadata
        
    Raises:
        ValidationError: If reaction SMILES is invalid
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    start_time = time.perf_counter()
    
    # Get raw ML recommendations
    raw_data = chem.recommend.conditions(
        reaction=req.reaction,
        reaction_type=req.reaction_type,
        k=req.k,
        limit=req.limit,
        relax=req.relax or {},
        constraints=req.constraints or {},
    )
    
    # Extract detection info
    detection = raw_data.get("detection", {})
    detected_type = detection.get("reaction_type", "Unknown")
    confidence = detection.get("auto", {}).get("confidence", 0.0) if detection.get("source") == "auto" else 1.0
    
    # Format recommendations using the UI formatter
    recommendations_data = []
    for rec in raw_data.get("recommendations", [])[:req.limit]:
        summary = rec.get("summary", {})
        chemicals = rec.get("chemicals", [])
        conditions_info = rec.get("conditions", {})
        
        # Build reagents list
        reagents = []
        for chemical in chemicals:
            reagents.append({
                "id": chemical.get("uid", chemical.get("cas")),
                "role": chemical.get("role", "reagent"),
                "name": chemical.get("name"),
                "abbreviation": None,
                "cas": chemical.get("cas"),
                "smiles": chemical.get("smiles"),
                "equivalents": None,  # Not in summary
            })
        
        # Format conditions
        conditions = {}
        if conditions_info.get("temperature"):
            conditions["temperature"] = {
                "value": conditions_info["temperature"],
                "unit": "°C"
            }
        if conditions_info.get("time"):
            conditions["time"] = {
                "value": conditions_info["time"],
                "unit": "hours"
            }
        if conditions_info.get("atmosphere"):
            conditions["atmosphere"] = conditions_info["atmosphere"]
        
        formatted_rec = {
            "rank": rec.get("rank", len(recommendations_data) + 1),
            "confidence": summary.get("confidence", 0.0) / 100.0 if summary.get("confidence") else 0.0,  # Convert % to decimal
            "reagents": reagents,
            "conditions": conditions,
            "precedent_count": summary.get("support", {}).get("count") if isinstance(summary.get("support"), dict) else summary.get("support", 0),
        }
        
        recommendations_data.append(formatted_rec)
    
    # Build compact output using UI formatter
    processing_time_ms = (time.perf_counter() - start_time) * 1000
    
    output = output_formatter.format_ml_output(
        reaction_smiles=req.reaction,
        requested_type=req.reaction_type,
        detected_type=detected_type,
        detection_confidence=confidence,
        recommendations_data=recommendations_data,
        processing_time_ms=processing_time_ms,
    )
    
    return output


def recommend_from_reaction(req: RecommendFromReactionRequest) -> Dict[str, Any]:
    """
    Get recommendations from a reaction using ML or rule-based strategies.
    
    Supports multiple reranking strategies:
    - 'ml': ML-based reranking (default)
    - 'rule': Rule-based reranking
    - 'drfp': DRFP similarity-based reranking
    
    Args:
        req: RecommendFromReactionRequest with reaction and parameters
        
    Returns:
        Standardized recommendation output with:
        - input: Input reaction info
        - detection: Reaction type detection
        - recommendations: List of recommended conditions
        - meta: Processing metadata
        
    Raises:
        ValidationError: If reaction SMILES is invalid
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    start_time = time.perf_counter()
    
    # Get raw recommendations
    raw_result = chem.recommend.recommend_from_reaction(
        req.reaction,
        k=req.k,
        relax=req.relax or {},
        constraint_rules=req.constraints or {},
        rerank_strategy=req.rerank_strategy,
        filter_unknown_reagents=req.filter_unknown_reagents,
    )
    
    processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)
    
    # Get formatted output
    formatted_source = raw_result.get("formatted") or raw_result
    
    # Ensure standard output format
    standard = output_formatter.ensure_standard_output(
        formatted_source,
        default_model="ML-precedent-knn",
        fallback_reaction_smiles=req.reaction,
        extras={"raw_recommendation": raw_result},
    )
    
    # Add processing time to metadata
    standard_meta = standard.setdefault("meta", {})
    standard_meta["processing_time_ms"] = processing_time_ms
    
    return standard


def design_plate(req: PlateDesignRequest) -> Dict[str, Any]:
    """
    Design an experiment plate from a reaction.
    
    Generates diverse conditions for experimental screening.
    
    Args:
        req: PlateDesignRequest with reaction and plate parameters
        
    Returns:
        Plate design with recommended conditions
        
    Raises:
        ValidationError: If reaction SMILES is invalid
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")
    
    return chem.recommend.design_plate_from_reaction(
        req.reaction,
        plate_size=req.plate_size,
        relax=req.relax or {},
        constraint_rules=req.constraints or {},
    )
