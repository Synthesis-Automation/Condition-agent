"""
Recommendation Service - HTE condition recommendations.

This service handles:
- HTE-based recommendations (primary backend)
- Recommendation formatting and post-processing
"""

from typing import Dict, Any, Optional
import time
from chemtools import chem
from chemtools.formatters import ensure_standard_output
from chemtools.core.contracts import (
    RecommendConditionsRequest,
    RecommendFromReactionRequest,
    PlateDesignRequest,
)
from chemtools.core.errors import ValidationError


def recommend_conditions(req: RecommendConditionsRequest) -> Dict[str, Any]:
    """
    Get ML-based condition recommendations.

    Returns an HTE-backed recommendation payload.
    """
    if not req.reaction or not req.reaction.strip():
        raise ValidationError("Reaction SMILES cannot be empty")

    raw_data = chem.recommend.conditions(
        reaction=req.reaction,
        reaction_type=req.reaction_type,
        k=req.k,
        limit=req.limit,
        relax=req.relax or {},
        constraints=req.constraints or {},
    )

    return raw_data


def recommend_from_reaction(req: RecommendFromReactionRequest) -> Dict[str, Any]:
    """
    Get HTE-backed recommendations from a reaction.

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
    
    # Prepare relax dict with cross-family parameters
    relax_dict = req.relax or {}
    if req.search_all_families:
        relax_dict.update({
            'reaction_type_threshold': req.reaction_type_threshold,
            'mechanism_similarity_threshold': req.mechanism_similarity_threshold,
            'mechanism_weight': req.mechanism_weight
        })
    
    # Get raw recommendations
    raw_result = chem.recommend.recommend_from_reaction(
        req.reaction,
        k=req.k,
        relax=relax_dict,
        constraint_rules=req.constraints or {},
        rerank_strategy=req.rerank_strategy,
        filter_unknown_reagents=req.filter_unknown_reagents,
        search_all_families=req.search_all_families,
    )
    
    processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)
    
    # Get formatted output
    formatted_source = raw_result.get("formatted") or raw_result
    
    # Ensure standard output format
    standard = ensure_standard_output(
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
    
    raise ValidationError("Plate design is not available in the unified recommendation system.")
