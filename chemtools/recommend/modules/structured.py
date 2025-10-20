"""
Structured output wrapper for API/UI consumers.

This module provides a wrapper around the core recommendation function
with enhanced metadata, timing, and structured formatting.
"""

from __future__ import annotations

import time
from datetime import datetime, timezone
from typing import Dict, Any, Optional


def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
    rerank_strategy: str = 'rule',
    filter_unknown_reagents: bool = False,
    search_all_families: bool = False,
    _recommend_fn=None,  # Injected to avoid circular imports
) -> Dict[str, Any]:
    """
    Recommend conditions with structured output format for API/UI consumers.
    
    This is a wrapper around recommend_from_reaction() that provides
    a more structured output format with proper metadata and precedent summaries.
    
    Args:
        reaction: Reaction SMILES string
        reaction_type: Optional reaction type override (ignored if search_all_families=True)
        k: Number of precedents (default: 50 for better coverage)
        limit: Maximum number of recommendations to return (default: 5)
        relax: Relaxation parameters for precedent search
        constraints: Constraint rules (inventory, blacklist, etc.)
        rerank_strategy: 'rule', 'analytics', or 'none' (default: 'rule')
        filter_unknown_reagents: Filter precedents with unknown reagents (default: False)
        search_all_families: Search across all reaction families (default: False)
    
    Returns:
        Dict with keys:
            - meta: Metadata (status, timestamp, strategy, result_count)
            - input: Normalized reaction SMILES and family info
            - detection: Family detection info with auto/rule sources
            - recommendations: List of condition variants (limited to `limit`)
            - alternatives: Alternative cores, bases, solvents
            - precedents: Precedent summary with top examples
    """
    # Start timing
    start_time = time.time()
    
    limit = max(1, int(limit or 1))
    cfg_relax = dict(relax or {})
    
    # Import here to avoid circular dependency
    if _recommend_fn is None:
        from .recommender import recommend_from_reaction as _recommend_fn
    
    result = _recommend_fn(
        reaction=reaction,
        k=int(k or 50),
        relax=cfg_relax,
        constraint_rules=constraints or {},
        family_override=reaction_type,
        max_variants=limit,
        rerank_strategy=rerank_strategy,
        filter_unknown_reagents=filter_unknown_reagents,
        search_all_families=search_all_families,
    )
    
    # Calculate processing time
    processing_time_ms = round((time.time() - start_time) * 1000, 2)
    
    # Extract formatted section
    formatted = dict(result.get("formatted") or {})
    recommendations = list(formatted.get("recommended_conditions") or [])
    recommendations = recommendations[:limit]
    
    # Ensure rank is set correctly
    for idx, rec in enumerate(recommendations, start=1):
        rec.setdefault("rank", idx)
        summary = rec.get("summary")
        if isinstance(summary, dict):
            summary.setdefault("rank", idx)
    
    # Build detection section with all fields
    detection = dict(formatted.get("detection") or {})
    detected_family = detection.get("family") or result.get("family") or "Unknown"
    detection_confidence = detection.get("confidence", 0.95)
    
    if reaction_type and not detection.get("source"):
        detection["source"] = "user_supplied"
    detection.setdefault("source", detection.get("source") or "auto")
    detection["detected_reaction_type"] = detected_family
    detection["confidence"] = detection_confidence
    detection.setdefault("method", "drfp_precedent_search")
    detection.setdefault("provided_reaction_type", reaction_type)
    
    # Build meta section with all fields
    meta = dict(formatted.get("meta") or {})
    meta["generated_at"] = datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')
    meta.setdefault("model", "drfp_similarity")
    meta.setdefault("status", "success")
    meta.setdefault("strategy", "drfp_similarity")
    meta["result_count"] = len(recommendations)
    meta["processing_time_ms"] = processing_time_ms
    
    # Build input section with all fields
    input_data = dict(formatted.get("input") or {})
    input_data["reaction_smiles"] = input_data.get("reaction_smiles") or reaction
    input_data["requested_reaction_type"] = reaction_type
    input_data.setdefault("family", detected_family)
    
    # Build precedent summary
    pack = result.get("precedent_pack") or {}
    precedents = list(pack.get("precedents") or [])
    top_precedents = [
        {
            "reaction_id": p.get("reaction_id"),
            "core": p.get("core"),
            "yield_pct": p.get("yield_pct"),
        }
        for p in precedents[:10]
        if p.get("reaction_id")
    ]
    
    core_family = result.get("family")
    core_support = len([p for p in precedents if p.get("core") == core_family])
    
    precedent_summary = {
        "total_considered": len(precedents),
        "core_family": core_family,
        "core_support": core_support,
        "top_precedents": top_precedents,
    }
    
    # Include detailed precedent information from formatted output
    precedents_used = formatted.get("precedents_used")
    
    return {
        "meta": meta,
        "input": input_data,
        "detection": detection,
        "recommendations": recommendations,
        "alternatives": result.get("alternatives"),
        "precedents": precedent_summary,
        "precedents_used": precedents_used,
        "constraint_filters": result.get("constraint_filters"),
    }
