"""
Precedent Service - Precedent search and core management.

This service handles:
- Precedent KNN search
- Condition core listing
- Core-based reaction search
- Constraint filtering
- Precedent explanation
"""

from typing import Dict, Any, Optional
from chemtools import chem
from chemtools.core.contracts import (
    PrecedentKNNRequest,
    ConstraintsFilterRequest,
    ExplainPrecedentsRequest,
    CoreSearchRequest,
)
from chemtools.core.errors import ValidationError


def knn_search(req: PrecedentKNNRequest) -> Dict[str, Any]:
    """
    Perform KNN search for precedents.
    
    Args:
        req: PrecedentKNNRequest with family, features, k, and optional relax params
        
    Returns:
        KNN search results with precedents
        
    Raises:
        ValidationError: If parameters are invalid
    """
    if not req.family or not req.family.strip():
        raise ValidationError("Family cannot be empty")
    if not req.features:
        raise ValidationError("Features cannot be empty")
    if req.k < 1:
        raise ValidationError(f"k must be >= 1, got {req.k}")
    
    return chem.precedent.knn(req.family, req.features, req.k, req.relax or {})


def list_cores(
    family: Optional[str] = None,
    limit: int = 200,
    include_counts: bool = True,
) -> Dict[str, Any]:
    """
    List available condition cores.
    
    Args:
        family: Optional family filter
        limit: Maximum number of cores to return
        include_counts: Whether to include usage counts
        
    Returns:
        Dictionary with cores list
    """
    if limit < 1:
        limit = 200
    
    data = chem.precedent.list_cores(
        family=family,
        top_n=int(limit),
        include_counts=bool(include_counts),
    )
    
    return {"cores": data}


def search_by_core(req: CoreSearchRequest) -> Dict[str, Any]:
    """
    Find reactions by condition core.
    
    Args:
        req: CoreSearchRequest with core and optional parameters
        
    Returns:
        Dictionary with:
        - query: Query parameters
        - count: Number of results
        - results: List of matching reactions
        
    Raises:
        ValidationError: If core is invalid
    """
    if not req.core or not req.core.strip():
        raise ValidationError("Core cannot be empty")
    
    limit = int(req.limit or 50)
    if limit < 1:
        limit = 50
    
    rows = chem.precedent.find_reactions_by_core(
        req.core,
        family=req.family,
        fuzzy=bool(req.fuzzy),
        limit=limit,
    )
    
    return {
        "query": {
            "core": req.core,
            "family": req.family,
            "fuzzy": bool(req.fuzzy),
            "limit": limit,
        },
        "count": len(rows),
        "results": rows,
    }


def filter_constraints(req: ConstraintsFilterRequest) -> Dict[str, Any]:
    """
    Filter candidates by constraint rules.
    
    Args:
        req: ConstraintsFilterRequest with candidates and rules
        
    Returns:
        Filtered candidates
        
    Raises:
        ValidationError: If candidates are invalid
    """
    if not req.candidates:
        raise ValidationError("Candidates list cannot be empty")
    
    return chem.constraints.filter(req.candidates, req.rules or {})


def explain_precedents(req: ExplainPrecedentsRequest) -> Dict[str, Any]:
    """
    Explain precedents for a pack.
    
    Args:
        req: ExplainPrecedentsRequest with pack and features
        
    Returns:
        Explanation of precedents
        
    Raises:
        ValidationError: If parameters are invalid
    """
    if not req.pack:
        raise ValidationError("Pack cannot be empty")
    if not req.features:
        raise ValidationError("Features cannot be empty")
    
    return chem.explain.for_pack(req.pack, req.features)


