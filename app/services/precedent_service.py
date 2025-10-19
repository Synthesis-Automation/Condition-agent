"""
Precedent Service - Precedent search and core management.

This service handles:
- Precedent KNN search
- Condition core listing
- Core-based reaction search
- Constraint filtering
- Precedent explanation
- Condition core parsing and validation
"""

from typing import Dict, Any, List, Optional
import json
import os
from chemtools import chem, condition_core
from chemtools.contracts import (
    PrecedentKNNRequest,
    ConstraintsFilterRequest,
    ExplainPrecedentsRequest,
    ConditionCoreParseRequest,
    ConditionCoreValidateRequest,
    CoreSearchRequest,
)
from chemtools.exceptions import ValidationError, DatabaseNotFoundError


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


def parse_condition_core(req: ConditionCoreParseRequest) -> Dict[str, Any]:
    """
    Parse reagents into condition core representation.
    
    Args:
        req: ConditionCoreParseRequest with reagents and optional text
        
    Returns:
        Parsed condition core
        
    Raises:
        ValidationError: If reagents are invalid
    """
    if not req.reagents:
        raise ValidationError("Reagents list cannot be empty")
    
    return condition_core.parse_core(req.reagents, req.text or "")


def validate_dataset(req: ConditionCoreValidateRequest) -> Dict[str, Any]:
    """
    Validate condition core parsing against a dataset.
    
    Loads a JSONL dataset and validates that condition_core parsing
    matches the ground truth cores in the dataset.
    
    Args:
        req: ConditionCoreValidateRequest with path and validation options
        
    Returns:
        Dictionary with:
        - records: Total records processed
        - matches: Number of matching cores
        - accuracy: Match percentage
        - mismatches: List of mismatches (up to show_mismatches limit)
        
    Raises:
        DatabaseNotFoundError: If dataset file not found
        ValidationError: If dataset is invalid
    """
    path = req.path
    limit = int(req.limit or 0)
    
    if not os.path.exists(path):
        raise DatabaseNotFoundError(f"Dataset not found: {path}")
    
    total = 0
    ok = 0
    mismatches = []
    
    def _norm_core(s: str) -> str:
        """Normalize core string for comparison."""
        s = (s or "").strip()
        return s[:-5] if s.endswith("/none") else s
    
    def _metal_part(s: str) -> str:
        """Extract metal part from core."""
        s = (s or "").strip()
        return s.split("/", 1)[0] if s else ""
    
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                
                # Extract reagents from record
                cat = rec.get("catalyst") or {}
                reagents = []
                
                # Add catalyst components
                for item in (cat.get("full_system") or cat.get("core") or []):
                    reagents.append({
                        "name": item.get("name"),
                        "uid": item.get("cas"),
                        "role": "CATALYST",
                    })
                
                # Add reagents
                for item in (rec.get("reagents") or []):
                    reagents.append({
                        "name": item.get("name"),
                        "uid": item.get("cas"),
                        "role": item.get("role") or "ADDITIVE",
                    })
                
                # Add solvents
                for item in (rec.get("solvents") or []):
                    reagents.append({
                        "name": item.get("name"),
                        "uid": item.get("cas"),
                        "role": "SOLVENT",
                    })
                
                # Parse condition core
                out = condition_core.parse(reagents, "")
                truth = (rec.get("condition_core") or "").strip()
                pred = (out.get("core") or "").strip()
                
                # Normalize for comparison
                t = _norm_core(truth)
                p = _norm_core(pred)
                
                # Check match (with metal_only_ok option)
                ok_flag = (t == p) or (
                    req.metal_only_ok and _metal_part(t) and _metal_part(t) == _metal_part(p)
                )
                
                total += 1
                if ok_flag:
                    ok += 1
                elif len(mismatches) < int(req.show_mismatches or 0):
                    mismatches.append({
                        "reaction_id": rec.get("reaction_id"),
                        "truth": truth,
                        "pred": pred,
                    })
                
                # Stop if limit reached
                if limit and total >= limit:
                    break
    
    except FileNotFoundError as exc:
        raise DatabaseNotFoundError(f"Dataset not found: {path}") from exc
    except Exception as exc:
        raise ValidationError(f"Error processing dataset: {exc}") from exc
    
    # Calculate accuracy
    acc = (ok / total) * 100.0 if total else 0.0
    
    return {
        "records": total,
        "matches": ok,
        "accuracy": round(acc, 2),
        "mismatches": mismatches,
    }
