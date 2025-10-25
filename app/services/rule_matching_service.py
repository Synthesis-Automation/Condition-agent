"""
Rule Matching Service - SCDB (SchemeConditionDB) rule-based matching.

This service handles:
- Loading and caching SCDB databases
- Rule-based reaction matching
- Catalyst class filtering
- Result formatting with trace information
"""

from typing import Dict, Any, Optional
import os
import time
from pathlib import Path
from chemtools import chem, output_formatter
from chemtools.contracts import SchemeMatchRequest
from chemtools.exceptions import ValidationError, DatabaseNotFoundError, SchemeConditionDBError

try:
    from chemtools.rule import loader as scdb_loader
    _HAS_SCDB = True
except Exception:
    _HAS_SCDB = False

# Default database path from environment
_SCDB_DEFAULT_DB = (
    os.environ.get("SCDB_MATCH_DB_PATH", "cn_coupling_pd_db.json").strip()
    or "cn_coupling_pd_db.json"
)


def match_reaction(req: SchemeMatchRequest) -> Dict[str, Any]:
    """
    Match a reaction to rule-based scheme database.
    
    Uses ChemTools v2.0 chem.rules.* API to:
    1. Load database (with caching)
    2. Match reaction to rules
    3. Apply optional catalyst filtering
    4. Format results
    
    Args:
        req: SchemeMatchRequest with reaction and optional parameters
        
    Returns:
        Formatted match result with:
        - input: Original reaction
        - detection: Match detection info
        - recommended_conditions: Matched conditions
        - meta: Processing metadata
        
    Raises:
        ResourceNotAvailableError: If SCDB is not available
        ValidationError: If reaction is invalid
        DatabaseNotFoundError: If database file not found
        SchemeConditionDBError: If matching fails
    """
    from chemtools.exceptions import ResourceNotAvailableError
    
    if not _HAS_SCDB:
        raise ResourceNotAvailableError("SchemeConditionDB matcher unavailable")
    
    reaction = (req.reaction or "").strip()
    if not reaction:
        raise ValidationError("Reaction must be a non-empty string")
    
    db_path = (req.db or _SCDB_DEFAULT_DB).strip()
    if not db_path:
        raise ValidationError("No database path configured for matching")
    
    start_time = time.perf_counter()
    
    try:
        # Use ChemTools v2.0 chem.rules.* API
        db = chem.rules.load_database(db_path, cache=True)
        result = chem.rules.match(db, reaction)
    except FileNotFoundError as exc:
        raise DatabaseNotFoundError(f"Database file not found: {db_path}") from exc
    except SchemeConditionDBError as exc:
        raise exc
    except Exception as exc:
        raise SchemeConditionDBError(f"Unexpected error while matching: {exc}") from exc
    
    # Convert result to JSON dict
    payload = result.to_json_dict()
    
    # Apply catalyst class filtering if provided via relax parameter
    catalyst_filter = None
    if req.relax and isinstance(req.relax, dict):
        catalyst_filter = req.relax.get("catalyst_class")
    
    if catalyst_filter and payload.get("recommended_conditions"):
        payload = _filter_by_catalyst_class(payload, catalyst_filter)
    
    # Remove trace if not requested
    if not req.include_trace:
        payload.pop("trace", None)
    
    # Calculate processing time
    processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)
    database_label = Path(db_path).name if db_path else "SchemeConditionDB"
    
    # Format output
    return output_formatter.format_rule_match_result(
        reaction_smiles=reaction,
        match_result=payload,
        requested_type=None,
        database_name=database_label,
        processing_time_ms=processing_time_ms,
    )


def _filter_by_catalyst_class(payload: Dict[str, Any], catalyst_filter: str) -> Dict[str, Any]:
    """
    Filter conditions by catalyst class.
    
    Args:
        payload: Match result payload
        catalyst_filter: Catalyst class to filter by (e.g., "Cu", "Pd")
        
    Returns:
        Updated payload with filtered conditions
    """
    recommended_conditions = payload.get("recommended_conditions", [])
    filtered_conditions = []
    
    for cond in recommended_conditions:
        # Extract catalyst class from condition (look in core or full_system)
        condition_core = cond.get("core", "")
        
        # Simple matching: check if catalyst filter (e.g., "Cu") appears in core
        catalyst_match = False
        if isinstance(condition_core, str):
            # Normalize for comparison
            core_lower = condition_core.lower()
            filter_lower = str(catalyst_filter).lower()
            
            # Match if filter appears in core (e.g., "Cu" in "Cu/L1")
            if filter_lower in core_lower or filter_lower in core_lower.split("/")[0]:
                catalyst_match = True
        
        if catalyst_match:
            filtered_conditions.append(cond)
    
    # Update payload with filtered conditions
    original_count = len(recommended_conditions)
    payload["recommended_conditions"] = filtered_conditions
    payload["_catalyst_filtered"] = {
        "requested": catalyst_filter,
        "original_count": original_count,
        "filtered_count": len(filtered_conditions),
    }
    
    return payload


def is_scdb_available() -> bool:
    """Check if SCDB matcher is available."""
    return _HAS_SCDB


def get_default_database_path() -> str:
    """Get the default SCDB database path."""
    return _SCDB_DEFAULT_DB
