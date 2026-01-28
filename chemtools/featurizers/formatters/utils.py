"""
Utility functions for feature formatting.

Common helpers for motif ID extraction, normalization, and score handling.
"""

from __future__ import annotations

from typing import Any, Iterable, List
import re


def extract_motif_ids(*motif_collections: Iterable[Any]) -> List[str]:
    """
    Extract compound IDs from one or more motif collections.
    
    Args:
        *motif_collections: Variable number of motif iterables
        
    Returns:
        List of compound ID strings
    """
    items: List[str] = []
    for collection in motif_collections:
        for motif in collection or []:
            if isinstance(motif, dict):
                cid = motif.get("compound_id")
                if cid:
                    items.append(str(cid))
            else:
                text = str(motif).strip()
                if text:
                    items.append(text)
    return items


def extract_scores(result: Any) -> List[float]:
    """Extract numerical scores from result dict or list."""
    scores: List[float] = []
    if isinstance(result, dict):
        for key in ("score", "value", "electronic_parameter"):
            value = result.get(key)
            if value is not None:
                try:
                    scores.append(float(value))
                except (TypeError, ValueError):
                    pass
    elif isinstance(result, (list, tuple)):
        for item in result:
            if isinstance(item, dict):
                scores.extend(extract_scores(item))
            else:
                try:
                    scores.append(float(item))
                except (TypeError, ValueError):
                    pass
    return scores


def normalize_motif_id(motif_id: str) -> str:
    """
    Normalize motif ID by removing double dashes and extra whitespace.
    
    Args:
        motif_id: Raw motif ID string
        
    Returns:
        Cleaned motif ID
    """
    text = str(motif_id).strip()
    # Replace multiple consecutive dashes with single dash
    text = re.sub(r"-{2,}", "-", text)
    # Remove leading/trailing dashes
    text = text.strip("-")
    return text


def group_id_from_motif_id(motif_id: str) -> str:
    """
    Extract group ID from a full motif compound ID.
    
    Examples:
        "Ar-Br" -> "Br"
        "RCH2-NH2" -> "NH2"
        "Pyridine" -> "Pyridine"
    
    Args:
        motif_id: Full compound ID
        
    Returns:
        Group ID portion
    """
    text = str(motif_id).strip()
    if not text:
        return ""
    if "-" in text:
        return text.split("-")[-1].strip()
    return text
