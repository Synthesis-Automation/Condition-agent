"""
Utility functions for feature formatting.

Common helpers for motif ID extraction, normalization, and score handling.
"""

from __future__ import annotations

from typing import Any, Iterable, List
import re

_HETERO_C_H_SCAFFOLD_GROUPS = {"HeteroAr", "AromN"}


def extract_motif_ids(*motif_collections: Iterable[Any]) -> List[str]:
    """
    Extract compound IDs from one or more motif collections.
    
    Args:
        *motif_collections: Variable number of motif iterables
        
    Returns:
        List of normalized compound ID strings
    """
    items: List[str] = []
    for collection in motif_collections:
        for motif in collection or []:
            if isinstance(motif, dict):
                cid = motif.get("compound_id") or motif.get("id")
                if cid:
                    items.append(normalize_motif_id(str(cid)))
            else:
                text = str(motif).strip()
                if text:
                    items.append(normalize_motif_id(text))
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
    Normalize motif ID by removing double dashes, extra whitespace, and aliases.
    
    Args:
        motif_id: Raw motif ID string
        
    Returns:
        Cleaned motif ID
    """
    text = str(motif_id).strip()
    # Replace multiple consecutive dashes with single dash
    text = re.sub(r"-{2,}", "-", text)
    # Canonicalize legacy pinacol boronate motif labels to generic boronate ester.
    if text.endswith("-Bpin"):
        text = text[:-5] + "-B(OR)2"
    elif text in {"Bpin", "-Bpin"}:
        text = "B(OR)2"
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
        scaffold, group_id = text.rsplit("-", 1)
        scaffold = scaffold.strip()
        group_id = group_id.strip()
        # Preserve heteroaromatic scaffold identity for C-H spectator motifs.
        if group_id == "H" and scaffold in _HETERO_C_H_SCAFFOLD_GROUPS:
            return scaffold
        return group_id
    return text
