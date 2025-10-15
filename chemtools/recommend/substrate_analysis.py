"""
Generic substrate analysis utilities.

Functions for analyzing molecular properties and functional groups.
NOTE: This module now uses chemtools.util.functional_groups for detection.
"""

from __future__ import annotations

from typing import Optional

# Import from centralized functional_groups utility
from ..util.functional_groups import (
    has_free_alcohol,
    has_phenol,
    has_sulfonamide,
    has_hydroxylamine,
)


def detect_functional_groups(smiles: Optional[str]) -> dict[str, bool]:
    """
    Detect common functional groups in molecule.
    
    NOTE: This is a compatibility wrapper. For comprehensive functional group
    detection, use chemtools.util.functional_groups directly or via the
    Context API: chem.functional_groups.detect(smiles)
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dict of functional group flags (limited set for backward compatibility)
    """
    return {
        "free_alcohol": has_free_alcohol(smiles),
        "phenol": has_phenol(smiles),
        "sulfonamide": has_sulfonamide(smiles),
        "hydroxylamine": has_hydroxylamine(smiles),
    }
