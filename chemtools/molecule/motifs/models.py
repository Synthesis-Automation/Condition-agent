"""Data models for motif detection system."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, List

# Regex for stripping atom map numbers from SMARTS
_MAP_RE = re.compile(r":\d+(?=\])")

# Default template patterns for compound SMARTS generation
_DEFAULT_TEMPLATES = {
    "single_bond": "{A}-{B}",
    "via_oxygen": "{A}O{B}",
}

# Skip these substituent types in discovery mode
_DISCOVERY_SKIP_SUBSTITUENTS = {
    "Alkyl_Subst",
    "Alkenyl_Subst",
}


@dataclass(frozen=True)
class CompoundPattern:
    """Compiled SMARTS pattern for motif detection.
    
    Attributes:
        compound_id: Unique identifier for the compound
        smarts: SMARTS pattern string
        query: Compiled RDKit query molecule
        group_a: Scaffold group identifier
        group_b: Substituent group identifier
        b_tags: Additional tags for the substituent
        priority: Detection priority (higher = more specific)
        complexity: Structural complexity score
        reactivity_weight: Reactivity importance weight
    """
    compound_id: str
    smarts: str
    query: Any
    group_a: str
    group_b: str
    b_tags: List[str]
    priority: int = 1
    complexity: int = 0
    reactivity_weight: float = 0.0


__all__ = [
    "CompoundPattern",
    "_MAP_RE",
    "_DEFAULT_TEMPLATES",
    "_DISCOVERY_SKIP_SUBSTITUENTS",
]
