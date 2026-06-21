"""Motif detection and classification system.

This subpackage provides taxonomy-driven motif (compound) detection using
compiled SMARTS patterns from organic groups and compounds definitions.

Public API
----------
**Data Models:**
- CompoundPattern: Compiled SMARTS pattern with metadata

**Registry:**
- build_compound_registry: Load and compile motif registry from taxonomy files

**Detection:**
- detect_motifs: Detect compound motifs in molecules
- discover_undocumented_motifs: Find undocumented scaffold-substituent combinations

**Classification:**
- classify_compound_smiles: Classify SMILES into compound motif IDs
- classify_compound_batch: Batch classification
- choose_best_compound_hit: Select most specific motif from matches

**Utilities:**
- calculate_smarts_complexity: Compute SMARTS specificity score
"""

from .classification import (
    choose_best_compound_hit,
    classify_compound_batch,
    classify_compound_smiles,
)
from .detection import detect_motifs, discover_undocumented_motifs
from .models import CompoundPattern
from .registry import build_compound_registry, _default_registry_paths
from .utils import calculate_smarts_complexity

__all__ = [
    # Models
    "CompoundPattern",
    # Registry
    "build_compound_registry",
    "_default_registry_paths",
    # Detection
    "detect_motifs",
    "discover_undocumented_motifs",
    # Classification
    "classify_compound_smiles",
    "classify_compound_batch",
    "choose_best_compound_hit",
    # Utilities
    "calculate_smarts_complexity",
]
