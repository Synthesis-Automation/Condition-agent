"""Precedent search module for ChemTools.

Split from monolithic precedent.py into focused submodules:
- loader: Dataset loading and transformation
- core_utils: Family name normalization and parsing utilities
- catalyst: Catalyst class detection and matching
- similarity: Feature similarity calculations
- search: Main k-NN search and core-based lookup

Public API (backward compatible):
- knn: Find k nearest neighbor precedents
- find_reactions_by_core: Find reactions by condition core
- list_cores: List unique condition cores

Internal API (for chemtools.context):
- _load_selective: Load datasets filtered by family
"""
from .search import knn, find_reactions_by_core, list_cores
from .loader import _load_selective

__all__ = [
    "knn",
    "find_reactions_by_core",
    "list_cores",
]
