"""
Backward compatibility shim for reaction_detection module.

THIS MODULE IS DEPRECATED. Import from chemtools.featurizers.detection instead.

Refactored into subpackage structure:
- detection/models.py: Data structures
- detection/bond_changes.py: Bond change analysis
- detection/matchers.py: Motif-based matching
- detection/coupling.py: Coupling confirmation
- detection/core.py: Main API
"""

# Re-export everything from the new detection subpackage for backward compatibility
from .detection import (
    detect_reaction_types,
    detect_reaction_types_from_smiles,
    detect_motif_ids_from_smiles,
    clear_bond_change_cache,
    ReactionMatch,
    ReactionDetectionResult,
)

__all__ = [
    "detect_reaction_types",
    "detect_reaction_types_from_smiles",
    "detect_motif_ids_from_smiles",
    "clear_bond_change_cache",
    "ReactionMatch",
    "ReactionDetectionResult",
]
