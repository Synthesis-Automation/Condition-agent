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

# Re-export coupling functions for backward compatibility
from .detection.coupling import (
    confirm_coupling_product_by_attachment,
)

# Re-export taxonomy loader for tests
from ..taxonomy.reaction_catalog import load_reaction_catalog


# Wrapper for old Suzuki-specific function signature
def confirm_suzuki_product_by_attachment(reactants, products):
    """Backward compatibility wrapper for old Suzuki confirmation function.
    
    Args:
        reactants: List of reactant SMILES
        products: List of product SMILES
        
    Returns:
        Tuple of (success: bool, reason: str)
    """
    return confirm_coupling_product_by_attachment(
        reactants, products, "Suzuki_miyaura"
    )


__all__ = [
    "detect_reaction_types",
    "detect_reaction_types_from_smiles",
    "detect_motif_ids_from_smiles",
    "clear_bond_change_cache",
    "ReactionMatch",
    "ReactionDetectionResult",
    "confirm_coupling_product_by_attachment",
    "confirm_suzuki_product_by_attachment",
    "load_reaction_catalog",
]
