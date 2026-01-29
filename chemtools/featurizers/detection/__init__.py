"""
Reaction type detection using motif-based matching and bond change analysis.

Refactored from monolithic reaction_detection.py (820 lines) into focused modules.
"""

from .models import ReactionMatch, ReactionDetectionResult
from .core import (
    detect_reaction_types,
    detect_reaction_types_from_smiles,
    detect_motif_ids_from_smiles,
    clear_bond_change_cache,
)
from .coupling import confirm_coupling_product_by_attachment
from ...taxonomy.reaction_catalog import load_reaction_catalog


# Backward compatibility wrapper for old Suzuki-specific function signature
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
