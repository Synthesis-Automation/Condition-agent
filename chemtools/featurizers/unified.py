"""
Unified feature bundles for molecules and reactions.

SIMPLIFIED: This module now imports from the formatters/ subpackage.
The implementation has been refactored into:
- formatters/molecule.py: Molecule featurization and bundling
- formatters/reaction.py: Reaction featurization and bundling
- formatters/aggregation.py: Feature aggregation logic
- formatters/utils.py: Shared utility functions
"""

from .formatters import (
    featurize_molecule,
    featurize_reaction,
    format_reaction_key,
    select_primary_reacted_motifs,
    select_primary_formed_motif,
)

__all__ = [
    "featurize_molecule",
    "featurize_reaction",
    "format_reaction_key",
    "select_primary_reacted_motifs",
    "select_primary_formed_motif",
]
