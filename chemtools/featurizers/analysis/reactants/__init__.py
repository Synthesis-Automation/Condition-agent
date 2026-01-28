"""Reactant classification and reaction type utilities.

This subpackage provides taxonomy-driven reactant classification and
reaction type lookups, with all definitions loaded from the canonical
taxonomy files (organic_compounds.v1.3.json and reaction_types.v4.0.json).

Public API
----------
- ReactantMatch: Result structure for reactant classification
- classify_reactant_smiles: Get best reactant type match for SMILES
- classify_reactant_category: Get category ID for SMILES
- iter_reactant_matches: Get all reactant type matches
- normalize_reactant_identifier: Resolve aliases to canonical IDs
- get_reactant_type_definitions: Load reactant taxonomy
- normalize_reaction_type: Resolve reaction type aliases
- get_reaction_type_definitions: Load reaction type taxonomy
- required_reactant_categories: Get required categories for reaction
"""

from .classification import (
    ReactantMatch,
    build_reactant_lookup,
    classify_reactant_batch,
    classify_reactant_category,
    classify_reactant_group,
    classify_reactant_smiles,
    get_all_reactant_matches,
    get_reactant_category_matches,
    iter_reactant_matches,
    reactant_category_for,
)
from .reactions import (
    build_reaction_lookup,
    clear_reaction_type_cache,
    describe_reaction_type,
    get_reaction_type_definitions,
    get_reaction_types_file,
    iter_reactions_for_category,
    list_reaction_type_ids,
    normalize_reaction_type,
    required_reactant_categories,
)
from .taxonomy import (
    clear_reactant_type_cache,
    get_reactant_type_definitions,
    get_reactant_types_file,
    normalize_reactant_identifier,
)

__all__ = [
    # Classification
    "ReactantMatch",
    "classify_reactant_smiles",
    "classify_reactant_category",
    "classify_reactant_group",
    "classify_reactant_batch",
    "iter_reactant_matches",
    "get_reactant_category_matches",
    "get_all_reactant_matches",
    "reactant_category_for",
    "build_reactant_lookup",
    # Taxonomy
    "get_reactant_type_definitions",
    "get_reactant_types_file",
    "clear_reactant_type_cache",
    "normalize_reactant_identifier",
    # Reactions
    "get_reaction_type_definitions",
    "get_reaction_types_file",
    "clear_reaction_type_cache",
    "normalize_reaction_type",
    "describe_reaction_type",
    "list_reaction_type_ids",
    "build_reaction_lookup",
    "iter_reactions_for_category",
    "required_reactant_categories",
]
