"""
Reagent registry and taxonomy management system.

This module provides:
- Runtime reagent lookup (lookup.py)
- Taxonomy store and heuristics (taxonomy_store.py)
- CLI tools for taxonomy management (taxonomy_cli.py)
- Utility functions (taxonomy_utils.py)

Public API:
    # Runtime lookup (most commonly used)
    from chemtools.reagent import find_reagent, enrich_reagent_info, enrich_conditions
    
    # Taxonomy management
    from chemtools.reagent import TaxonomyStore, RoleHeuristics
    from chemtools.reagent import normalize_cas, resolve_identity_from_cas
"""

# Runtime lookup functions
from .lookup import (
    find_reagent,
    enrich_reagent_info,
    enrich_conditions,
    load_reagent_database,
    normalize_name,
    get_all_reagent_types,
    get_all_reagents_by_type,
    count_reagents_by_type,
    is_reagent_in_database,
    filter_precedents_by_database_availability,
)

# Taxonomy management classes
from .taxonomy_store import (
    TaxonomyStore,
    RoleHeuristics,
    build_embedding_text,
)

# Constants
from .constants import (
    ROLE_FILES,
    DEFAULT_FAMILY_BY_ROLE,
    ROLE_PRIORITY,
    DEFAULT_TAXONOMY_DIR,
    DEFAULT_REGISTRY_DIR,
)

# Utility functions
from .taxonomy_utils import (
    normalize_cas,
    resolve_identity_from_cas,
    dedupe_synonyms,
    tokenize_all,
    build_entry,
)

# Validation functions
from .validator import (
    validate_entry,
    validate_role_file,
    validate_database,
    print_validation_summary,
    print_critical_errors_summary,
)

# Analytics functions
from .analytics import (
    get_database_statistics,
    get_family_statistics,
    find_reagents_by_family,
    get_missing_data_report,
    print_database_summary,
    print_role_summary,
    print_missing_data_report,
)

# Reactant/reaction type taxonomy utilities
from .types import (
    ReactantMatch,
    get_reactant_type_definitions,
    get_reactant_types_file,
    clear_reactant_type_cache,
    build_reactant_lookup,
    normalize_reactant_identifier,
    reactant_category_for,
    classify_reactant_smiles,
    classify_reactant_category,
    classify_reactant_group,
    classify_reactant_batch,
    get_reactant_category_matches,
    get_all_reactant_matches,
    get_reaction_type_definitions,
    get_reaction_types_file,
    clear_reaction_type_cache,
    list_reaction_type_ids,
    build_reaction_lookup,
    describe_reaction_type,
    normalize_reaction_type,
    iter_reactions_for_category,
    required_reactant_categories,
)

__all__ = [
    # Runtime lookup
    "find_reagent",
    "enrich_reagent_info", 
    "enrich_conditions",
    "load_reagent_database",
    "normalize_name",
    "get_all_reagent_types",
    "get_all_reagents_by_type",
    "count_reagents_by_type",
    "is_reagent_in_database",
    "filter_precedents_by_database_availability",
    # Taxonomy management
    "TaxonomyStore",
    "RoleHeuristics",
    "build_embedding_text",
    # Constants
    "ROLE_FILES",
    "DEFAULT_FAMILY_BY_ROLE",
    "ROLE_PRIORITY",
    "DEFAULT_TAXONOMY_DIR",
    "DEFAULT_REGISTRY_DIR",
    # Utilities
    "normalize_cas",
    "resolve_identity_from_cas",
    "dedupe_synonyms",
    "tokenize_all",
    "build_entry",
    # Validation
    "validate_entry",
    "validate_role_file",
    "validate_database",
    "print_validation_summary",
    "print_critical_errors_summary",
    # Analytics
    "get_database_statistics",
    "get_family_statistics",
    "find_reagents_by_family",
    "get_missing_data_report",
    "print_database_summary",
    "print_role_summary",
    "print_missing_data_report",
    # Reactant/reaction taxonomy
    "ReactantMatch",
    "get_reactant_type_definitions",
    "get_reactant_types_file",
    "clear_reactant_type_cache",
    "build_reactant_lookup",
    "normalize_reactant_identifier",
    "reactant_category_for",
    "classify_reactant_smiles",
    "classify_reactant_category",
    "classify_reactant_group",
    "classify_reactant_batch",
    "get_reactant_category_matches",
    "get_all_reactant_matches",
    "get_reaction_type_definitions",
    "get_reaction_types_file",
    "clear_reaction_type_cache",
    "list_reaction_type_ids",
    "build_reaction_lookup",
    "describe_reaction_type",
    "normalize_reaction_type",
    "iter_reactions_for_category",
    "required_reactant_categories",
]
