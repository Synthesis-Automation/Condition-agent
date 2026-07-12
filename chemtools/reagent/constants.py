"""
Constants for reagent taxonomy system.

Defines role mappings, priority, and family patterns used throughout the
reagent management system.
"""

from typing import Dict, Tuple

# Registry CSV filename (single source for all roles)
REGISTRY_CSV = "substances.v1.csv"

# Mapping of reagent roles to the shared registry CSV
ROLE_FILES: Dict[str, str] = {
    "acid": REGISTRY_CSV,
    "additive": REGISTRY_CSV,
    "ligand": REGISTRY_CSV,
    "metal_catalyst": REGISTRY_CSV,
    "base": REGISTRY_CSV,
    "condensation_agent": REGISTRY_CSV,
    "oxidant": REGISTRY_CSV,
    "reductant": REGISTRY_CSV,
    "solvent": REGISTRY_CSV,
    "other_reagent": REGISTRY_CSV,
    "enzyme": REGISTRY_CSV,
    "organo_catalyst": REGISTRY_CSV,
}

# Default family assignment when role is inferred but family cannot be determined
DEFAULT_FAMILY_BY_ROLE: Dict[str, str] = {
    "acid": "mineral_acids",
    "additive": "quaternary_ammonium_ptc",
    "ligand": "trialkyl_triaryl_phosphines",
    "metal_catalyst": "pd_ii_salts",
    "base": "tertiary_amines_aliphatic",
    "condensation_agent": "carbodiimides",
    "oxidant": "o2_gas",
    "reductant": "metal_powders",
    "solvent": "hydrocarbons_aromatic",
    "other_reagent": "misc_general",
}

# Priority order for role assignment (lower number = higher priority)
ROLE_PRIORITY: Dict[str, int] = {
    "ligand": 0,
    "metal_catalyst": 1,
    "base": 2,
    "acid": 3,
    "condensation_agent": 4,
    "oxidant": 5,
    "reductant": 6,
    "additive": 7,
    "solvent": 8,
    "other_reagent": 9,
}

# Taxonomy directory path relative to data folder
DEFAULT_TAXONOMY_DIR = "compound_taxonomy"

# Registry directory path relative to data folder  
DEFAULT_REGISTRY_DIR = "condition_registry/definitions"

# Legacy role aliases that should normalize to the canonical role keys
ROLE_ALIASES: Dict[str, str] = {
    "metal_precursor": "metal_catalyst",
    "preformed_metal_catalyst": "metal_catalyst",
}

