"""
Constants for reagent taxonomy system.

Defines role mappings, priority, and family patterns used throughout the
reagent management system.
"""

from typing import Dict, Tuple

# Mapping of reagent roles to their registry JSON files
ROLE_FILES: Dict[str, str] = {
    "acid": "acid.json",
    "additive": "additive.json",
    "ligand": "ligand.json",
    "metal_catalyst": "metal_catalyst.json",
    "base": "base.json",
    "condensation_agent": "condensation_agent.json",
    "oxidant": "oxidant.json",
    "reductant": "reductant.json",
    "solvent": "solvent.json",
    "other_reagent": "other_reagent.json",
    "enzyme": "enzyme.json",
    "organo_catalyst": "organo_catalyst.json",
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
DEFAULT_REGISTRY_DIR = "reagent_db"

# Legacy role aliases that should normalize to the canonical role keys
ROLE_ALIASES: Dict[str, str] = {
    "metal_precursor": "metal_catalyst",
    "preformed_metal_catalyst": "metal_catalyst",
}

