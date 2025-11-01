"""
Constants for reagent taxonomy system.

Defines role mappings, priority, and family patterns used throughout the
reagent management system.
"""

from typing import Dict, Tuple

# Mapping of reagent roles to their taxonomy JSON files
ROLE_FILES: Dict[str, str] = {
    "acid": "taxonomy_acid.json",
    "additive": "taxonomy_additive.json",
    "ligand": "taxonomy_ligand.json",
    "metal_catalyst": "taxonomy_catalysts_precursor.json",
    "base": "taxonomy_base.json",
    "condensation_agent": "taxonomy_condensation_agent.json",
    "oxidant": "taxonomy_oxidant.json",
    "reductant": "taxonomy_reductant.json",
    "solvent": "taxonomy_solvent.json",
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

