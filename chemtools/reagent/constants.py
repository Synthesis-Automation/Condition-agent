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
    "metal_precursor": "taxonomy_catalysts_precursor.json",
    "base": "taxonomy_base.json",
    "coupling_reagent": "taxonomy_coupling_reagent.json",
    "oxidant": "taxonomy_oxidant.json",
    "reductant": "taxonomy_reductant.json",
    "solvent": "taxonomy_solvent.json",
}

# Default family assignment when role is inferred but family cannot be determined
DEFAULT_FAMILY_BY_ROLE: Dict[str, str] = {
    "acid": "mineral_acids",
    "additive": "quaternary_ammonium_ptc",
    "ligand": "trialkyl_triaryl_phosphines",
    "metal_precursor": "pd_ii_salts",
    "base": "tertiary_amines_aliphatic",
    "coupling_reagent": "carbodiimides",
    "oxidant": "o2_gas",
    "reductant": "metal_powders",
    "solvent": "hydrocarbons_aromatic",
    "other_reagent": "misc_general",
}

# Priority order for role assignment (lower number = higher priority)
ROLE_PRIORITY: Dict[str, int] = {
    "ligand": 0,
    "metal_precursor": 1,
    "base": 2,
    "acid": 3,
    "coupling_reagent": 4,
    "oxidant": 5,
    "reductant": 6,
    "additive": 7,
    "solvent": 8,
    "other_reagent": 9,
}

# Taxonomy directory path relative to data folder
DEFAULT_TAXONOMY_DIR = "compound_taxonomy"

# Registry directory path relative to data folder  
DEFAULT_REGISTRY_DIR = "reagents"
