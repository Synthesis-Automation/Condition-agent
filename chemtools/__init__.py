"""
ChemTools - Chemistry Toolbox for Automated Synthesis

A unified, object-oriented API for chemistry operations with intelligent
resource management, lazy loading, and performance optimization.

Primary Interface:
    from chemtools import chem
    
    # SMILES operations
    result = chem.smiles.normalize("CCO")
    
    # Precedent search
    precedents = chem.precedent.knn("C_N_Coupling", features={...}, k=5)
    
    # ML recommendations
    recommendations = chem.recommend.conditions(reaction="...", k=5)

Custom Instances:
    from chemtools import ChemTools
    
    # For API servers - preload everything
    api_chem = ChemTools(
        datasets=["C_N_Coupling"],
        preload=True
    )
    
    # For CLI tools - minimal config
    cli_chem = ChemTools(
        datasets=None,
        preload=False
    )

Performance:
    - 50-100x faster with selective loading
    - 70% less memory with targeted datasets
    - Thread-safe resource management
    - Intelligent caching with lazy loading

See CHEMTOOLS_QUICKSTART.md for more examples.
"""

from .context import ChemTools, chem, ResourceConfig
from .detection import detect_reaction
from ._atom_mapping import (
    add_atom_mapping,
    analyze_bond_changes,
    analyze_with_mcs,
    analyze_bond_changes_hybrid,
    is_available as rxnmapper_available,
)
from .util.reaction_center_detector import (
    identify_changed_atoms_from_mapped_smiles,
    compare_unmapped_reaction_to_find_changes,
)
from .mechanism import classify_mechanism_simple
from .mechanism.electron_flow import predict_electron_flow
from .mechanism.intermediates import predict_intermediates
from .visualization import render_molecule_image, render_reaction_image

# Primary exports - ChemTools is the recommended interface
__all__ = [
    "chem",          # Global singleton instance (recommended)
    "ChemTools",     # Class for custom instances
    "ResourceConfig", # Configuration dataclass
    "detect_reaction", # Unified reaction detection API
    # Atom mapping and bond analysis (NEW)
    "add_atom_mapping",  # Add atom mapping to reaction SMILES
    "analyze_bond_changes",  # High-level bond analysis (RXNMapper)
    "analyze_with_mcs",  # MCS-based bond change estimation
    "analyze_bond_changes_hybrid",  # Hybrid RXNMapper + MCS (recommended)
    "identify_changed_atoms_from_mapped_smiles",  # Low-level bond analysis
    "compare_unmapped_reaction_to_find_changes",  # Auto-map then analyze
    "rxnmapper_available",  # Check if RXNMapper is installed
    "classify_mechanism_simple",  # Lightweight mechanism classifier
    "predict_electron_flow",  # Rule-based electron flow heuristics
    "predict_intermediates",  # Mechanistic intermediate heuristics
    "render_molecule_image",  # Molecule drawing helper (RDKit)
    "render_reaction_image",  # Reaction drawing helper (RDKit)
]

# Version info
__version__ = "2.0.0"
__author__ = "Synthesis Automation Team"
