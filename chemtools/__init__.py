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
from .detection import detect_reaction_type, detect_reaction_types, DetectionResult, ReactionMatch
from .reaction_inference import (
    analyze_reaction_general,
    GeneralReactionAnalysis,
    ReactionDecision,
    ReactionValidation,
)
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
try:
    from .mechanism import classify_mechanism_simple
    from .mechanism.electron_flow import predict_electron_flow
    from .mechanism.intermediates import predict_intermediates
    _MECHANISM_AVAILABLE = True
except Exception:  # pragma: no cover - optional dependency
    _MECHANISM_AVAILABLE = False

    def _missing_mechanism(name: str):
        raise ImportError(
            f"{name} unavailable because chemtools.mechanism is not installed."
        )

    def classify_mechanism_simple(*args, **kwargs):
        return _missing_mechanism("classify_mechanism_simple")

    def predict_electron_flow(*args, **kwargs):
        return _missing_mechanism("predict_electron_flow")

    def predict_intermediates(*args, **kwargs):
        return _missing_mechanism("predict_intermediates")
from .visualization import render_molecule_image, render_reaction_image

# Primary exports - ChemTools is the recommended interface
__all__ = [
    "chem",          # Global singleton instance (recommended)
    "ChemTools",     # Class for custom instances
    "ResourceConfig", # Configuration dataclass
    # Reaction type detection (NEW - clean API)
    "detect_reaction_type",  # Main detection function
    "detect_reaction_types", # Alias for detect_reaction_type
    "DetectionResult",       # Detection result dataclass
    "ReactionMatch",         # Individual match dataclass
    "analyze_reaction_general",  # General taxonomy-first reaction inference
    "GeneralReactionAnalysis",
    "ReactionDecision",
    "ReactionValidation",
    # Atom mapping and bond analysis
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
