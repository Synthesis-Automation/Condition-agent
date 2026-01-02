"""
Compatibility shim for molecule and reaction featurizers.
"""

from __future__ import annotations

from .molecule import analyze_smiles, featurize_molecule
from .reaction import featurize_reaction

__all__ = ["featurize_molecule", "featurize_reaction", "analyze_smiles"]
