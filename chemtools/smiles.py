"""
Backwards-compatible wrapper for the SMILES helpers now housed under
``chemtools.featurizers.analysis.smiles``.
"""

from __future__ import annotations

from .featurizers.analysis.smiles import normalize, normalize_reaction

__all__ = ["normalize", "normalize_reaction"]
