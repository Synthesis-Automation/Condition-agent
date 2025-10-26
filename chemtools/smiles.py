"""
Backwards-compatible wrapper for the SMILES helpers now housed under
``chemtools.analysis.smiles``.
"""

from __future__ import annotations

from .analysis.smiles import normalize, normalize_reaction

__all__ = ["normalize", "normalize_reaction"]

