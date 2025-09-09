"""
Molecular featurizer (general representation).

This module re-exports the existing featurize() from the legacy
`ullmann` featurizer to provide a more general entry point while
keeping backwards compatibility.
"""

from typing import Dict, Any

from .ullmann import featurize as featurize  # noqa: F401

__all__ = ["featurize"]

