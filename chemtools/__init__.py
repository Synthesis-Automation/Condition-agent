"""Chemtools package.

Keep package import light. Submodules are available via explicit imports, e.g.:

    from chemtools import registry
    from chemtools import featurizers

This avoids importing heavy optional dependencies (e.g., RDKit) unless needed.
"""

__all__ = [
    "smiles",
    "router",
    "featurizers",
    "condition_core",
    "properties",
    "precedent",
    "constraints",
    "explain",
    "registry",
]
