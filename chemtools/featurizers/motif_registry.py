"""Backward compatibility shim for motif registry.

.. deprecated::
    This module is a compatibility shim. Import from ``chemtools.featurizers.motifs``
    subpackage instead (split into models.py, registry.py, classification.py, utils.py).

All functionality has been moved to the motifs/ subpackage:
- motifs.models: CompoundPattern dataclass, constants
- motifs.registry: Registry loading and compilation
- motifs.classification: Compound classification functions
- motifs.utils: SMARTS utilities and validation
"""

# Re-export everything from the new subpackage for backward compatibility
from .motifs import *  # noqa: F401,F403

# Maintain the same public API
from .motifs import __all__  # noqa: F401
