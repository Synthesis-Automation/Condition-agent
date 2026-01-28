"""Backward compatibility shim for motif detection.

.. deprecated::
    This module is a compatibility shim. Import from ``chemtools.featurizers.motifs``
    subpackage instead (split into detection.py, models.py).

All functionality has been moved to the motifs/ subpackage:
- motifs.detection: Motif detection functions (detect_motifs, discover_undocumented_motifs)
- motifs.models: CompoundPattern dataclass and constants
"""

# Re-export everything from the new subpackage for backward compatibility
from .motifs.detection import *  # noqa: F401,F403
from .motifs.models import CompoundPattern, _DISCOVERY_SKIP_SUBSTITUENTS  # noqa: F401

__all__ = [
    "detect_motifs",
    "discover_undocumented_motifs",
    "CompoundPattern",
    "_DISCOVERY_SKIP_SUBSTITUENTS",
]
