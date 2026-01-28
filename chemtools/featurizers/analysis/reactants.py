"""Backward compatibility shim for reactant classification and reaction type utilities.

.. deprecated::
    This module is a compatibility shim. Import from ``chemtools.featurizers.analysis.reactants``
    subpackage instead (split into taxonomy.py, classification.py, reactions.py).

All functionality has been moved to the reactants/ subpackage:
- reactants.taxonomy: Reactant type definitions and alias resolution
- reactants.classification: SMARTS-based reactant classification  
- reactants.reactions: Reaction type definitions and lookups
"""

# Re-export everything from the new subpackage for backward compatibility
from .reactants import *  # noqa: F401,F403

# Maintain the same public API
from .reactants import __all__  # noqa: F401
