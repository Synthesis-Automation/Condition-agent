"""
chemtools.retro — Retrosynthetic analysis engine.

Provides SMARTS-based retron detection and disconnection generation
for LLM-augmented retrosynthesis in ChemCoworker.

Public API:
    from chemtools.retro.retron_patterns import RETRONS
    from chemtools.retro.disconnector import rank_disconnections, find_retrons

Template extraction & integration:
    from chemtools.retro.extract_templates import (
        load_reactions_from_csv,
        extract_templates_from_reactions,
        quality_score,
    )
    from chemtools.retro.integrate_templates import integrate
"""
from __future__ import annotations
