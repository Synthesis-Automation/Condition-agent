"""
chemtools.retro — Retrosynthetic analysis engine.

Provides SMARTS-based retron detection and disconnection generation
for LLM-augmented retrosynthesis in ChemCoworker.

Public API:
    from chemtools.retro.retron_patterns import RETRONS
    from chemtools.retro.disconnector import rank_disconnections, find_retrons
"""
from __future__ import annotations
