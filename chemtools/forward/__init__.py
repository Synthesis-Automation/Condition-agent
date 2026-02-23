"""
chemtools.forward — Forward synthesis prediction engine.

Given two reactant SMILES (or a single reactant for unimolecular reactions),
enumerate possible products and rank them by yield likelihood.

Mirrors the chemtools.retro pipeline in the opposite direction:

    retro:   target_product  → [retro SMARTS] → precursor pair
    forward: reactant pair   → [fwd SMARTS]   → product candidates

Public API
----------
    from chemtools.forward import ReactantAnalyzer, ForwardReactor, score_products

    analyzer = ReactantAnalyzer()
    analysis = analyzer.analyze_pair("Brc1ccccc1", "Nc1ccccc1")
    # → {"compatible_reactions": [...], "fg_a": [...], "fg_b": [...]}

    reactor   = ForwardReactor()
    products  = reactor.generate("Brc1ccccc1", "Nc1ccccc1")
    # → List[ProductPrediction]

    ranked    = score_products(products, smiles_a="Brc1ccccc1", smiles_b="Nc1ccccc1")
    # → List[ProductPrediction] sorted by overall_score desc
"""
from __future__ import annotations

from .contracts import ForwardTemplateMatch, ProductPrediction, ReactantProfile
from .reactor import ForwardReactor, ReactantAnalyzer
from .scoring import score_products

__all__ = [
    "ForwardTemplateMatch",
    "ProductPrediction",
    "ReactantProfile",
    "ForwardReactor",
    "ReactantAnalyzer",
    "score_products",
]
