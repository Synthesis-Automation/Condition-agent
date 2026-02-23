"""
Product scoring for forward synthesis predictions.

score_products() takes a list of ProductPrediction objects and enriches them
with scoring components, then returns them sorted by overall_score descending.

Scoring components (mirror the retrosynthesis scorer):
  retro:   complexity_delta + fragment_balance + difficulty
  forward: hte_yield_proxy + chemoselectivity_penalty + structural_complexity

Design
------
- hte_yield_proxy   : avg yield for this reaction class in the HTE database
- chemoselectivity_penalty : if a reactant matches multiple templates,
                              lower-ranked templates are penalised
- structural_complexity : BertzCT of the predicted product (informational,
                          high = higher-value target but not necessarily better)
- overall_score = normalize(yield_proxy) - chemoselectivity_penalty
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

from .contracts import ProductPrediction

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HTE yield proxy: cached avg yield per HTE family
# ---------------------------------------------------------------------------

_YIELD_CACHE: Dict[str, float] = {}


def _get_avg_yield(hte_families: List[str]) -> float:
    """Return the average yield (0–100) for the first resolvable HTE family."""
    if not hte_families:
        return 50.0   # neutral prior

    for family in hte_families:
        if family in _YIELD_CACHE:
            return _YIELD_CACHE[family]
        try:
            from chem_coworker.tools.retrosynthesis import _fast_load_hte_family_cached
            rows = list(_fast_load_hte_family_cached(family))
            if rows:
                yields = [r.get("yield_value") or 0.0 for r in rows]
                valid_yields = [y for y in yields if y > 0]
                avg = sum(valid_yields) / len(valid_yields) if valid_yields else 50.0
                avg = round(avg, 1)
                _YIELD_CACHE[family] = avg
                return avg
        except Exception:
            pass

    return 50.0


def _bertz_ct(smiles: str) -> float:
    """BertzCT complexity of a molecule."""
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem.GraphDescriptors import BertzCT
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
        return round(BertzCT(mol), 1)
    except Exception:
        return 0.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def score_products(
    predictions: List[ProductPrediction],
    smiles_a: str = "",
    smiles_b: str = "",
) -> List[ProductPrediction]:
    """Score and rank a list of ProductPrediction objects.

    Enriches each prediction with:
    - hte_yield_proxy       (from HTE database avg yield for the reaction family)
    - structural_complexity (BertzCT of predicted product)
    - chemoselectivity_penalty (if multiple templates fire for the same reactants)
    - confidence_label      ("high" / "medium" / "low")
    - competing_templates   (names of other matching templates)
    - overall_score

    Returns predictions sorted by overall_score descending.

    Args:
        predictions: List of ProductPrediction objects from ForwardReactor.generate().
        smiles_a: First reactant SMILES (used for competing template detection).
        smiles_b: Second reactant SMILES.
    """
    if not predictions:
        return []

    # ── 1. HTE yield proxy ─────────────────────────────────────────────────
    for pred in predictions:
        pred.hte_yield_proxy = _get_avg_yield(pred.hte_families)

    # ── 2. Structural complexity ────────────────────────────────────────────
    for pred in predictions:
        pred.structural_complexity = _bertz_ct(pred.product_smiles)

    # ── 3. Chemoselectivity penalty ─────────────────────────────────────────
    # If multiple templates fired for the same reactant pair, penalise all but
    # the highest-yield one.
    template_names = [p.template_name for p in predictions]
    for pred in predictions:
        pred.competing_templates = [n for n in template_names if n != pred.template_name]
        # Penalty: 5 points per competing template (max 25)
        n_competitors = len(pred.competing_templates)
        pred.chemoselectivity_penalty = min(25.0, n_competitors * 5.0)

    # ── 4. Difficulty adjustment ────────────────────────────────────────────
    # Easier reactions (lower difficulty) get a +10 bonus
    for pred in predictions:
        difficulty_bonus = (1.0 - pred.difficulty) * 10.0

        # ── 5. Overall score ────────────────────────────────────────────────
        # Scale yield_proxy to 0–100, then subtract penalties, add bonuses
        pred.overall_score = round(
            pred.hte_yield_proxy             # 0–100 base
            - pred.chemoselectivity_penalty  # –0 to –25
            + difficulty_bonus               # +0 to +10
        , 1)

    # ── 6. Confidence label ─────────────────────────────────────────────────
    for pred in predictions:
        score = pred.overall_score
        if score >= 70:
            pred.confidence_label = "high"
        elif score >= 45:
            pred.confidence_label = "medium"
        else:
            pred.confidence_label = "low"

    # ── 7. Sort descending ──────────────────────────────────────────────────
    predictions.sort(key=lambda p: p.overall_score, reverse=True)
    return predictions
