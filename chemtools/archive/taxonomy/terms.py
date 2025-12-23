"""
Chemistry term evaluation.

This module adds a maintainable layer for "fuzzy" organic chemistry language
such as "electron-poor aldehyde" or "sterically bulky carbonyl".

Terms are stored in `chemtools/archive/taxonomy/data/chem_terms.json` and expressed as
boolean rules over existing feature tokens (primarily from
`calculable_features.json`).
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, Iterable, Mapping, Optional

from . import load_registry
from .models import ChemTerm
from ..util.boolean_expr import evaluate as evaluate_rule


@lru_cache(maxsize=1)
def _load_terms() -> Dict[str, ChemTerm]:
    registry = load_registry()
    return dict(registry.chem_terms)


def get_terms() -> Dict[str, ChemTerm]:
    """Return the loaded chem terms keyed by term id."""
    return dict(_load_terms())


def evaluate_terms(
    features: Mapping[str, Any],
    *,
    term_ids: Optional[Iterable[str]] = None,
) -> Dict[str, bool]:
    """
    Evaluate chemistry terms against a mapping of feature tokens.

    Args:
        features: Mapping of feature token -> value (bool/int/etc).
        term_ids: Optional iterable restricting which term ids to evaluate.

    Returns:
        Dict of term_id -> boolean match.
    """
    terms = _load_terms()
    selected = list(term_ids) if term_ids is not None else list(terms.keys())

    results: Dict[str, bool] = {}
    for term_id in selected:
        term = terms.get(term_id)
        if term is None:
            results[term_id] = False
            continue
        rule = (term.rule or "").strip()
        results[term_id] = evaluate_rule(rule, features) if rule else False
    return results


def evaluate_terms_from_smiles(
    smiles: str,
    *,
    term_ids: Optional[Iterable[str]] = None,
) -> Dict[str, bool]:
    """
    Convenience helper: compute calculable features, then evaluate terms.

    This imports `chemtools.featurizers.calculable` lazily to avoid creating a
    hard dependency from the taxonomy subsystem back into feature extraction.
    """
    from ..featurizers import calculable

    features = calculable.detect_all_features(smiles)
    return evaluate_terms(features, term_ids=term_ids)
