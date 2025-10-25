"""
Automatic Reactant Type Classifier using SMARTS patterns.

This module now delegates to :mod:`chemtools.reagent.types` so that downstream
HTE notebooks and scripts can continue to import ``classify_reactant`` while the
canonical implementation lives inside ChemTools.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

from chemtools.reagent import (
    ReactantMatch,
    classify_reactant_batch as _classify_batch,
    classify_reactant_category as _classify_category,
    classify_reactant_group as _classify_group,
    classify_reactant_smiles as _classify_single,
    get_all_reactant_matches as _get_all_matches,
    get_reactant_category_matches as _get_category_matches,
    get_reactant_type_definitions,
)


def load_reactant_types() -> Dict:
    """Return the SMARTS taxonomy used for classification."""
    return get_reactant_type_definitions()


def _match_to_dict(match: Optional[ReactantMatch]) -> Optional[Dict]:
    if match is None:
        return None
    return {
        "category": match.category,
        "member_type": match.member_type,
        "name": match.name,
        "group": match.group,
        "description": match.description,
        "smarts": match.smarts,
        "category_smarts": match.category_smarts,
        "specificity": match.specificity,
        "is_general": match.is_general,
    }


def classify_reactant(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[Dict]:
    """Classify a reactant based on SMILES using the shared ChemTools taxonomy."""
    match = _classify_single(smiles, reactant_types)
    return _match_to_dict(match)


def classify_batch(
    smiles_list: Iterable[str], reactant_types: Optional[Dict] = None
) -> List[Optional[Dict]]:
    """Batch classification convenience wrapper returning legacy dict payloads."""
    matches = _classify_batch(smiles_list, reactant_types)
    return [_match_to_dict(match) for match in matches]


def get_all_matches(smiles: str, reactant_types: Optional[Dict] = None) -> List[Dict]:
    """Return all SMARTS matches as dictionaries (legacy structure)."""
    matches = _get_all_matches(smiles, reactant_types)
    return [
        {
            "category": match.category,
            "member_type": match.member_type,
            "name": match.name,
            "group": match.group,
            "smarts": match.smarts,
        }
        for match in matches
    ]


def get_category_matches(smiles: str, reactant_types: Optional[Dict] = None) -> List[str]:
    """Return all matching category identifiers for the molecule."""
    return _get_category_matches(smiles, reactant_types)


def classify_by_category(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[str]:
    """Return only the category identifier for the best match."""
    return _classify_category(smiles, reactant_types)


def classify_by_group(smiles: str, reactant_types: Optional[Dict] = None) -> Optional[str]:
    """Return only the functional group name for the best match."""
    return _classify_group(smiles, reactant_types)


if __name__ == "__main__":  # pragma: no cover - manual smoke test retained
    import sys

    smiles_inputs = sys.argv[1:] or [
        "c1ccccc1Br",
        "CCBr",
        "c1cccnc1Br",
        "C=CCBr",
        "c1ccccc1B(O)O",
    ]

    rt = load_reactant_types()
    for smiles in smiles_inputs:
        result = classify_reactant(smiles, rt)
        print(f"{smiles:>25s} -> {result if result else 'No match'}")
