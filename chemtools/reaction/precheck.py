"""
Precheck helpers for reaction-level featurization.

These utilities keep lightweight, deterministic normalization heuristics out of
`reaction.py` so the main formatter stays easier to read and extend.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Iterable, List, Optional, Tuple

from chemtools.core import rdkit as rdkit_helpers

_HALOGEN_ELEMENTS = {"F", "Cl", "Br", "I"}


def _count_elements(smiles_list: Iterable[str]) -> Counter[str]:
    counts: Counter[str] = Counter()
    if not rdkit_helpers.rdkit_available():
        return counts
    for smiles in smiles_list or []:
        mol = rdkit_helpers.parse_smiles(smiles)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts


def infer_reactant_repeats_from_stoichiometry(
    reactant_smiles: List[str],
    product_smiles: List[str],
) -> Dict[str, Any]:
    """
    Conservatively infer missing repeated reactants from atom-count imbalance.

    This heuristic only applies when one unique reactant cleanly explains the
    non-halogen atom surplus on the product side.
    """
    if not reactant_smiles or not product_smiles or not rdkit_helpers.rdkit_available():
        return {"applied": False}

    reactant_counts = _count_elements(reactant_smiles)
    product_counts = _count_elements(product_smiles)
    if not reactant_counts or not product_counts:
        return {"applied": False}

    delta: Counter[str] = Counter()
    for element in set(reactant_counts) | set(product_counts):
        delta[element] = int(product_counts.get(element, 0)) - int(reactant_counts.get(element, 0))

    positive_non_halogen = {
        el: cnt
        for el, cnt in delta.items()
        if cnt > 0 and el not in _HALOGEN_ELEMENTS and el != "H"
    }
    if not positive_non_halogen:
        return {"applied": False}

    # Keep the heuristic conservative: do not attempt correction when there are
    # deficits on non-halogen heavy atoms.
    if any(
        cnt < 0 and el not in _HALOGEN_ELEMENTS and el != "H"
        for el, cnt in delta.items()
    ):
        return {"applied": False, "reason": "non_halogen_deficit_present"}

    if any(delta.get(hal, 0) > 0 for hal in _HALOGEN_ELEMENTS):
        return {"applied": False, "reason": "positive_halogen_delta_present"}

    candidates: List[Tuple[int, int]] = []
    for idx, smiles in enumerate(reactant_smiles):
        mol_counts = _count_elements([smiles])
        if not mol_counts:
            continue
        if sum(mol_counts.get(hal, 0) for hal in _HALOGEN_ELEMENTS) <= 0:
            continue

        non_halogen_formula = {
            el: cnt
            for el, cnt in mol_counts.items()
            if cnt > 0 and el not in _HALOGEN_ELEMENTS and el != "H"
        }
        if not non_halogen_formula:
            continue
        if set(non_halogen_formula) != set(positive_non_halogen):
            continue

        repeat_count: Optional[int] = None
        valid = True
        for element, atom_count in non_halogen_formula.items():
            surplus = positive_non_halogen.get(element, 0)
            if atom_count <= 0 or surplus <= 0 or surplus % atom_count != 0:
                valid = False
                break
            multiplier = surplus // atom_count
            if multiplier <= 0:
                valid = False
                break
            if repeat_count is None:
                repeat_count = int(multiplier)
            elif repeat_count != int(multiplier):
                valid = False
                break
        if not valid or repeat_count is None:
            continue

        # Conservative upper bound to avoid over-correction.
        if repeat_count > 2:
            continue
        candidates.append((idx, int(repeat_count)))

    if len(candidates) != 1:
        if not candidates:
            return {"applied": False}
        return {"applied": False, "reason": "ambiguous_repeat_candidates"}

    idx, repeat_count = candidates[0]
    if repeat_count <= 0:
        return {"applied": False}

    return {
        "applied": True,
        "reactant_index": idx,
        "reactant_smiles": reactant_smiles[idx],
        "repeat_count": repeat_count,
        "reason": "stoichiometry_non_halogen_surplus_match",
        "delta": dict(delta),
    }


__all__ = ["infer_reactant_repeats_from_stoichiometry"]

