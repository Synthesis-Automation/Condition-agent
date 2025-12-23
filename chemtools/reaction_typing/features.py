"""
Feature extraction helpers for reaction typing.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, Iterable, List, Mapping

from chemtools.util.smarts_cache import compile_smarts

from .motif_index import type_molecule

_FG_SMARTS = {
    "amide": "C(=O)N",
    "acid_or_ester": "C(=O)O",
    "aldehyde": "[#6X3H1](=O)",
}


@lru_cache(maxsize=1)
def _fg_queries() -> Dict[str, Any]:
    queries: Dict[str, Any] = {}
    for key, smarts in _FG_SMARTS.items():
        queries[key] = compile_smarts(smarts, validate=False)
    return queries


def count_elements(mol: Any) -> Dict[str, int]:
    """
    Count element symbols in a molecule.
    """
    counts: Dict[str, int] = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
    return counts


def count_functional_groups(mol: Any) -> Dict[str, int]:
    """
    Count functional group SMARTS matches in a molecule.
    """
    counts: Dict[str, int] = {}
    for key, query in _fg_queries().items():
        if query is None:
            continue
        counts[key] = len(mol.GetSubstructMatches(query, uniquify=True))
    return counts


def aggregate_side(
    mols: Iterable[Any],
    query_index: Mapping[str, Iterable[Any]],
) -> Dict[str, Dict[str, int]]:
    """
    Aggregate motif, element, and FG counts for a list of molecules.
    """
    motif_counts: Dict[str, int] = {}
    element_counts: Dict[str, int] = {}
    fg_counts: Dict[str, int] = {}

    for mol in mols:
        for hit in type_molecule(mol, query_index):
            motif_counts[hit] = motif_counts.get(hit, 0) + 1

        for key, value in count_elements(mol).items():
            element_counts[key] = element_counts.get(key, 0) + value

        for key, value in count_functional_groups(mol).items():
            fg_counts[key] = fg_counts.get(key, 0) + value

    return {
        "motifs": motif_counts,
        "elements": element_counts,
        "fg": fg_counts,
    }


def delta_counts(prod: Dict[str, Dict[str, int]], reac: Dict[str, Dict[str, int]]) -> Dict[str, Dict[str, int]]:
    """
    Compute product-reactant deltas for elements and functional groups.
    """
    out: Dict[str, Dict[str, int]] = {}
    for domain in ("elements", "fg"):
        delta: Dict[str, int] = {}
        keys = set(reac.get(domain, {})) | set(prod.get(domain, {}))
        for key in keys:
            delta[key] = prod.get(domain, {}).get(key, 0) - reac.get(domain, {}).get(key, 0)
        out[domain] = delta
    return out
