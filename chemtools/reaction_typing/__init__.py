"""
Reaction typing via motif hits and handle deltas (taxonomy v2).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from .motif_index import build_query_index
from .reaction_parser import parse_rxn, split_mols, mols_from_smiles_list
from .rules import classify_reaction


def classify_reaction_smiles(
    rxn_line: str,
    *,
    query_index: Optional[Mapping[str, List[Any]]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
) -> Dict[str, Any]:
    """
    Classify a reaction SMILES using motif hits and handle deltas.
    """
    result: Dict[str, Any] = {
        "rxn": rxn_line,
        "ok": False,
        "predicted": "Unknown",
        "confidence": 0.0,
        "evidence": {},
    }
    try:
        left, right = parse_rxn(rxn_line)
    except ValueError as exc:
        result["error"] = str(exc)
        return result

    reactant_smiles = split_mols(left)
    product_smiles = split_mols(right)
    reactant_mols, invalid_reactants = mols_from_smiles_list(reactant_smiles)
    product_mols, invalid_products = mols_from_smiles_list(product_smiles)
    if invalid_reactants or invalid_products:
        result["error"] = "invalid_smiles"
        result["invalid_reactants"] = invalid_reactants
        result["invalid_products"] = invalid_products
        return result

    if query_index is None:
        query_index = build_query_index(registry_paths=registry_paths)

    classification = classify_reaction(
        reactant_mols,
        product_mols,
        query_index,
    )
    result.update(classification)
    result["ok"] = True
    return result


def classify_reaction_batch(
    rxn_lines: Iterable[str],
    *,
    query_index: Optional[Mapping[str, List[Any]]] = None,
    registry_paths: Optional[Mapping[str, str | Path]] = None,
) -> List[Dict[str, Any]]:
    """
    Classify a batch of reaction SMILES strings.
    """
    if query_index is None:
        query_index = build_query_index(registry_paths=registry_paths)
    return [
        classify_reaction_smiles(
            rxn_line,
            query_index=query_index,
            registry_paths=registry_paths,
        )
        for rxn_line in rxn_lines
    ]


__all__ = ["classify_reaction_smiles", "classify_reaction_batch", "build_query_index"]
