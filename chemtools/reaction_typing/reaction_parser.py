"""
Reaction SMILES parsing utilities.
"""

from __future__ import annotations

import re
from typing import List, Tuple

from chemtools.util.rdkit_helpers import parse_smiles

_LABEL_RE = re.compile(r"\s*\([^)]*\)\s*$")


def strip_label(rxn_line: str) -> str:
    """
    Remove a trailing "(...)" label from a reaction line.
    """
    return _LABEL_RE.sub("", rxn_line or "").strip()


def parse_rxn(rxn_line: str) -> Tuple[str, str]:
    """
    Parse "left>>right" reaction strings (labels optional).
    """
    rxn_line = strip_label(rxn_line)
    if ">>" not in rxn_line:
        raise ValueError("missing_reaction_arrow")
    left, right = rxn_line.split(">>", 1)
    return left.strip(), right.strip()


def split_mols(side_smiles: str) -> List[str]:
    """
    Split dot-separated SMILES into individual molecules.
    """
    side_smiles = (side_smiles or "").strip()
    if not side_smiles:
        return []
    return [item for item in side_smiles.split(".") if item]


def mols_from_smiles_list(smiles_list: List[str]):
    """
    Convert a list of SMILES to RDKit mols; return (mols, invalid_smiles).
    """
    mols = []
    invalid = []
    for smiles in smiles_list:
        mol = parse_smiles(smiles)
        if mol is None:
            invalid.append(smiles)
            continue
        mols.append(mol)
    return mols, invalid
