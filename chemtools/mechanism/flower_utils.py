"""
Adapted lightweight helpers from the FlowER project for electron bookkeeping.

Source: Fong Mun Hong et al., FlowER (MIT License).
We vendor a minimal subset (bond-electron matrix utilities) to avoid the heavy
training/inference stack while benefiting from their RDKit-based heuristics.
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


BT_TO_ELECTRON = {
    rdchem.BondType.SINGLE: 2,
    rdchem.BondType.DOUBLE: 4,
    rdchem.BondType.TRIPLE: 6,
    rdchem.BondType.AROMATIC: 3,
}

_PERIODIC_TABLE = Chem.GetPeriodicTable()
_SMI_PARAMS = Chem.SmilesParserParams()
_SMI_PARAMS.removeHs = False
_SMI_PARAMS.sanitize = True


def _count_lone_pairs(atom: rdchem.Atom) -> int:
    valence = _PERIODIC_TABLE.GetNOuterElecs(atom.GetAtomicNum())
    charge = atom.GetFormalCharge()
    bonds = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    hydrogens = atom.GetTotalNumHs()
    return int(valence - charge - bonds - hydrogens)


def get_be_matrix(smiles: str) -> np.ndarray:
    """
    Return the bond-electron matrix used by FlowER for an atom-mapped SMILES.

    Diagonal entries store lone-pair counts; off-diagonals store half of the
    bonding electrons between atom pairs so that the matrix sums to zero.
    """

    mol = Chem.MolFromSmiles(smiles, _SMI_PARAMS)
    if mol is None:
        raise ValueError("RDKit failed to parse SMILES for BE matrix.")
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # keep aromatic flags if kekulization fails
        pass

    atoms = list(mol.GetAtoms())
    if not atoms or not atoms[0].HasProp("molAtomMapNumber"):
        raise ValueError("Atom-mapped SMILES required for BE matrix.")

    size = len(atoms)
    matrix = np.zeros((size, size))

    for atom in atoms:
        idx = atom.GetIntProp("molAtomMapNumber") - 1
        matrix[idx, idx] = _count_lone_pairs(atom)

    for bond in mol.GetBonds():
        start = bond.GetBeginAtom().GetIntProp("molAtomMapNumber") - 1
        end = bond.GetEndAtom().GetIntProp("molAtomMapNumber") - 1
        electrons = BT_TO_ELECTRON.get(bond.GetBondType())
        if electrons is None:
            continue
        matrix[start, end] = matrix[end, start] = electrons / 2.0

    return matrix


def _matrix_to_atom_bond_dicts(matrix: np.ndarray) -> Tuple[Dict[int, float], Dict[Tuple[int, int], float]]:
    diag = matrix.diagonal()
    atom_dict = {idx + 1: float(value) for idx, value in enumerate(diag) if value != 0}

    iu, ju = np.triu_indices_from(matrix, k=1)
    bond_values = matrix[iu, ju] + matrix[ju, iu]
    bond_dict = {
        (int(i + 1), int(j + 1)): float(val)
        for i, j, val in zip(iu, ju, bond_values)
        if val != 0
    }
    return atom_dict, bond_dict


def compute_electron_balance(reaction_smiles: str) -> Dict[str, List[Dict[str, float]]]:
    """
    Compute per-atom electron differences between reactants and products.

    Returns a dictionary with atom and bond deltas suitable for evidence logs.
    """

    try:
        reactants_section, _, products_section = _split_reaction_sections(reaction_smiles)
        reactant_matrix = get_be_matrix(reactants_section)
        product_matrix = get_be_matrix(products_section)
    except Exception as exc:
        raise ValueError(f"Unable to compute electron balance: {exc}") from exc

    react_atoms, react_bonds = _matrix_to_atom_bond_dicts(reactant_matrix)
    prod_atoms, prod_bonds = _matrix_to_atom_bond_dicts(product_matrix)

    atom_keys = sorted(set(react_atoms) | set(prod_atoms))
    bond_keys = sorted(set(react_bonds) | set(prod_bonds))

    atom_balance = []
    for key in atom_keys:
        delta = prod_atoms.get(key, 0.0) - react_atoms.get(key, 0.0)
        if delta != 0:
            atom_balance.append({"atom": int(key), "delta_lone_pairs": float(delta)})

    bond_balance = []
    for key in bond_keys:
        delta = prod_bonds.get(key, 0.0) - react_bonds.get(key, 0.0)
        if delta != 0:
            bond_balance.append({"bond": list(map(int, key)), "delta_electrons": float(delta)})

    return {
        "atom_balance": atom_balance,
        "bond_balance": bond_balance,
    }


def _split_reaction_sections(reaction_smiles: str) -> Tuple[str, str, str]:
    """Split a reaction SMILES into reactants, agents, products sections."""
    if ">>" not in reaction_smiles:
        raise ValueError("Reaction SMILES must contain '>>'.")
    parts = reaction_smiles.split(">")
    if len(parts) == 2:
        reactants, products = parts[0], parts[1]
        agents = ""
    else:
        reactants, agents, products = parts[0], parts[1], parts[2]
    return reactants.strip(), agents.strip(), products.strip()


__all__ = ["get_be_matrix", "compute_electron_balance"]
