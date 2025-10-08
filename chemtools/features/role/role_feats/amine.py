from __future__ import annotations

from typing import Any, List, Tuple

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles


def _amine_class(mol) -> int:
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return -1
    patt_tert = Chem.MolFromSmarts("[NX3]([#6])([#6])[#6]")
    patt_sec = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")
    patt_prim = Chem.MolFromSmarts("[NX3;H2][#6]")
    if mol.HasSubstructMatch(patt_tert):
        return 3
    if mol.HasSubstructMatch(patt_sec):
        return 2
    if mol.HasSubstructMatch(patt_prim):
        return 1
    return -1


def amine_features(mol_or_smiles: Any, center_atoms: List[int]) -> Tuple[List[float], List[str]]:
    fields = [
        "amine.present",
        "amine.class_ps3",
        "amine.aniline_flag",
        "amine.alpha_branch",
        "amine.formal_charge",
        "amine.is_aromatic_N",
        "amine.aryl_neighbors",
        "amine.alpha_heavy_neighbors",
        "amine.h_count_on_N",
    ]
    present = 1 if center_atoms else 0
    if not present:
        return [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    if not rdkit_available():
        # Present but RDKit unavailable: mark presence, rest defaults
        return [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
    if mol is None:
        return [0.0, -1.0, 0.0, 0.0, 0.0], fields
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    # Class
    cls = float(_amine_class(mol))
    # Aniline flag
    patt_aniline = Chem.MolFromSmarts("[$([NX3;H1][c]),$([NX3]([c])[H])]")
    aniline = 1.0 if mol.HasSubstructMatch(patt_aniline) else 0.0
    # Alpha branching: max heavy neighbor count of any N
    alpha = 0
    for ai in center_atoms:
        try:
            at = mol.GetAtomWithIdx(int(ai))
        except Exception:
            continue
        hn = sum(1 for nb in at.GetNeighbors() if nb.GetAtomicNum() > 1)
        alpha = max(alpha, hn)
    alpha_flag = 1.0 if alpha >= 3 else (0.5 if alpha == 2 else 0.0)
    alpha_heavy = float(alpha)
    # Formal charge on N (first center)
    charge = 0.0
    if center_atoms:
        try:
            charge = float(mol.GetAtomWithIdx(int(center_atoms[0])).GetFormalCharge())
        except Exception:
            charge = 0.0
    # Aromatic N and aryl neighbor count
    is_arom_N = 0.0
    aryl_neighbors = 0.0
    try:
        if center_atoms:
            n = mol.GetAtomWithIdx(int(center_atoms[0]))
            is_arom_N = 1.0 if n.GetIsAromatic() else 0.0
            count = 0
            for nb in n.GetNeighbors():
                if nb.GetIsAromatic() and nb.GetAtomicNum() == 6:
                    count += 1
            aryl_neighbors = float(count)
    except Exception:
        pass
    # Hydrogens on N (total)
    hcount = 0.0
    try:
        if center_atoms:
            hcount = float(mol.GetAtomWithIdx(int(center_atoms[0])).GetTotalNumHs())
    except Exception:
        hcount = 0.0
    return [1.0, cls, aniline, alpha_flag, charge, is_arom_N, aryl_neighbors, alpha_heavy, hcount], fields
