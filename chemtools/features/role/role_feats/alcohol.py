from __future__ import annotations

from typing import Any, List, Tuple

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles


def alcohol_features(mol_or_smiles: Any, center_atoms: List[int]) -> Tuple[List[float], List[str]]:
    fields = [
        "alcohol.present",
        "alcohol.class_ps3",
        "alcohol.benzylic",
        "alcohol.allylic",
        "alcohol.phenol_flag",
        "alcohol.alpha_degree",
        "alcohol.alpha_heavy_neighbors",
        "alcohol.aryl_neighbors",
    ]
    present = 1 if center_atoms else 0
    if not present:
        return [0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0], fields
    if not rdkit_available():
        return [1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0], fields
    mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
    if mol is None:
        return [0.0, -1.0, 0.0, 0.0], fields
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return [1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0], fields
    # Class by alpha carbon degree for first O center
    cls = -1.0
    benzylic = 0.0
    allylic = 0.0
    phenol = 0.0
    alpha_deg = -1.0
    alpha_heavy = 0.0
    aryl_neighbors = 0.0
    for ai in center_atoms:
        a = mol.GetAtomWithIdx(int(ai))
        # find carbon neighbor
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 6:
                deg = nb.GetTotalDegree()
                alpha_deg = float(deg)
                cls = 1.0 if deg <= 2 else (2.0 if deg == 3 else 3.0)
                # benzylic if carbon is attached to aromatic carbon
                benzylic = 1.0 if any(n2.GetIsAromatic() for n2 in nb.GetNeighbors()) else 0.0
                # aryl neighbor count (around alpha carbon)
                try:
                    aryl_neighbors = float(sum(1 for n2 in nb.GetNeighbors() if n2.GetIsAromatic() and n2.GetAtomicNum() == 6))
                except Exception:
                    aryl_neighbors = 0.0
                # allylic if carbon has a double-bonded neighbor carbon
                allylic = 0.0
                for b in nb.GetBonds():
                    if b.GetBondTypeAsDouble() == 2.0:
                        j = b.GetOtherAtomIdx(nb.GetIdx())
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 6:
                            allylic = 1.0
                            break
                # alpha heavy neighbors (excluding O)
                alpha_heavy = float(sum(1 for n2 in nb.GetNeighbors() if n2.GetAtomicNum() > 1 and n2.GetSymbol() != 'O'))
                break
        break
    # Phenol: O attached to aromatic carbon
    try:
        for ai in center_atoms:
            a = mol.GetAtomWithIdx(int(ai))
            if any(nb.GetIsAromatic() and nb.GetAtomicNum() == 6 for nb in a.GetNeighbors()):
                phenol = 1.0
                break
    except Exception:
        phenol = 0.0
    return [1.0, cls, benzylic, allylic, phenol, alpha_deg, alpha_heavy, aryl_neighbors], fields
