from __future__ import annotations

from typing import Any, List

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles

try:
    import numpy as np  # type: ignore
except Exception:  # pragma: no cover
    np = None  # type: ignore


def centered_ecfp(mol_or_smiles: Any, center_atoms: List[int], bits: int = 512, radius: int = 2) -> List[float]:
    if not center_atoms or not rdkit_available():
        return [0.0] * int(bits)
    mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
    if mol is None:
        return [0.0] * int(bits)
    try:
        from rdkit.Chem import AllChem  # type: ignore
    except Exception:
        return [0.0] * int(bits)
    try:
        bv = AllChem.GetMorganFingerprintAsBitVect(mol, int(radius), nBits=int(bits), useFeatures=False, fromAtoms=[int(i) for i in center_atoms])
        arr = np.zeros((int(bits),), dtype=int) if np is not None else [0] * int(bits)
        if np is not None:
            from rdkit import DataStructs  # type: ignore
            DataStructs.ConvertToNumpyArray(bv, arr)  # type: ignore
            return [float(x) for x in arr.tolist()]
        # Fallback: try to iterate bits
        return [1.0 if bv.GetBit(i) else 0.0 for i in range(int(bits))]
    except Exception:
        return [0.0] * int(bits)

