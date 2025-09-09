from __future__ import annotations

from typing import Any, List, Tuple

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles


def compute_global(mol_or_smiles: Any) -> Tuple[List[float], List[str]]:
    fields = [
        "global.MW",
        "global.logP",
        "global.TPSA",
        "global.rotB",
        "global.HBD",
        "global.HBA",
        "global.aromatic_rings",
        "global.heteroatoms",
    ]
    if not rdkit_available():
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    try:
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors  # type: ignore
    except Exception:
        return [0.0] * len(fields), fields
    mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
    if mol is None:
        return [0.0] * len(fields), fields
    try:
        mw = float(Descriptors.MolWt(mol))
    except Exception:
        mw = 0.0
    try:
        from rdkit.Chem import Crippen  # type: ignore
        logp = float(Crippen.MolLogP(mol))
    except Exception:
        logp = 0.0
    try:
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    except Exception:
        tpsa = 0.0
    try:
        rotb = float(Lipinski.NumRotatableBonds(mol))
    except Exception:
        rotb = 0.0
    try:
        hbd = float(Lipinski.NumHDonors(mol))
    except Exception:
        hbd = 0.0
    try:
        hba = float(Lipinski.NumHAcceptors(mol))
    except Exception:
        hba = 0.0
    try:
        ring_info = mol.GetRingInfo()
        aromatic_rings = 0
        for ring in ring_info.AtomRings():
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                aromatic_rings += 1
    except Exception:
        aromatic_rings = 0
    try:
        hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    except Exception:
        hetero = 0
    vals: List[float] = [mw, logp, tpsa, rotb, hbd, hba, float(aromatic_rings), float(hetero)]
    return vals, fields

