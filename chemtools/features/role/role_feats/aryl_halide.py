from __future__ import annotations

from typing import Any, List, Tuple

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles


def _halide_code(mol) -> int:
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return 0
    for sym, code in (("I", 4), ("Br", 3), ("Cl", 2), ("F", 1)):
        patt = Chem.MolFromSmarts(f"c[{sym}]")
        if mol.HasSubstructMatch(patt):
            return code
    return 0


def aryl_halide_features(mol_or_smiles: Any, center_atoms: List[int]) -> Tuple[List[float], List[str]]:
    fields = [
        "aryl_halide.present",
        "aryl_halide.halide",
        "aryl_halide.ortho_block",
        "aryl_halide.ipso_degree",
        "aryl_halide.para_EWG",
        "aryl_halide.heteroaryl",
        "aryl_halide.ring_aromatic_count",
        "aryl_halide.triflate_flag",
    ]
    present = 1 if center_atoms else 0
    if not present:
        return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    if not rdkit_available():
        return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
    if mol is None:
        return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], fields
    try:
        hal = float(_halide_code(mol))
    except Exception:
        hal = 0.0
    # Ortho block approx via two ortho positions substituted
    ortho = 0.0
    ipso_deg = 0.0
    heteroaryl = 0.0
    ring_aromatic_count = 0.0
    para_ewg = 0.0
    triflate = 0.0
    try:
        for ai in center_atoms:
            at = mol.GetAtomWithIdx(int(ai))
            ipso_deg = float(at.GetTotalDegree())
            # count ortho substituents on aromatic ring
            ring_info = at.GetOwningMol().GetRingInfo()
            rings = [r for r in ring_info.AtomRings() if int(ai) in r and len(r) == 6]
            if rings:
                ring = rings[0]
                # heteroaryl if any non-C/H atom in ring
                heteroaryl = 1.0 if any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6) for i in ring) else 0.0
                pos = ring.index(int(ai))
                o1 = ring[(pos - 1) % 6]
                o2 = ring[(pos + 1) % 6]
                def subs(idx: int) -> int:
                    a = mol.GetAtomWithIdx(idx)
                    return sum(1 for nb in a.GetNeighbors() if nb.GetIdx() not in ring and nb.GetAtomicNum() > 1)
                c = 0
                if subs(o1) > 0:
                    c += 1
                if subs(o2) > 0:
                    c += 1
                ortho = float(c)
                break
        # count aromatic rings in molecule
        ring_info = mol.GetRingInfo()
        for r in ring_info.AtomRings():
            try:
                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
                    ring_aromatic_count += 1.0
            except Exception:
                continue
    except Exception:
        pass
    # EWG and triflate flags (anywhere)
    try:
        from rdkit import Chem  # type: ignore
        ewg_patts = [
            Chem.MolFromSmarts("[N+](=O)[O-]"),
            Chem.MolFromSmarts("C#N"),
            Chem.MolFromSmarts("C(F)(F)F"),
            Chem.MolFromSmarts("C(=O)[!O]"),
        ]
        if any(mol.HasSubstructMatch(p) for p in ewg_patts if p is not None):
            para_ewg = 1.0
        patt_triflate = Chem.MolFromSmarts("OS(=O)(=O)C(F)(F)F")
        if patt_triflate is not None and mol.HasSubstructMatch(patt_triflate):
            triflate = 1.0
    except Exception:
        pass
    return [1.0, hal, ortho, ipso_deg, para_ewg, heteroaryl, ring_aromatic_count, triflate], fields
