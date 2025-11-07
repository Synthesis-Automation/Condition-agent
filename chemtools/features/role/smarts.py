from __future__ import annotations

from typing import Any, Dict, List

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles
from chemtools.util.smarts_cache import compile_smarts


# SMARTS patterns for role-based features
_SMARTS_PATTERNS = {
    # Amines (primary/secondary/tertiary)
    "amine_prim": "[NX3;H2][#6]",
    "amine_sec": "[NX3;H1]([#6])[#6]",
    "amine_tert": "[NX3]([#6])([#6])[#6]",
    # Aniline-like
    "aniline": "[$([NX3;H1][c]),$([NX3]([c])[H])]",
    # Alcohol OH
    "alcohol": "[OX2H]",
    # Aryl halide: aromatic carbon attached to halogen
    "aryl_F": "cF",
    "aryl_Cl": "cCl",
    "aryl_Br": "cBr",
    "aryl_I": "cI",
}


def find_centers(mol_or_smiles: Any) -> Dict[str, Dict[str, List[int]]]:
    """Locate reactive centers per role. Returns atom indices per role.

    Fallback to text heuristics if RDKit is unavailable.
    """
    if rdkit_available():
        mol = mol_or_smiles if getattr(mol_or_smiles, "GetAtoms", None) else parse_smiles(str(mol_or_smiles))
        if mol is None:
            return {"amine": {"atoms": []}, "alcohol": {"atoms": []}, "aryl_halide": {"atoms": []}}
        
        am_atoms: List[int] = []
        for key in ("amine_tert", "amine_sec", "amine_prim"):
            smarts = _SMARTS_PATTERNS.get(key)
            if smarts:
                patt = compile_smarts(smarts, validate=False)
                if patt is not None:
                    for match in mol.GetSubstructMatches(patt):
                        am_atoms.append(match[0])
        
        al_atoms: List[int] = []
        alcohol_smarts = _SMARTS_PATTERNS.get("alcohol")
        if alcohol_smarts:
            patt_al = compile_smarts(alcohol_smarts, validate=False)
            if patt_al is not None:
                for match in mol.GetSubstructMatches(patt_al):
                    al_atoms.append(match[0])
        
        ar_atoms: List[int] = []
        for key in ("aryl_I", "aryl_Br", "aryl_Cl", "aryl_F"):
            smarts = _SMARTS_PATTERNS.get(key)
            if smarts:
                patt = compile_smarts(smarts, validate=False)
                if patt is not None and mol.HasSubstructMatch(patt):
                    for m in mol.GetSubstructMatches(patt):
                        # ipso carbon index is first in patterns above
                        ar_atoms.append(m[0])
                    break
        
        return {
            "amine": {"atoms": sorted(set(am_atoms))},
            "alcohol": {"atoms": sorted(set(al_atoms))},
            "aryl_halide": {"atoms": sorted(set(ar_atoms))},
        }
    
    # Fallback heuristics when RDKit is not available
    t = str(mol_or_smiles or "").lower()
    am = bool("n" in t)
    al = bool("o" in t and "n" not in t)
    ar = any(x in t for x in ["cl", "br", " i", "os(=o)(=o)c(f)(f)f", "otf"]) and ("c1" in t or "c[" in t)
    return {
        "amine": {"atoms": [0] if am else []},
        "alcohol": {"atoms": [0] if al else []},
        "aryl_halide": {"atoms": [0] if ar else []},
    }

