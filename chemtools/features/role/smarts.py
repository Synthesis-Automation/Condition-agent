from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, List, Sequence, Tuple

from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles
from chemtools.util.smarts_cache import compile_smarts
from chemtools.util.functional_groups import FUNCTIONAL_GROUP_SMARTS


_ROLE_TO_SPEC = {
    "amine_prim": ("amine_primary",),
    "amine_sec": ("amine_secondary",),
    "amine_tert": ("amine_tertiary",),
    "aniline": ("aniline",),
    "alcohol": ("alcohol", "phenol"),
    "aryl_F": ("aryl_fluoride",),
    "aryl_Cl": ("aryl_chloride",),
    "aryl_Br": ("aryl_bromide",),
    "aryl_I": ("aryl_iodide",),
}

_FALLBACK_SMARTS = {
    "amine_prim": "[NX3;H2][#6]",
    "amine_sec": "[NX3;H1]([#6])[#6]",
    "amine_tert": "[NX3]([#6])([#6])[#6]",
    "aniline": "[$([NX3;H1][c]),$([NX3]([c])[H])]",
    "alcohol": "[OX2H]",
    "aryl_F": "cF",
    "aryl_Cl": "cCl",
    "aryl_Br": "cBr",
    "aryl_I": "cI",
}


def _as_tuple(value: Sequence[str] | str | None) -> Tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    return tuple(value)


@lru_cache(maxsize=None)
def _patterns_for(key: str) -> Tuple[str, ...]:
    patterns: List[str] = []
    for spec_name in _ROLE_TO_SPEC.get(key, ()):
        spec_patterns = FUNCTIONAL_GROUP_SMARTS.get(spec_name)
        if spec_patterns:
            patterns.extend(_as_tuple(spec_patterns))
    if not patterns:
        patterns.extend(_as_tuple(_FALLBACK_SMARTS.get(key)))
    # Remove duplicates while preserving order
    seen = set()
    deduped: List[str] = []
    for smarts in patterns:
        if smarts in seen:
            continue
        seen.add(smarts)
        deduped.append(smarts)
    return tuple(deduped)


def _match_atom_indices(mol, key: str) -> List[int]:
    atoms: List[int] = []
    for smarts in _patterns_for(key):
        patt = compile_smarts(smarts, validate=False)
        if patt is None:
            continue
        try:
            matches = mol.GetSubstructMatches(patt)
        except Exception:
            continue
        for match in matches:
            if match:
                atoms.append(match[0])
    return atoms


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
            am_atoms.extend(_match_atom_indices(mol, key))
        
        al_atoms: List[int] = []
        al_atoms.extend(_match_atom_indices(mol, "alcohol"))
        
        ar_atoms: List[int] = []
        for key in ("aryl_I", "aryl_Br", "aryl_Cl", "aryl_F"):
            matches = _match_atom_indices(mol, key)
            if matches:
                ar_atoms.extend(matches)
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

