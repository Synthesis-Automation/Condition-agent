from __future__ import annotations

from typing import Any, Dict, List, Tuple, Optional

from .registry import REGISTRY
from . import smarts as _smarts
from .global_feats import compute_global
from .role_feats.amine import amine_features
from .role_feats.alcohol import alcohol_features
from .role_feats.aryl_halide import aryl_halide_features
from .fingerprints import centered_ecfp

try:
    import numpy as np  # type: ignore
except Exception:  # pragma: no cover
    np = None  # type: ignore


def _ensure_numpy_array(xs: List[float]):
    if np is not None:
        return np.asarray(xs, dtype=float)
    return xs  # type: ignore


def _role_order() -> List[str]:
    return ["amine", "alcohol", "aryl_halide"]


def _role_present_mask(found: Dict[str, Any]) -> Dict[str, int]:
    return {r: (1 if (found.get(r) and found[r].get("present")) else 0) for r in _role_order()}


def featurize_mol(mol_or_smiles: Any, roles: Optional[List[str]] = None) -> Dict[str, Any]:
    """Role-aware featurization for a single molecule.

    Returns {vector, fields, masks, meta}. Deterministic and safe without RDKit.
    """
    if roles is None:
        roles = _role_order()

    # Find reactive centers per role
    centers = _smarts.find_centers(mol_or_smiles)

    # Global descriptor block
    gvals, gfields = compute_global(mol_or_smiles)

    # Role-specific blocks
    found: Dict[str, Dict[str, Any]] = {}
    role_blocks: List[Tuple[List[float], List[str]]] = []

    for r in roles:
        if r == "amine":
            vals, fields = amine_features(mol_or_smiles, centers.get("amine", {}).get("atoms", []))
        elif r == "alcohol":
            vals, fields = alcohol_features(mol_or_smiles, centers.get("alcohol", {}).get("atoms", []))
        elif r == "aryl_halide":
            vals, fields = aryl_halide_features(mol_or_smiles, centers.get("aryl_halide", {}).get("atoms", []))
        else:
            vals, fields = ([], [])
        # track present flags for masks
        present = 1 if (centers.get(r, {}).get("atoms")) else 0
        found[r] = {"present": present}
        role_blocks.append((vals, fields))

    # Centered fingerprints per role
    fp_blocks: List[Tuple[List[float], List[str]]] = []
    for r in roles:
        cfg = REGISTRY.get("fingerprints", {}).get(r, {"bits": 512, "radius": 2})
        bits = int(cfg.get("bits", 512))
        radius = int(cfg.get("radius", 2))
        atoms = centers.get(r, {}).get("atoms", [])
        fp = centered_ecfp(mol_or_smiles, atoms, bits=bits, radius=radius)
        fields = [f"fp.{r}.{i}" for i in range(bits)]
        fp_blocks.append((fp, fields))

    # Assemble vector in fixed order: global -> roles -> fingerprints
    vec: List[float] = []
    fields: List[str] = []
    vec.extend(gvals); fields.extend(gfields)
    for vals, flds in role_blocks:
        vec.extend(vals); fields.extend(flds)
    for vals, flds in fp_blocks:
        vec.extend(vals); fields.extend(flds)

    return {
        "vector": _ensure_numpy_array(vec),
        "fields": fields,
        "masks": _role_present_mask(found),
        "meta": {"centers": centers, "roles": roles},
    }


def featurize_reaction(rxn_smiles: str) -> Dict[str, Any]:
    """Featurize all reactants in a reaction SMILES string.

    Returns {reactants: [ {vector, fields, masks, meta, smiles} ]}
    """
    parts = (rxn_smiles or "").split(">")
    if len(parts) == 2 and ">>" in (rxn_smiles or ""):
        reactants = parts[0]
    elif len(parts) >= 1:
        reactants = parts[0]
    else:
        reactants = rxn_smiles
    mols = [s for s in (reactants or "").split(".") if s]
    out = []
    for smi in mols:
        out.append({**featurize_mol(smi), "smiles": smi})
    return {"reactants": out}

