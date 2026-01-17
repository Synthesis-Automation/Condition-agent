"""
General molecular featurizer for electrophile/nucleophile pairs.

Extracts structural and chemical features across cross-coupling and related
reaction spaces (C-N, C-O, C-S, etc.). Features include leaving-group type,
electrophile classification, nucleophile basicity, steric properties, and
optional role-aware/descriptor attachments.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional, Sequence, Tuple
from functools import lru_cache
import os
import re

from ..util.rdkit_helpers import rdkit_available, parse_smiles
from ..util.smarts_cache import compile_smarts

try:  # Optional role-aware vectors
    from chemtools.features.role import featurize_mol as _role_feat  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:  # pragma: no cover
    _role_feat = None  # type: ignore
    _HAS_ROLE_FEATS = False


_HALOGEN_GROUPS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("I", ("aryl_iodide", "vinyl_iodide", "alkyl_iodide", "iodide")),
    ("Br", ("aryl_bromide", "vinyl_bromide", "alkyl_bromide")),
    ("Cl", ("aryl_chloride", "vinyl_chloride", "alkyl_chloride")),
)

_ARYL_GROUPS = ("aryl_halide", "aryl_chloride", "aryl_bromide", "aryl_iodide")
_VINYL_GROUPS = ("vinyl_halide", "vinyl_chloride", "vinyl_bromide", "vinyl_iodide")
_ALKYL_GROUPS = ("alkyl_halide", "alkyl_chloride", "alkyl_bromide", "alkyl_iodide")

_EWG_GROUPS = ("nitro", "nitrile", "carbonyl", "trifluoromethyl")

_NUCLEOPHILE_PRIORITY: Tuple[Tuple[Tuple[str, ...], str, str], ...] = (
    (("indole",), "indole", "deactivated"),
    (("aniline",), "aniline", "aromatic_primary"),
    (("phenol",), "phenol", "deactivated"),
    (("amide", "amide_primary", "amide_secondary"), "amide_deactivated", "deactivated"),
    (("amine_secondary",), "amine_secondary", "secondary"),
    (("amine_primary",), "amine_primary", "aliphatic_primary"),
    (("amine_tertiary",), "amine_tertiary", "tertiary"),
)


def _get_env_bool(name: str, default: bool = False) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v.strip().lower() in {"1", "true", "yes", "on"}


def _maybe_structural(smiles: str) -> Optional[Dict[str, Any]]:
    try:
        from .structural import featurize_molecule
    except Exception:
        return None
    try:
        return featurize_molecule(smiles)
    except Exception:
        return None


def _guess_lg_text(s: str) -> str:
    """Heuristically guess leaving group from SMILES text."""
    t = (s or "").lower()
    if "os(=o)(=o)c(f)(f)f" in t or "otf" in t:
        return "OTf"
    if "br" in t:
        return "Br"
    if "cl" in t:
        return "Cl"
    if "i" in t:
        return "I"
    return "UNK"


def _iter_matches(mol, name: str, fallback: Optional[Sequence[str]] = None):
    patterns = tuple(fallback or ())
    for smarts in patterns:
        patt = compile_smarts(smarts, validate=False)
        if patt is None:
            continue
        try:
            matches = mol.GetSubstructMatches(patt)
        except Exception:
            continue
        for match in matches:
            yield match


def _find_halide_ipso(mol, halogen: str):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != halogen:
            continue
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 6:
                return nb
    return None


def _find_triflate_ipso(mol):
    for match in _iter_matches(mol, "triflate", ("OS(=O)(=O)C(F)(F)F",)):
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "O":
                continue
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 6 and nb.GetIdx() not in match:
                    return nb
    return None


def _infer_leaving_group(fg_hits: Dict[str, bool], smiles: str) -> Tuple[str, Optional[str]]:
    text_guess = _guess_lg_text(smiles)
    if fg_hits.get("triflate"):
        return "OTf", "aryl"

    for symbol, names in _HALOGEN_GROUPS:
        for name in names:
            if fg_hits.get(name):
                if "aryl" in name:
                    return symbol, "aryl"
                if "vinyl" in name:
                    return symbol, "vinyl"
                if "alkyl" in name:
                    return symbol, "alkyl"
                return symbol, None

    if fg_hits.get("aryl_halide"):
        return text_guess, "aryl"
    if fg_hits.get("vinyl_halide"):
        return text_guess, "vinyl"
    if fg_hits.get("alkyl_halide"):
        return text_guess, "alkyl"

    return text_guess, None


def _has_ewg(fg_hits: Dict[str, bool], smiles: str) -> bool:
    if any(fg_hits.get(name) for name in _EWG_GROUPS):
        return True
    t = (smiles or "").lower()
    return any(token in t for token in ("[n+](=o)[o-]", "c#n", "c(f)(f)f", "s(=o)(=o)"))


def _has_boron(smiles: str) -> bool:
    if not smiles:
        return False
    return re.search(r"B(?!r)", smiles) is not None


def _infer_class_from_hits(fg_hits: Dict[str, bool], fallback: Optional[str]) -> str:
    if fallback:
        return fallback
    if any(fg_hits.get(name) for name in _ARYL_GROUPS):
        return "aryl"
    if any(fg_hits.get(name) for name in _VINYL_GROUPS):
        return "vinyl"
    if any(fg_hits.get(name) for name in _ALKYL_GROUPS):
        return "alkyl"
    return "aryl"


def _detect_electrophile(
    mol,
    fg_hits: Optional[Dict[str, bool]] = None,
    smiles_text: str = ""
) -> Dict[str, Any]:
    """Extract electrophile features (LG, class, ortho substitution, EWG, heteroaryl)."""
    fg_hits = fg_hits or {}
    lg, class_hint = _infer_leaving_group(fg_hits, smiles_text)
    elec_class = _infer_class_from_hits(fg_hits, class_hint)
    para_ewg = _has_ewg(fg_hits, smiles_text)

    result = {
        "LG": lg,
        "elec_class": elec_class,
        "ortho_count": "0",
        "para_EWG": para_ewg if (lg != "UNK" and elec_class in ("aryl", "vinyl", "alkenyl")) else False,
        "heteroaryl": False,
    }

    if mol is None or not rdkit_available():
        return result

    ipso_atom = None
    if lg == "OTf":
        ipso_atom = _find_triflate_ipso(mol)
    elif lg in {"Cl", "Br", "I"}:
        ipso_atom = _find_halide_ipso(mol, lg)

    if ipso_atom is None:
        for symbol, _ in _HALOGEN_GROUPS:
            candidate = _find_halide_ipso(mol, symbol)
            if candidate is not None:
                if lg == "UNK":
                    lg = symbol
                ipso_atom = candidate
                break

    ortho_count = 0
    heteroaryl = False
    if ipso_atom is not None and ipso_atom.GetIsAromatic():
        elec_class = "aryl"
        ri = ipso_atom.GetOwningMol().GetRingInfo()
        atom_idx = ipso_atom.GetIdx()
        rings = [r for r in ri.AtomRings() if atom_idx in r and len(r) == 6]
        if rings:
            ring = rings[0]
            pos = ring.index(atom_idx)
            ortho_atoms = [
                mol.GetAtomWithIdx(ring[(pos - 1) % 6]),
                mol.GetAtomWithIdx(ring[(pos + 1) % 6]),
            ]
            heteroaryl = any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6, 1) for i in ring)

            def is_substituted(ar_atom):
                count = 0
                for nb in ar_atom.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetAtomicNum() > 1:
                        count += 1
                return count > 0

            ortho_count = sum(1 for atom in ortho_atoms if is_substituted(atom))
    else:
        if any(fg_hits.get(name) for name in _VINYL_GROUPS):
            elec_class = "vinyl"
        elif any(fg_hits.get(name) for name in _ALKYL_GROUPS):
            elec_class = "alkyl"

    para_ewg_flag = para_ewg if (lg != "UNK" and elec_class in ("aryl", "vinyl", "alkenyl")) else False

    result.update({
        "LG": lg,
        "elec_class": elec_class,
        "ortho_count": "2+" if ortho_count >= 2 else ("1" if ortho_count == 1 else "0"),
        "para_EWG": para_ewg_flag,
        "heteroaryl": heteroaryl,
    })
    return result


def _nuc_class_text(s: str) -> str:
    """Heuristically classify nucleophile from SMILES text."""
    t = (s or "").lower()
    if "indole" in t:
        return "indole"
    if re.search(r"c[^)]*n", t):
        return "aniline"
    # Fallback aromatic detection when RDKit is unavailable: detect common ring
    # patterns like ``c1ccccc1`` combined with an amine token.
    if "n" in t and any(marker in t for marker in ("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8")):
        return "aniline"
    if "n(" in t:
        return "amine_secondary"
    if "n" in t:
        return "amine_primary"
    if "o" in t and "n" not in t:
        return "phenol"
    return "amine"


def _nucleophile_features(
    mol,
    text: str,
    fg_hits: Optional[Dict[str, bool]] = None
) -> Dict[str, Any]:
    """Extract nucleophile features (class, basicity, sterics)."""
    fg_hits = fg_hits or {}
    nuc_class: Optional[str] = None
    n_basicity = "unknown"
    steric_alpha = "low"

    for group_names, cls, basicity in _NUCLEOPHILE_PRIORITY:
        if any(fg_hits.get(name) for name in group_names):
            nuc_class = cls
            n_basicity = basicity
            break

    if nuc_class is None:
        nuc_class = _nuc_class_text(text)
        lowered = (text or "").lower()
        if nuc_class == "aniline":
            n_basicity = "aromatic_primary"
        elif nuc_class == "amine_primary":
            n_basicity = "aliphatic_primary"
        elif nuc_class == "amine_secondary":
            n_basicity = "secondary"
        elif nuc_class == "phenol":
            n_basicity = "deactivated"
        elif "n" in lowered:
            nuc_class = "amine_primary"
            n_basicity = "aliphatic_primary"
        elif "o" in lowered:
            nuc_class = "phenol"
            n_basicity = "deactivated"

    if mol is None or not rdkit_available():
        return {"nuc_class": nuc_class, "n_basicity": n_basicity, "steric_alpha": steric_alpha}

    # Sterics at alpha: approximate via heavy neighbor count of the nucleophilic heteroatom (N/O)
    steric_level = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ("N", "O"):
            # Count heavy atom neighbors (exclude H)
            hn = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
            steric_level = max(steric_level, hn)
    steric_alpha = "high" if steric_level >= 3 else ("med" if steric_level == 2 else "low")

    return {"nuc_class": nuc_class, "n_basicity": n_basicity, "steric_alpha": steric_alpha}


@lru_cache(maxsize=1024)
def _featurize_cached(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    """Core featurization logic with caching."""
    # Parse molecules if RDKit available
    emol = parse_smiles(electrophile)
    nmol = parse_smiles(nucleophile)

    elec = _detect_electrophile(emol, {}, electrophile)
    nuc = _nucleophile_features(nmol, nucleophile, {})

    # Compose bin key (coarse)
    lg = elec.get("LG", "UNK")
    nuc_class = nuc.get("nuc_class", "unknown")
    bin_key = f"LG:{lg}|NUC:{nuc_class}"

    return {
        **elec,
        **nuc,
        "bin": bin_key,
    }


def featurize_flat(
    electrophile: str,
    nucleophile: str,
    *,
    include_calculable: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return the legacy flat feature dictionary for a reactant pair."""
    base = _featurize_cached(electrophile, nucleophile)
    out = dict(base)

    # Attach calculable features when explicitly enabled via env or parameter.
    if include_calculable is None:
        include_calculable = _get_env_bool("CHEMTOOLS_INCLUDE_CALCULABLE_FEATURES", False)
    
    if include_calculable:
        try:
            from . import calculable
            calculable_data = {}
            if electrophile:
                calculable_data["electrophile"] = calculable.detect_all_features(electrophile)
            if nucleophile:
                calculable_data["nucleophile"] = calculable.detect_all_features(nucleophile)
            if calculable_data:
                out["calculable"] = calculable_data
        except Exception:
            pass  # Graceful degradation if calculable module unavailable

    # Attach role-aware only when explicitly enabled via env (default off for speed).
    attach_flag = (os.environ.get("CHEMTOOLS_ATTACH_ROLE_AWARE", "").strip().lower() in {"1", "true", "yes", "on"})
    if attach_flag and _HAS_ROLE_FEATS and _role_feat is not None:
        try:
            elec_ra = _role_feat(electrophile, roles=["aryl_halide"])  # type: ignore[arg-type]
            nuc_ra = _role_feat(nucleophile, roles=["amine", "alcohol"])  # type: ignore[arg-type]
            out["role_aware"] = {
                "electrophile": elec_ra,
                "nucleophile": nuc_ra,
            }
        except Exception:
            pass

    # Standardized substrate-tag tokens (expandable).
    ortho = out.get("ortho_count")
    if isinstance(ortho, str):
        out["ortho_0"] = ortho == "0"
        out["ortho_1"] = ortho == "1"
        out["ortho_2plus"] = ortho == "2+"
        out["ortho_hindered"] = ortho == "2+"

    ewg_on_electrophile = _has_ewg({}, electrophile)
    out["electron_poor_aryl"] = bool(
        out.get("elec_class") == "aryl"
        and (
            bool(out.get("heteroaryl"))
            or bool(out.get("para_EWG"))
            or ewg_on_electrophile
        )
    )

    boron_on_nucleophile = _has_boron(nucleophile)
    ewg_on_nucleophile = _has_ewg({}, nucleophile)
    out["electron_poor_boronate"] = bool(boron_on_nucleophile and ewg_on_nucleophile)
    out["electron_poor_boronate_present"] = out["electron_poor_boronate"]

    return out


def featurize_pair(
    electrophile: str,
    nucleophile: str,
    *,
    include_calculable: Optional[bool] = None,
    include_structural: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return canonical, structured features for an electrophile/nucleophile pair."""
    errors: List[str] = []
    if not (electrophile or nucleophile):
        errors.append("empty_inputs")

    flat = featurize_flat(
        electrophile,
        nucleophile,
        include_calculable=include_calculable,
    )

    if include_structural is None:
        include_structural = _get_env_bool("CHEMTOOLS_INCLUDE_STRUCTURAL_FEATURES", False)

    resolved_include_calculable = (
        include_calculable if include_calculable is not None else ("calculable" in flat)
    )

    elec_structural = _maybe_structural(electrophile) if include_structural and electrophile else None
    nuc_structural = _maybe_structural(nucleophile) if include_structural and nucleophile else None

    elec_features = {
        "LG": flat.get("LG"),
        "elec_class": flat.get("elec_class"),
        "ortho_count": flat.get("ortho_count"),
        "para_EWG": flat.get("para_EWG"),
        "heteroaryl": flat.get("heteroaryl"),
    }
    nuc_features = {
        "nuc_class": flat.get("nuc_class"),
        "n_basicity": flat.get("n_basicity"),
        "steric_alpha": flat.get("steric_alpha"),
    }

    calculable_block = flat.get("calculable")
    role_aware_block = flat.get("role_aware")

    elec_payload: Dict[str, Any] = {
        "smiles": electrophile,
        "features": elec_features,
    }
    nuc_payload: Dict[str, Any] = {
        "smiles": nucleophile,
        "features": nuc_features,
    }

    if isinstance(calculable_block, dict):
        elec_calc = calculable_block.get("electrophile")
        nuc_calc = calculable_block.get("nucleophile")
        if elec_calc:
            elec_payload["calculable"] = elec_calc
        if nuc_calc:
            nuc_payload["calculable"] = nuc_calc

    if isinstance(role_aware_block, dict):
        elec_ra = role_aware_block.get("electrophile")
        nuc_ra = role_aware_block.get("nucleophile")
        if elec_ra:
            elec_payload["role_aware"] = elec_ra
        if nuc_ra:
            nuc_payload["role_aware"] = nuc_ra

    if elec_structural:
        elec_payload["structural"] = elec_structural
    if nuc_structural:
        nuc_payload["structural"] = nuc_structural

    pair_features = {
        "LG": flat.get("LG"),
        "elec_class": flat.get("elec_class"),
        "ortho_count": flat.get("ortho_count"),
        "para_EWG": flat.get("para_EWG"),
        "heteroaryl": flat.get("heteroaryl"),
        "nuc_class": flat.get("nuc_class"),
        "n_basicity": flat.get("n_basicity"),
        "steric_alpha": flat.get("steric_alpha"),
    }

    return {
        "schema_version": "v2",
        "electrophile": elec_payload,
        "nucleophile": nuc_payload,
        "pair": {
            "bin": flat.get("bin"),
            "features": pair_features,
        },
        "flat": flat,
        "meta": {
            "rdkit_available": rdkit_available(),
            "errors": errors,
            "options": {
                "include_calculable": resolved_include_calculable,
                "include_structural": include_structural,
            },
        },
    }


def featurize(
    electrophile: str,
    nucleophile: str,
    *,
    include_calculable: Optional[bool] = None,
    include_structural: Optional[bool] = None,
) -> Dict[str, Any]:
    """Alias for featurize_pair."""
    return featurize_pair(
        electrophile,
        nucleophile,
        include_calculable=include_calculable,
        include_structural=include_structural,
    )


__all__ = ["featurize", "featurize_pair", "featurize_flat"]
