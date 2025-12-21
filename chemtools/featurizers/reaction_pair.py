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
from ..util.functional_groups import (
    detect_all as detect_functional_groups,
    FUNCTIONAL_GROUP_SMARTS,
)
from ..util.smarts_cache import compile_smarts

try:  # Optional role-aware vectors
    from chemtools.features.role import featurize_mol as _role_feat  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:  # pragma: no cover
    _role_feat = None  # type: ignore
    _HAS_ROLE_FEATS = False


_MOLPIPELINE_SENTINEL = object()
_MOLPIPELINE_HELPERS: Any = _MOLPIPELINE_SENTINEL  # type: ignore[misc]
_DEFAULT_MOLPIPELINE_DESCRIPTOR_LIST = [
    "HeavyAtomMolWt",
    "TPSA",
    "MolLogP",
    "MolMR",
    "NumHAcceptors",
    "NumHDonors",
]

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


def _env_int(name: str, default: int) -> int:
    try:
        value = int(os.environ.get(name, str(default)))
    except Exception:
        return default
    return value if value > 0 else default


def _maybe_import_molpipeline():
    global _MOLPIPELINE_HELPERS
    if _MOLPIPELINE_HELPERS is None:
        return None
    if _MOLPIPELINE_HELPERS is _MOLPIPELINE_SENTINEL:
        try:
            from . import molpipeline as _mp_helpers  # type: ignore
        except Exception:
            _MOLPIPELINE_HELPERS = None
            return None
        _MOLPIPELINE_HELPERS = _mp_helpers
    return _MOLPIPELINE_HELPERS


def _to_float_list(arr: Any) -> Optional[List[float]]:
    if arr is None:
        return None
    try:
        import numpy as np

        return np.asarray(arr, dtype=float).ravel().tolist()
    except Exception:
        try:
            return [float(x) for x in arr]  # type: ignore[arg-type]
        except Exception:
            return None


def _molpipeline_vectors(smiles: str) -> Optional[Dict[str, Any]]:
    helpers = _maybe_import_molpipeline()
    if helpers is None:
        return None
    settings = {
        "morgan_bits": _env_int("CHEMTOOLS_MOLPIPE_MORGAN_BITS", 1024),
        "morgan_radius": _env_int("CHEMTOOLS_MOLPIPE_MORGAN_RADIUS", 2),
        "physchem_descriptors": _DEFAULT_MOLPIPELINE_DESCRIPTOR_LIST,
    }
    descriptor_env = os.environ.get("CHEMTOOLS_MOLPIPE_PHYS_DESC")
    if descriptor_env:
        parsed = [part.strip() for part in descriptor_env.split(",") if part.strip()]
        if parsed:
            settings["physchem_descriptors"] = parsed

    try:
        morgan = helpers.morgan_fingerprint(
            smiles,
            n_bits=settings["morgan_bits"],
            radius=settings["morgan_radius"],
            return_sparse=False,
        )
        morgan_list = _to_float_list(morgan)
        physchem = helpers.physchem_features(
            smiles,
            descriptor_list=settings["physchem_descriptors"],
        )
        physchem_list = _to_float_list(physchem)
    except Exception:
        return None

    result: Dict[str, Any] = {}
    if morgan_list is not None:
        result["morgan_fp"] = morgan_list
    if physchem_list is not None:
        descriptor_names = settings.get("physchem_descriptors", [])
        result["physchem"] = physchem_list
        if descriptor_names and len(descriptor_names) == len(physchem_list):
            result["physchem_map"] = dict(zip(descriptor_names, physchem_list))
    if not result:
        return None
    result["_settings"] = settings
    return result


def _present_groups(groups: Dict[str, bool]) -> List[str]:
    return sorted(name for name, present in groups.items() if present)


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


def _fg_smarts(name: str, fallback: Optional[Sequence[str]] = None) -> Tuple[str, ...]:
    patterns = FUNCTIONAL_GROUP_SMARTS.get(name)
    if patterns:
        return tuple(patterns)
    if fallback:
        return tuple(fallback)
    return ()


def _iter_matches(mol, name: str, fallback: Optional[Sequence[str]] = None):
    for smarts in _fg_smarts(name, fallback):
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

    elec_groups = detect_functional_groups(electrophile)
    nuc_groups = detect_functional_groups(nucleophile)

    elec = _detect_electrophile(emol, elec_groups, electrophile)
    nuc = _nucleophile_features(nmol, nucleophile, nuc_groups)

    # Compose bin key (coarse)
    lg = elec.get("LG", "UNK")
    nuc_class = nuc.get("nuc_class", "unknown")
    bin_key = f"LG:{lg}|NUC:{nuc_class}"

    return {
        **elec,
        **nuc,
        "bin": bin_key,
        "_elec_groups": elec_groups,
        "_nuc_groups": nuc_groups,
    }


def _enrich_with_functional_groups(
    features: Dict[str, Any],
    electrophile: str,
    nucleophile: str,
    elec_groups: Optional[Dict[str, bool]] = None,
    nuc_groups: Optional[Dict[str, bool]] = None,
) -> Dict[str, Any]:
    """
    Enrich feature dict with comprehensive functional group detection.
    
    Adds functional-group detection tokens (already normalized as *_present)
    plus special aliases for SNAr and reductive amination reactions.
    
    Args:
        features: Base feature dictionary to enrich
        electrophile: SMILES string of electrophile
        nucleophile: SMILES string of nucleophile
        
    Returns:
        Enriched feature dictionary with functional group tokens
    """
    enriched = dict(features)
    
    # Detect functional groups in both reactants
    if elec_groups is None:
        elec_groups = detect_functional_groups(electrophile)
    if nuc_groups is None:
        nuc_groups = detect_functional_groups(nucleophile)

    for group, present in elec_groups.items():
        if present:
            enriched[group] = True
    
    for group, present in nuc_groups.items():
        if present:
            enriched[group] = True
    
    # Add convenience aliases for common patterns
    # Primary/secondary/tertiary amine detection
    if enriched.get("amine_primary_present"):
        enriched["primary_amine_present"] = True
    if enriched.get("amine_secondary_present"):
        enriched["secondary_amine_present"] = True
    if enriched.get("amine_tertiary_present"):
        enriched["tertiary_amine_present"] = True
    
    # Aniline is both an amine and aromatic
    if enriched.get("aniline_present"):
        enriched["primary_amine_present"] = True
    
    # Amide convenience tokens
    if enriched.get("amide_primary_present"):
        enriched["primary_amide_present"] = True
    if enriched.get("amide_secondary_present"):
        enriched["secondary_amide_present"] = True
    
    # General amine nucleophile token (any N-nucleophile)
    if any(enriched.get(k) for k in [
        "amine_primary_present", "amine_secondary_present", "amine_tertiary_present",
        "aniline_present", "primary_amine_present", "secondary_amine_present"
    ]):
        enriched["amine_nucleophile_present"] = True
    
    # Alcohol/phenol/thiol tokens
    if enriched.get("alcohol_present") or enriched.get("phenol_present"):
        enriched["o_nucleophile_present"] = True
    if enriched.get("thiol_present"):
        enriched["s_nucleophile_present"] = True
    
    # Carbonyl tokens for reductive amination
    if enriched.get("aldehyde_present") or enriched.get("ketone_present"):
        enriched["carbonyl_present"] = True
    
    # SNAr-specific electrophile detection
    # Activated aryl halides or heteroaryls suitable for SNAr
    is_heteroaryl = enriched.get("heteroaryl", False) or enriched.get("heteroaryl_present", False)
    has_aryl_halide = any(enriched.get(k) for k in [
        "aryl_chloride_present", "aryl_bromide_present", 
        "aryl_fluoride_present", "aryl_iodide_present",
        "aryl_halide_present"  # Also check general aryl halide token
    ])
    has_strong_ewg = any(enriched.get(k) for k in [
        "nitro_present", "nitrile_present", "carbonyl_present"
    ])
    
    # Aromatic electrophile present if heteroaryl or aryl halide
    if is_heteroaryl or has_aryl_halide:
        enriched["aromatic_electrophile_present"] = True
    
    # SNAr applicable if heteroaryl (intrinsically activated) or aryl halide with EWG
    if is_heteroaryl or (has_aryl_halide and has_strong_ewg):
        enriched["snar_applicable_electrophile_present"] = True
    # Special case: heteroaryl halides are ALWAYS SNAr-applicable even without additional EWGs
    if is_heteroaryl and has_aryl_halide:
        enriched["snar_applicable_electrophile_present"] = True
    
    # Aryl halide general token
    if has_aryl_halide:
        enriched["aryl_halide_present"] = True
    
    # sp2 halide token (for Suzuki/cross-coupling)
    if has_aryl_halide or enriched.get("vinyl_halide_present"):
        enriched["sp2_halide_present"] = True
    
    # Boron reagent tokens (for Suzuki)
    if any(enriched.get(k) for k in ["boronic_acid_present", "boronic_ester_present"]):
        enriched["sp2_boron_present"] = True
    
    return enriched


def featurize_flat(
    electrophile: str,
    nucleophile: str,
    *,
    include_molpipeline: Optional[bool] = None,
    include_calculable: Optional[bool] = None,
) -> Dict[str, Any]:
    """Return the legacy flat feature dictionary for a reactant pair."""
    base = _featurize_cached(electrophile, nucleophile)
    out = dict(base)

    elec_groups = out.pop("_elec_groups", None)
    nuc_groups = out.pop("_nuc_groups", None)
    if elec_groups is None:
        elec_groups = detect_functional_groups(electrophile)
    if nuc_groups is None:
        nuc_groups = detect_functional_groups(nucleophile)

    # Enrich with comprehensive functional group detection.
    out = _enrich_with_functional_groups(out, electrophile, nucleophile, elec_groups, nuc_groups)

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

    if include_molpipeline is None:
        include_molpipeline = _get_env_bool(
            "CHEMTOOLS_INCLUDE_MOLPIPELINE_FEATURES",
            _maybe_import_molpipeline() is not None,
        )
    if include_molpipeline:
        molpipeline_data: Dict[str, Any] = {}
        settings_snapshot: Optional[Dict[str, Any]] = None

        elec_vec = _molpipeline_vectors(electrophile) if electrophile else None
        if elec_vec is not None:
            settings_snapshot = dict(elec_vec.pop("_settings", {}))
            molpipeline_data["electrophile"] = elec_vec

        nuc_vec = _molpipeline_vectors(nucleophile) if nucleophile else None
        if nuc_vec is not None:
            nuc_settings = dict(nuc_vec.pop("_settings", {}))
            if settings_snapshot is None:
                settings_snapshot = nuc_settings
            else:
                for key, value in nuc_settings.items():
                    settings_snapshot.setdefault(key, value)
            molpipeline_data["nucleophile"] = nuc_vec

        if molpipeline_data:
            if settings_snapshot is not None:
                molpipeline_data["settings"] = settings_snapshot
            out["molpipeline"] = molpipeline_data

    # Standardized substrate-tag tokens (expandable).
    ortho = out.get("ortho_count")
    if isinstance(ortho, str):
        out["ortho_0"] = ortho == "0"
        out["ortho_1"] = ortho == "1"
        out["ortho_2plus"] = ortho == "2+"
        out["ortho_hindered"] = ortho == "2+"

    elec_group_map = elec_groups if isinstance(elec_groups, dict) else {}
    ewg_on_electrophile = any(
        bool(elec_group_map.get(token))
        for token in ("nitro_present", "nitrile_present", "trifluoromethyl_present", "carbonyl_present")
    )
    out["electron_poor_aryl"] = bool(
        out.get("elec_class") == "aryl"
        and (
            bool(out.get("heteroaryl"))
            or bool(out.get("para_EWG"))
            or ewg_on_electrophile
        )
    )

    nuc_group_map = nuc_groups if isinstance(nuc_groups, dict) else {}
    boron_on_nucleophile = any(
        bool(nuc_group_map.get(token))
        for token in ("sp2_boron_present", "boronic_acid_present", "boronic_ester_present")
    )
    ewg_on_nucleophile = any(
        bool(nuc_group_map.get(token))
        for token in ("nitro_present", "nitrile_present", "trifluoromethyl_present", "carbonyl_present")
    )
    out["electron_poor_boronate"] = bool(boron_on_nucleophile and ewg_on_nucleophile)
    out["electron_poor_boronate_present"] = out["electron_poor_boronate"]

    return out


def featurize_pair(
    electrophile: str,
    nucleophile: str,
    *,
    include_molpipeline: Optional[bool] = None,
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
        include_molpipeline=include_molpipeline,
        include_calculable=include_calculable,
    )

    elec_groups = detect_functional_groups(electrophile)
    nuc_groups = detect_functional_groups(nucleophile)

    if include_structural is None:
        include_structural = _get_env_bool("CHEMTOOLS_INCLUDE_STRUCTURAL_FEATURES", False)

    resolved_include_molpipeline = (
        include_molpipeline if include_molpipeline is not None else ("molpipeline" in flat)
    )
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

    molpipeline_block = flat.get("molpipeline")
    calculable_block = flat.get("calculable")
    role_aware_block = flat.get("role_aware")

    elec_payload: Dict[str, Any] = {
        "smiles": electrophile,
        "features": elec_features,
        "functional_groups": _present_groups(elec_groups),
        "functional_group_map": elec_groups,
    }
    nuc_payload: Dict[str, Any] = {
        "smiles": nucleophile,
        "features": nuc_features,
        "functional_groups": _present_groups(nuc_groups),
        "functional_group_map": nuc_groups,
    }

    if isinstance(molpipeline_block, dict):
        elec_mp = molpipeline_block.get("electrophile")
        nuc_mp = molpipeline_block.get("nucleophile")
        if elec_mp:
            elec_payload["molpipeline"] = elec_mp
        if nuc_mp:
            nuc_payload["molpipeline"] = nuc_mp

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
                "include_molpipeline": resolved_include_molpipeline,
                "include_calculable": resolved_include_calculable,
                "include_structural": include_structural,
            },
        },
    }


def featurize(
    electrophile: str,
    nucleophile: str,
    *,
    include_molpipeline: Optional[bool] = None,
    include_calculable: Optional[bool] = None,
    include_structural: Optional[bool] = None,
) -> Dict[str, Any]:
    """Alias for featurize_pair."""
    return featurize_pair(
        electrophile,
        nucleophile,
        include_molpipeline=include_molpipeline,
        include_calculable=include_calculable,
        include_structural=include_structural,
    )


__all__ = ["featurize", "featurize_pair", "featurize_flat"]
