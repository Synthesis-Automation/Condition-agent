"""
Molecular featurizer for C-N coupling substrates.

Extracts structural and chemical features from electrophile/nucleophile pairs
for C-N coupling reactions (both Ullmann and Buchwald-Hartwig). Features include
leaving group type, electrophile classification, nucleophile basicity, and
steric properties. Conditionally attaches role-aware vectors when available.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional
from functools import lru_cache
import os
import re

from ..util.rdkit_helpers import rdkit_available, parse_smiles
from ..util.functional_groups import detect_all as detect_functional_groups

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


def _detect_electrophile(mol) -> Dict[str, Any]:
    """Extract electrophile features (LG, class, ortho substitution, EWG, heteroaryl)."""
    # Default values
    res = {
        "LG": "UNK",
        "elec_class": "aryl",
        "ortho_count": "0",
        "para_EWG": False,
        "heteroaryl": False,
    }
    if mol is None or not rdkit_available():
        return res
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return res

    # SMARTS patterns
    patt_aryl_halide = Chem.MolFromSmarts("[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]")
    patt_triflate = Chem.MolFromSmarts("OS(=O)(=O)C(F)(F)F")
    patt_vinyl_halide = Chem.MolFromSmarts("C=C[Cl,Br,I]")
    patt_alkyl_halide = Chem.MolFromSmarts("[CX4][Cl,Br,I]")

    # EWG patterns (anywhere on ring; approximate)
    ewg_smarts = [
        Chem.MolFromSmarts("[N+](=O)[O-]"),
        Chem.MolFromSmarts("C#N"),
        Chem.MolFromSmarts("C(F)(F)F"),
        Chem.MolFromSmarts("C(=O)[!O]")
    ]

    # Determine LG and class
    lg = "UNK"
    elec_class = "aryl"
    ipso_atom = None
    # Triflate attached to oxygen; find directly
    if mol.HasSubstructMatch(patt_triflate):
        lg = "OTf"
        # class heuristic: assume aryl/vinyl via presence of aromatic neighbor
        for match in mol.GetSubstructMatches(patt_triflate):
            # match gives atoms of triflate; find the oxygen index connected outward
            for ai in match:
                a = mol.GetAtomWithIdx(ai)
                if a.GetSymbol() == 'O':
                    for nb in a.GetNeighbors():
                        if nb.GetAtomicNum() == 6:
                            if nb.GetIsAromatic():
                                elec_class = "aryl"
                                ipso_atom = nb
                            else:
                                # approximate vinyl/alkyl
                                elec_class = "vinyl" if any(b.GetAtomicNum()==6 and b.GetIsAromatic()==False and any(x.GetSymbol()=="C" and x.GetIsAromatic()==False for x in nb.GetNeighbors()) for b in [nb]) else "alkyl"
                            break
                    if ipso_atom:
                        break
            if ipso_atom:
                break
    # Halides
    if lg == "UNK":
        for sym in ("I", "Br", "Cl"):
            patt = Chem.MolFromSmarts(f"[c,C][{sym}]")
            if mol.HasSubstructMatch(patt):
                lg = sym
                # pick a match
                mi = mol.GetSubstructMatch(patt)
                if mi:
                    c_idx = mi[0]
                    c_atom = mol.GetAtomWithIdx(c_idx)
                    ipso_atom = c_atom
                    elec_class = "aryl" if c_atom.GetIsAromatic() else "alkenyl" if any(b.GetIsAromatic()==False and b.GetTotalDegree()==3 for b in c_atom.GetNeighbors()) else "alkyl"
                break

    # Ortho count and heteroaryl using a 6-member aromatic ring if present
    ortho_count = 0
    heteroaryl = False
    if ipso_atom is not None and ipso_atom.GetIsAromatic():
        ri = ipso_atom.GetOwningMol().GetRingInfo()
        atom_idx = ipso_atom.GetIdx()
        rings = [r for r in ri.AtomRings() if atom_idx in r and len(r) == 6]
        if rings:
            ring = rings[0]
            pos = ring.index(atom_idx)
            # ortho positions are pos-1 and pos+1
            ortho_atoms = [mol.GetAtomWithIdx(ring[(pos - 1) % 6]), mol.GetAtomWithIdx(ring[(pos + 1) % 6])]
            # heteroaryl if any atom in ring is hetero
            heteroaryl = any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6, 1) for i in ring)
            def is_substituted(ar_atom):
                # count non-ring heavy neighbors not in ring
                cnt = 0
                for nb in ar_atom.GetNeighbors():
                    if nb.GetIdx() not in ring and nb.GetAtomicNum() > 1:
                        cnt += 1
                return cnt > 0
            ortho_count = sum(1 for a in ortho_atoms if is_substituted(a))

    # para EWG approximation: check presence of any EWG anywhere on molecule
    para_ewg = False
    for patt in ewg_smarts:
        try:
            if mol.HasSubstructMatch(patt):
                para_ewg = True
                break
        except Exception:
            continue

    res.update({
        "LG": lg,
        "elec_class": "aryl" if (ipso_atom is not None and ipso_atom.GetIsAromatic()) else ("vinyl" if mol.HasSubstructMatch(patt_vinyl_halide) else ("alkyl" if mol.HasSubstructMatch(patt_alkyl_halide) else elec_class)),
        "ortho_count": "2+" if ortho_count >= 2 else ("1" if ortho_count == 1 else "0"),
        "para_EWG": para_ewg,
        "heteroaryl": heteroaryl,
    })
    return res


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


def _nucleophile_features(mol, text: str) -> Dict[str, Any]:
    """Extract nucleophile features (class, basicity, sterics)."""
    nuc_class = _nuc_class_text(text)
    n_basicity = "unknown"
    steric_alpha = "low"

    if mol is None or not rdkit_available():
        if nuc_class == "aniline":
            n_basicity = "aromatic_primary"
        elif nuc_class == "amine_primary":
            n_basicity = "aliphatic_primary"
        elif nuc_class == "amine_secondary":
            n_basicity = "secondary"
        elif nuc_class == "phenol":
            n_basicity = "deactivated"
        return {"nuc_class": nuc_class, "n_basicity": n_basicity, "steric_alpha": steric_alpha}

    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return {"nuc_class": nuc_class, "n_basicity": n_basicity, "steric_alpha": steric_alpha}

    # SMARTS patterns
    patt_aniline = Chem.MolFromSmarts("[$([NX3;H1][c]),$([NX3]([c])[H])]")
    patt_indole = Chem.MolFromSmarts("[nH]")
    patt_phenol = Chem.MolFromSmarts("c[OH]")
    patt_amide = Chem.MolFromSmarts("N[C;X3](=O)")
    patt_sec_amine = Chem.MolFromSmarts("[NX3;H1]([!#6])([#6])")
    patt_prim_amine = Chem.MolFromSmarts("[NX3;H2][#6]")

    if mol.HasSubstructMatch(patt_indole):
        nuc_class = "indole"
        n_basicity = "deactivated"
    elif mol.HasSubstructMatch(patt_aniline):
        nuc_class = "aniline"
        n_basicity = "aromatic_primary"
    elif mol.HasSubstructMatch(patt_phenol):
        nuc_class = "phenol"
        n_basicity = "deactivated"
    elif mol.HasSubstructMatch(patt_amide):
        nuc_class = "amide_deactivated"
        n_basicity = "deactivated"
    elif mol.HasSubstructMatch(patt_sec_amine):
        nuc_class = "amine_secondary"
        n_basicity = "secondary"
    elif mol.HasSubstructMatch(patt_prim_amine):
        nuc_class = "amine_primary"
        n_basicity = "aliphatic_primary"
    else:
        # Fallback categorization
        if "n" in text.lower():
            nuc_class = "amine_primary"
            n_basicity = "aliphatic_primary"
        elif "o" in text.lower():
            nuc_class = "phenol"
            n_basicity = "deactivated"

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

    # Electrophile features
    if emol is not None and rdkit_available():
        elec = _detect_electrophile(emol)
    else:
        elec = {
            "LG": _guess_lg_text(electrophile),
            "elec_class": "aryl",
            "ortho_count": "0",
            "para_EWG": any(x in (electrophile or '').lower() for x in ("[n+](=o)[o-]", "c#n", "c(f)(f)f")),
            "heteroaryl": False,
        }

    # Nucleophile features
    nuc = _nucleophile_features(nmol, nucleophile)

    # Compose bin key (coarse)
    lg = elec.get("LG", "UNK")
    nuc_class = nuc.get("nuc_class", "unknown")
    bin_key = f"LG:{lg}|NUC:{nuc_class}"

    return {
        **elec,
        **nuc,
        "bin": bin_key,
    }


def _enrich_with_functional_groups(
    features: Dict[str, Any],
    electrophile: str,
    nucleophile: str
) -> Dict[str, Any]:
    """
    Enrich feature dict with comprehensive functional group detection.
    
    Adds <group_name>_present boolean tokens for all detected functional groups,
    plus special tokens for SNAr and reductive amination reactions.
    
    Args:
        features: Base feature dictionary to enrich
        electrophile: SMILES string of electrophile
        nucleophile: SMILES string of nucleophile
        
    Returns:
        Enriched feature dictionary with functional group tokens
    """
    enriched = dict(features)
    
    # Detect functional groups in both reactants
    if electrophile:
        elec_groups = detect_functional_groups(electrophile)
        for group, present in elec_groups.items():
            if present:
                enriched[f"{group}_present"] = True
    
    if nucleophile:
        nuc_groups = detect_functional_groups(nucleophile)
        for group, present in nuc_groups.items():
            if present:
                enriched[f"{group}_present"] = True
    
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


def featurize(
    electrophile: str,
    nucleophile: str,
    *,
    include_molpipeline: Optional[bool] = None,
    include_calculable: Optional[bool] = None,
) -> Dict[str, Any]:
    """Featurize electrophile/nucleophile pair for C-N coupling reactions.
    
    Extracts structural and chemical features including:
    - Leaving group type (Br, Cl, I, OTf)
    - Electrophile classification (aryl, vinyl, alkyl)
    - Ortho substitution count
    - Para electron-withdrawing groups
    - Nucleophile class and basicity
    - Steric hindrance
    
    Optionally attaches role-aware feature vectors when CHEMTOOLS_ATTACH_ROLE_AWARE=1.
    Optionally includes calculable features when CHEMTOOLS_INCLUDE_CALCULABLE_FEATURES=1.
    
    Args:
        electrophile: SMILES string of electrophile (aryl halide, triflate, etc.)
        nucleophile: SMILES string of nucleophile (amine, alcohol, etc.)
        include_molpipeline: Include MolPipeline features (fingerprints, descriptors)
        include_calculable: Include calculable features from calculable_features.json
        
    Returns:
        Dictionary with extracted features and optional role-aware vectors
    """
    base = _featurize_cached(electrophile, nucleophile)
    out = dict(base)

    # Enrich with comprehensive functional group detection
    # This adds <group>_present tokens needed by rule databases
    out = _enrich_with_functional_groups(out, electrophile, nucleophile)

    # Attach calculable features when explicitly enabled via env or parameter
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

    # Attach role-aware only when explicitly enabled via env (default off for speed)
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

    return out


__all__ = ["featurize"]
