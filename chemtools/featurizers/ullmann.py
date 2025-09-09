from typing import Dict, Any
from functools import lru_cache
import re

from ..util.rdkit_helpers import rdkit_available, parse_smiles

# Optional role-aware featurization (graceful fallback if package missing)
try:
    from chem_feats import featurize_mol as _role_feat
    _HAS_ROLE_FEATS = True
except Exception:
    _HAS_ROLE_FEATS = False


def _guess_lg_text(s: str) -> str:
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
    t = (s or "").lower()
    if "indole" in t:
        return "indole"
    if re.search(r"c[^)]*n", t):
        return "aniline"
    if "n(" in t:
        return "amine_secondary"
    if "n" in t:
        return "amine_primary"
    if "o" in t and "n" not in t:
        return "phenol"
    return "amine"


def _nucleophile_features(mol, text: str) -> Dict[str, Any]:
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


def featurize(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    # Wrapper around cached implementation; returns a copy to avoid accidental mutation of cache entry
    feat = _featurize_cached(electrophile, nucleophile)
    out = dict(feat)
    # Attach role-aware vectors for downstream models/UI if available
    if _HAS_ROLE_FEATS:
        try:
            elec_ra = _role_feat(electrophile, roles=["aryl_halide"])  # type: ignore[arg-type]
            # nucleophile may be amine or alcohol; compute both roles
            nuc_ra = _role_feat(nucleophile, roles=["amine", "alcohol"])  # type: ignore[arg-type]
            out["role_aware"] = {
                "electrophile": elec_ra,
                "nucleophile": nuc_ra,
            }
        except Exception:
            pass
    return out
