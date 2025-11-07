from __future__ import annotations

"""Helpers for converting rule feature packs into selector matcher payloads."""

from dataclasses import dataclass
from typing import Any, Dict, Mapping

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Mol

from .util.smarts_cache import compile_smarts


# SMARTS patterns for functional group detection
# Compiled lazily via centralized cache
_PHENOL_SMARTS = "[cX3][OX2H]"
_FREE_ALCOHOL_SMARTS = "[CX4;!$([CX3](=O)[OX2H0-,OX1-])][OX2H]"
_CARBOXYLIC_ACID_SMARTS = "C(=O)[OX2H1]"


@dataclass(slots=True)
class _AlphaAmineInfo:
    is_alpha_amino_acid: bool
    internal_amine_free: bool
    has_beta_heteroatom: bool


def _to_mapping(data: Any) -> Mapping[str, Any]:
    return data if isinstance(data, Mapping) else {}


def _mol_from_smiles(smiles: Any) -> Mol | None:
    text = str(smiles).strip() if isinstance(smiles, str) else ""
    if not text:
        return None
    try:
        return Chem.MolFromSmiles(text)
    except Exception:
        return None


def _alpha_amino_context(mol: Mol | None) -> _AlphaAmineInfo:
    if mol is None:
        return _AlphaAmineInfo(False, False, False)

    # Compile carboxylic acid pattern lazily
    carboxylic_acid_pattern = compile_smarts(_CARBOXYLIC_ACID_SMARTS, validate=False)
    if carboxylic_acid_pattern is None:
        return _AlphaAmineInfo(False, False, False)

    is_alpha_amino = False
    amine_free = False
    beta_hetero = False

    try:
        matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    except Exception:
        matches = ()

    for match in matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        for alpha_atom in carbonyl_atom.GetNeighbors():
            if alpha_atom.GetIdx() in match[1:]:
                continue
            if alpha_atom.GetAtomicNum() != 6:
                continue
            # Check for hetero atoms attached at beta position
            for beta_atom in alpha_atom.GetNeighbors():
                if beta_atom.GetIdx() == carbonyl_idx:
                    continue
                if beta_atom.GetAtomicNum() in {7, 8, 15, 16}:
                    beta_hetero = True
                if beta_atom.GetAtomicNum() != 7:
                    continue
                is_alpha_amino = True
                if _is_free_internal_amine(beta_atom):
                    amine_free = True
    return _AlphaAmineInfo(is_alpha_amino, amine_free, beta_hetero)


def _is_free_internal_amine(atom: Atom) -> bool:
    if atom.GetAtomicNum() != 7:
        return False
    # Exclude amide / carbamate nitrogens attached to carbonyl or sulfonyl groups
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 6:
            for sub_bond in other.GetBonds():
                if sub_bond.GetOtherAtom(other).GetAtomicNum() == 8 and sub_bond.GetBondTypeAsDouble() >= 2.0:
                    return False
        if other.GetAtomicNum() == 16:
            return False
    if atom.GetFormalCharge() > 1:
        return False
    # Consider ammonium (NH3+) as free; require at least one hydrogen
    return atom.GetTotalNumHs(includeNeighbors=True) > 0


def _any_substructure(mol: Mol | None, pattern: Mol | None) -> bool:
    if mol is None or pattern is None:
        return False
    try:
        return mol.HasSubstructMatch(pattern)
    except Exception:
        return False


def _has_phenol(smiles: str) -> bool:
    pattern = compile_smarts(_PHENOL_SMARTS, validate=False)
    return _any_substructure(_mol_from_smiles(smiles), pattern)


def _has_free_alcohol(smiles: str) -> bool:
    pattern = compile_smarts(_FREE_ALCOHOL_SMARTS, validate=False)
    return _any_substructure(_mol_from_smiles(smiles), pattern)


def _choose_acid_class(acid_ctx: Mapping[str, Any], rule_features: Mapping[str, Any]) -> str:
    classes = {str(item).lower() for item in acid_ctx.get("classes", []) if isinstance(item, str)}
    sterics = str(acid_ctx.get("sterics") or "").lower()
    if bool(acid_ctx.get("alpha_chiral")):
        return "alpha-chiral"
    if sterics == "hindered":
        return "hindered"
    if bool(acid_ctx.get("heteroaromatic")):
        return "heteroaromatic"
    if "electron" in " ".join(classes):
        return "electron-poor"
    if bool(acid_ctx.get("aromatic")):
        return "benzoic"
    if bool(acid_ctx.get("aliphatic")):
        return "aliphatic"
    return "aliphatic"


def _looks_sulfonamide(smiles: str) -> bool:
    text = smiles.lower()
    return "s(=o)(=o)" in text or "[s](=o)(=o)" in text


def _looks_hydroxylamine(smiles: str) -> bool:
    text = smiles.replace(" ", "").upper()
    return "NOH" in text or "HON" in text


def _select_amine_atom(smiles: str) -> Atom | None:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not atoms:
        return None
    return max(atoms, key=lambda atom: (atom.GetTotalDegree(), atom.GetIdx()))


def _has_aromatic_neighbor(atom: Atom | None) -> bool:
    if atom is None:
        return False
    return any(neigh.GetIsAromatic() and neigh.GetAtomicNum() == 6 for neigh in atom.GetNeighbors())


def _neighbor_carbon_with_aromatic(atom: Atom | None) -> bool:
    if atom is None:
        return False
    for neigh in atom.GetNeighbors():
        if neigh.GetAtomicNum() != 6:
            continue
        for grand in neigh.GetNeighbors():
            if grand.GetIdx() == atom.GetIdx():
                continue
            if grand.GetIsAromatic() and grand.GetAtomicNum() == 6:
                return True
    return False


def _choose_amine_class(rule_features: Mapping[str, Any], amine_ctx: Mapping[str, Any], smiles: str) -> str:
    nuc_class = str(rule_features.get("nuc_class") or amine_ctx.get("nuc_class") or "").lower()
    if bool(amine_ctx.get("diamine")):
        return "diamine"
    if nuc_class == "aniline" and _has_aromatic_neighbor(_select_amine_atom(smiles)):
        return "aniline"
    if _looks_sulfonamide(smiles) or "sulfonamide" in " ".join({str(c).lower() for c in amine_ctx.get("classes", []) if isinstance(c, str)}):
        return "sulfonamide"
    if nuc_class == "hydroxylamine" or _looks_hydroxylamine(smiles):
        return "hydroxylamine"
    atom = _select_amine_atom(smiles)
    if atom is not None:
        carbon_neighbors = sum(1 for neigh in atom.GetNeighbors() if neigh.GetAtomicNum() == 6)
        if carbon_neighbors >= 2:
            return "secondary_aliphatic"
    return "primary_aliphatic"


def _amine_sterics(rule_features: Mapping[str, Any], amine_ctx: Mapping[str, Any]) -> str:
    steric = str(rule_features.get("steric_alpha") or amine_ctx.get("sterics") or "").lower()
    if steric in {"very_high", "very high", "high", "hindered"}:
        return "hindered"
    if steric in {"med", "medium", "moderate"}:
        return "moderate"
    return "unhindered"


def _amine_nucleophilicity(amine_class: str, sterics: str, rule_features: Mapping[str, Any], amine_ctx: Mapping[str, Any]) -> str:
    base = str(rule_features.get("n_basicity") or amine_ctx.get("basicity") or "").lower()
    if amine_class in {"aniline", "sulfonamide", "hydroxylamine"}:
        level = "poor"
    elif amine_class == "diamine":
        level = "moderate"
    else:
        level = "good"
    if base == "deactivated":
        level = "poor"
    if sterics == "hindered" and level == "good":
        return "moderate"
    return level


def _amine_is_benzylic(smiles: str) -> bool:
    atom = _select_amine_atom(smiles)
    if atom is None:
        return False
    if _has_aromatic_neighbor(atom):
        return False  # anilines are handled separately
    return _neighbor_carbon_with_aromatic(atom)


def _functional_group_flags(rule_features: Mapping[str, Any], acid_smiles: str, amine_smiles: str) -> Dict[str, bool]:
    acid_flags = _to_mapping(rule_features.get("functional_groups"))
    alcohol = _has_free_alcohol(acid_smiles) or _has_free_alcohol(amine_smiles)
    phenol = _has_phenol(acid_smiles) or _has_phenol(amine_smiles)
    acid_sensitive = bool(acid_flags.get("base_sensitive"))
    return {
        "alcohols_free": bool(alcohol),
        "phenols_free": bool(phenol),
        "acid_sensitive_motifs": acid_sensitive,
    }


def _infer_need_green(rule_features: Mapping[str, Any]) -> bool:
    category = str(rule_features.get("category") or "").lower()
    goal = str(rule_features.get("goal") or "").lower()
    if any(token in category for token in ("biocatalytic", "enzymatic", "aqueous")):
        return True
    if "green" in goal or "biocatalytic" in goal:
        return True
    return False


def _infer_water_tolerance(rule_features: Mapping[str, Any]) -> bool:
    water_plan = str(rule_features.get("water_management") or "").lower()
    if any(token in water_plan for token in ("aqueous", "buffer", "water")):
        return True
    category = str(rule_features.get("category") or "").lower()
    return "biocatalytic" in category


def build_amide_selector_payload(
    rule_features: Mapping[str, Any],
    *,
    scale_g: float | int | None = None,
    need_green: bool | None = None,
    water_tolerance_needed: bool | None = None,
) -> Dict[str, Any]:
    """Translate amide rule features into a selector payload respecting the schema."""

    acid_ctx = _to_mapping(rule_features.get("acid"))
    amine_ctx = _to_mapping(rule_features.get("amine"))
    electrophile_ctx = _to_mapping(rule_features.get("electrophile"))
    nucleophile_ctx = _to_mapping(rule_features.get("nucleophile"))

    acid_smiles = str(acid_ctx.get("smiles") or electrophile_ctx.get("smiles") or "")
    amine_smiles = str(amine_ctx.get("smiles") or nucleophile_ctx.get("smiles") or "")

    alpha_info = _alpha_amino_context(_mol_from_smiles(acid_smiles))

    acid_section = {
        "class": _choose_acid_class(acid_ctx or electrophile_ctx, rule_features),
        "is_alpha_amino_acid": alpha_info.is_alpha_amino_acid,
        "internal_amine_free": alpha_info.internal_amine_free,
    }
    if alpha_info.has_beta_heteroatom:
        acid_section["has_beta_heteroatom"] = True

    amine_class = _choose_amine_class(rule_features, amine_ctx, amine_smiles)
    amine_sterics = _amine_sterics(rule_features, amine_ctx)
    amine_section = {
        "class": amine_class,
        "sterics": amine_sterics,
        "nucleophilicity": _amine_nucleophilicity(amine_class, amine_sterics, rule_features, amine_ctx),
        "is_benzylic": _amine_is_benzylic(amine_smiles),
    }

    functional_groups = _functional_group_flags(rule_features, acid_smiles, amine_smiles)

    process_section = {
        "scale_g": float(scale_g) if scale_g is not None else 10.0,
        "need_green": bool(need_green) if need_green is not None else _infer_need_green(rule_features),
        "water_tolerance_needed": bool(water_tolerance_needed)
        if water_tolerance_needed is not None
        else _infer_water_tolerance(rule_features),
    }

    return {
        "acid": acid_section,
        "amine": amine_section,
        "functional_groups": functional_groups,
        "process": process_section,
    }


__all__ = ["build_amide_selector_payload"]
