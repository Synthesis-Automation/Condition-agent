"""
Compatibility detection API for legacy imports.

Tests and older callers import ``chemtools.featurizers.detection`` and expect:
- detect_reaction_types(...)
- optional bond-change matching via use_bond_changes=True
- load_reaction_catalog / clear_bond_change_cache
- coupling product confirmation helpers

This module provides that interface while delegating motif-based detection to
``chemtools.detection``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.detection import detect_reaction_types as _detect_motif_reaction_types
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog
from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.util.smarts_cache import compile_smarts
from chemtools.featurizers.formatters.reaction import format_bond_change_key


@dataclass
class ReactionMatch:
    reaction_type: str
    confidence: float
    slot_evidence: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DetectionResult:
    matches: List[ReactionMatch]
    reacted_motifs: List[str] = field(default_factory=list)
    formed_motifs: List[str] = field(default_factory=list)
    reaction_key: str = ""
    error: Optional[str] = None

    @property
    def top_match(self) -> Optional[ReactionMatch]:
        return self.matches[0] if self.matches else None


def _normalize_reaction_type_name(name: Any) -> str:
    text = str(name or "").strip()
    return text.lower() if text else text


def _normalize_bond_token(token: str) -> str:
    raw = str(token or "").strip()
    if not raw or "-" not in raw:
        return raw
    left, right = [part.strip() for part in raw.split("-", 1)]
    ordered = sorted([left, right])
    return f"{ordered[0]}-{ordered[1]}"


def _parse_bond_key_sections(bond_key: Optional[str]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    if not bond_key:
        return (), ()
    formed: Set[str] = set()
    broken: Set[str] = set()
    for part in str(bond_key).split("|"):
        section = part.strip()
        if section.lower().startswith("form:"):
            payload = section.split(":", 1)[1]
            for token in payload.split(";"):
                norm = _normalize_bond_token(token)
                if norm:
                    formed.add(norm)
        elif section.lower().startswith("break:"):
            payload = section.split(":", 1)[1]
            for token in payload.split(";"):
                norm = _normalize_bond_token(token)
                if norm:
                    broken.add(norm)
    return tuple(sorted(formed)), tuple(sorted(broken))


@lru_cache(maxsize=2048)
def _bond_change_signature(reaction_smiles: str) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    bond_key = format_bond_change_key(reaction_smiles)
    return _parse_bond_key_sections(bond_key)


def _bond_signatures_match(
    left: Tuple[Tuple[str, ...], Tuple[str, ...]],
    right: Tuple[Tuple[str, ...], Tuple[str, ...]],
) -> bool:
    return bool(left[0] or left[1]) and left == right


def clear_bond_change_cache() -> None:
    _bond_change_signature.cache_clear()


def _detect_with_bond_changes(reaction_smiles: str) -> DetectionResult:
    query_sig = _bond_change_signature(reaction_smiles)
    definitions, _alias = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for reaction_id, definition in definitions.items():
        references = list(getattr(definition, "reference_reactions", []) or [])
        for ref in references:
            ref_sig = _bond_change_signature(str(ref))
            if not _bond_signatures_match(query_sig, ref_sig):
                continue
            matches.append(
                ReactionMatch(
                    reaction_type=str(reaction_id),
                    confidence=1.0,
                    slot_evidence={
                        "bond_changes": {
                            "formed": list(query_sig[0]),
                            "broken": list(query_sig[1]),
                        }
                    },
                )
            )
            break

    return DetectionResult(matches=matches)


def detect_reaction_types(
    reaction_smiles: str,
    *,
    use_bond_changes: bool = False,
) -> DetectionResult:
    if use_bond_changes:
        result = _detect_with_bond_changes(reaction_smiles)
        if result.matches:
            return result

    motif_result = _detect_motif_reaction_types(reaction_smiles)
    compat_matches: List[ReactionMatch] = []
    for match in motif_result.matches:
        slot_evidence: Dict[str, Any] = {}
        electrophile = getattr(match, "electrophile", None) or []
        nucleophile = getattr(match, "nucleophile", None) or []
        product = getattr(match, "product", None) or []
        if electrophile:
            slot_evidence["electrophile"] = list(electrophile)
        if nucleophile:
            slot_evidence["nucleophile"] = list(nucleophile)
        if product:
            slot_evidence["product"] = list(product)
        compat_matches.append(
            ReactionMatch(
                reaction_type=_normalize_reaction_type_name(getattr(match, "reaction_type", "")),
                confidence=float(getattr(match, "confidence", 0.0)),
                slot_evidence=slot_evidence,
            )
        )
    return DetectionResult(
        matches=compat_matches,
        reacted_motifs=list(getattr(motif_result, "reacted_motifs", []) or []),
        formed_motifs=list(getattr(motif_result, "formed_motifs", []) or []),
        reaction_key=str(getattr(motif_result, "reaction_key", "") or ""),
        error=getattr(motif_result, "error", None),
    )


def detect_reaction_type(
    reaction_smiles: str,
    *,
    use_bond_changes: bool = False,
) -> DetectionResult:
    return detect_reaction_types(reaction_smiles, use_bond_changes=use_bond_changes)


def _iter_valid_mols(smiles_list: Iterable[str]) -> List[Any]:
    mols: List[Any] = []
    for smiles in smiles_list:
        mol = parse_smiles(str(smiles or "").strip())
        if mol is not None:
            mols.append(mol)
    return mols


def _has_any_pattern(mol: Any, smarts_patterns: Iterable[str]) -> bool:
    for smarts in smarts_patterns:
        pattern = compile_smarts(smarts, validate=False)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False


def _has_suzuki_partners(reactants: Iterable[str]) -> bool:
    reactant_mols = _iter_valid_mols(reactants)
    has_halide = any(
        _has_any_pattern(mol, ("[c][Cl,Br,I]", "[CX4][Cl,Br,I]"))
        for mol in reactant_mols
    )
    has_boron = any(
        _has_any_pattern(mol, ("[B](O)O", "[B](O)(O)", "[B](OC)(OC)"))
        for mol in reactant_mols
    )
    return has_halide and has_boron


def _has_biaryl_c_c_link(mol: Any) -> bool:
    ring_info = mol.GetRingInfo()
    ring_sets = [set(ring) for ring in ring_info.AtomRings()]
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() != 6 or b.GetAtomicNum() != 6:
            continue
        if not (a.GetIsAromatic() and b.GetIsAromatic()):
            continue
        a_idx = a.GetIdx()
        b_idx = b.GetIdx()
        same_ring = any((a_idx in ring and b_idx in ring) for ring in ring_sets)
        if not same_ring:
            return True
    return False


def confirm_suzuki_product_by_attachment(
    reactants: Iterable[str],
    products: Iterable[str],
) -> Tuple[bool, str]:
    if not rdkit_available():
        return False, "candidate_build_failed"
    if not _has_suzuki_partners(reactants):
        return False, "candidate_build_failed"
    for product_mol in _iter_valid_mols(products):
        if _has_biaryl_c_c_link(product_mol):
            return True, "substructure_match"
    return False, "no_substructure_match"


def _has_cn_link(mol: Any) -> bool:
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        pair = {a.GetAtomicNum(), b.GetAtomicNum()}
        if pair != {6, 7}:
            continue
        carbon = a if a.GetAtomicNum() == 6 else b
        if carbon.GetIsAromatic():
            return True
    return False


def confirm_coupling_product_by_attachment(
    reactants: Iterable[str],
    products: Iterable[str],
    reaction_type: str,
) -> Tuple[bool, str]:
    if not rdkit_available():
        return False, "candidate_build_failed"
    rtype = str(reaction_type or "").lower()
    if "suzuki" in rtype:
        return confirm_suzuki_product_by_attachment(reactants, products)
    if "c_n" in rtype or "c-n" in rtype or "cn" in rtype:
        for product_mol in _iter_valid_mols(products):
            if _has_cn_link(product_mol):
                return True, "substructure_match"
        return False, "no_substructure_match"
    return False, "candidate_build_failed"


__all__ = [
    "ReactionMatch",
    "DetectionResult",
    "detect_reaction_type",
    "detect_reaction_types",
    "load_reaction_catalog",
    "clear_bond_change_cache",
    "confirm_suzuki_product_by_attachment",
    "confirm_coupling_product_by_attachment",
]
