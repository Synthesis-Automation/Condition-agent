"""
Bond change signature analysis for reaction type detection.

Uses atom mapping and bond change patterns to match reactions against reference templates.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

from .models import ReactionMatch, ReactionDetectionResult


BondChangeSignature = Dict[str, Set[Tuple[str, str]]]


def _has_atom_mapping(smiles: str) -> bool:
    """Check if SMILES string contains atom mapping."""
    import re
    return bool(re.search(r"\[\w+:\d+\]", smiles))


def _build_mapnum_to_element(mapped_smiles: str) -> Dict[int, str]:
    """Extract mapping from atom map numbers to element symbols."""
    from rdkit import Chem

    mapping: Dict[int, str] = {}
    for side in mapped_smiles.split(">>"):
        for frag in side.split("."):
            frag = frag.strip()
            if not frag:
                continue
            mol = Chem.MolFromSmiles(frag)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0 and map_num not in mapping:
                    mapping[map_num] = atom.GetSymbol()
    return mapping


def _normalize_element_label(label: str) -> str:
    """Normalize element label to symbol only."""
    return label.split()[0].strip()


def _bond_change_signature(
    *,
    broken_bonds: Iterable[Tuple[Any, Any]],
    formed_bonds: Iterable[Tuple[Any, Any]],
    mapped_smiles: str,
) -> Optional[BondChangeSignature]:
    """Build a bond change signature from broken/formed bonds."""
    map_to_element = _build_mapnum_to_element(mapped_smiles)
    if not map_to_element:
        return None

    def resolve_element(value: Any) -> Optional[str]:
        if isinstance(value, int):
            return map_to_element.get(value)
        if isinstance(value, str):
            return _normalize_element_label(value)
        return None

    def normalize_pair(a: Any, b: Any) -> Optional[Tuple[str, str]]:
        elem_a = resolve_element(a)
        elem_b = resolve_element(b)
        if not elem_a or not elem_b:
            return None
        return tuple(sorted((elem_a, elem_b)))

    formed: Set[Tuple[str, str]] = set()
    broken: Set[Tuple[str, str]] = set()

    for bond in formed_bonds:
        pair = normalize_pair(*bond)
        if pair:
            formed.add(pair)
    for bond in broken_bonds:
        pair = normalize_pair(*bond)
        if pair:
            broken.add(pair)

    if not formed and not broken:
        return None

    return {"formed": formed, "broken": broken}


def _bond_change_similarity(
    left: BondChangeSignature,
    right: BondChangeSignature,
) -> Tuple[float, Set[Tuple[str, str]]]:
    """Calculate Jaccard similarity between two bond change signatures."""
    formed_intersection = left["formed"] & right["formed"]
    broken_intersection = left["broken"] & right["broken"]
    union = left["formed"] | left["broken"] | right["formed"] | right["broken"]
    intersection = formed_intersection | broken_intersection
    if not union:
        return 0.0, set()
    return len(intersection) / len(union), intersection


def _format_bond_changes(
    formed: Iterable[Tuple[str, str]],
    broken: Iterable[Tuple[str, str]],
) -> List[str]:
    """Format bond changes as human-readable evidence strings."""
    evidence = [f"formed:{a}-{b}" for a, b in sorted(formed)]
    evidence.extend([f"broken:{a}-{b}" for a, b in sorted(broken)])
    return evidence


def _signature_from_reaction(
    reaction_smiles: str
) -> Tuple[Optional[BondChangeSignature], Optional[str], Optional[str]]:
    """Extract bond change signature from reaction SMILES."""
    if ">>" not in reaction_smiles:
        return None, None, "missing_products"

    from chemtools._atom_mapping import analyze_bond_changes_hybrid

    analysis = analyze_bond_changes_hybrid(
        reaction_smiles,
        use_rxnmapper=True,
        use_mcs=False,
        use_manual=True,
        auto_map=True,
    )
    if not analysis.get("success"):
        return None, None, analysis.get("error") or "bond_change_unavailable"

    if analysis.get("method") == "mcs_only":
        return None, None, "bond_change_unavailable"

    result = analysis.get("recommended_result") or {}
    mapped_smiles = None
    if _has_atom_mapping(reaction_smiles):
        mapped_smiles = reaction_smiles
    else:
        mapped_smiles = result.get("mapped_smiles")

    if not mapped_smiles:
        return None, None, "mapped_smiles_missing"

    signature = _bond_change_signature(
        broken_bonds=result.get("broken_bonds", []),
        formed_bonds=result.get("formed_bonds", []),
        mapped_smiles=mapped_smiles,
    )
    if not signature:
        return None, mapped_smiles, "empty_signature"

    return signature, mapped_smiles, None


@lru_cache(maxsize=1)
def _load_bond_change_references() -> Dict[str, List[Tuple[BondChangeSignature, str]]]:
    """Load reference bond change signatures from reaction catalog."""
    definitions, _ = load_reaction_catalog()
    references: Dict[str, List[Tuple[BondChangeSignature, str]]] = {}

    for definition in definitions.values():
        if not definition.reference_reactions:
            continue
        for reference in definition.reference_reactions:
            signature, _, error = _signature_from_reaction(reference)
            if signature is None or error:
                continue
            references.setdefault(definition.id, []).append((signature, reference))

    return references


def clear_bond_change_cache() -> None:
    """Clear the cached bond change references."""
    _load_bond_change_references.cache_clear()


def detect_reaction_types_by_bond_changes(
    reaction_smiles: str,
    *,
    min_similarity: float = 0.4,
) -> ReactionDetectionResult:
    """Detect reaction types by comparing bond change signatures."""
    signature, _, error = _signature_from_reaction(reaction_smiles)
    if signature is None:
        return ReactionDetectionResult(matches=[], error=error or "bond_change_unavailable")

    references = _load_bond_change_references()
    if not references:
        return ReactionDetectionResult(matches=[], error="no_reference_reactions")

    definitions, _ = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for reaction_id, ref_signatures in references.items():
        best_similarity = 0.0
        best_reference = None
        best_signature: Optional[BondChangeSignature] = None
        best_intersection: Set[Tuple[str, str]] = set()

        for ref_signature, ref_smiles in ref_signatures:
            similarity, intersection = _bond_change_similarity(signature, ref_signature)
            if similarity > best_similarity:
                best_similarity = similarity
                best_reference = ref_smiles
                best_signature = ref_signature
                best_intersection = intersection

        if best_signature is None or best_similarity < min_similarity:
            continue

        union = (
            signature["formed"]
            | signature["broken"]
            | best_signature["formed"]
            | best_signature["broken"]
        )
        if not union:
            continue

        matched_slots = len(best_intersection)
        required_slots = len(union)
        evidence = _format_bond_changes(
            signature["formed"] & best_signature["formed"],
            signature["broken"] & best_signature["broken"],
        )
        slot_evidence: Dict[str, List[str]] = {"bond_changes": evidence}
        if best_reference:
            slot_evidence["reference_reaction"] = [best_reference]

        definition = definitions.get(reaction_id)
        matches.append(
            ReactionMatch(
                reaction_type=reaction_id,
                name=definition.name if definition else reaction_id,
                category=definition.category if definition else None,
                slot_evidence=slot_evidence,
                matched_slots=matched_slots,
                required_slots=required_slots,
            )
        )

    matches.sort(key=lambda m: (-m.confidence, m.reaction_type))
    return ReactionDetectionResult(matches=matches, error=error if not matches else None)
