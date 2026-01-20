"""
Motif-based reaction type detection for taxonomy v2.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available

from .motif_detect import detect_motifs
from .motif_registry import build_compound_registry
from ..taxonomy.reaction_catalog import ReactionTypeDefinition, load_reaction_catalog


@dataclass(frozen=True)
class ReactionMatch:
    reaction_type: str
    name: str
    category: Optional[str]
    slot_evidence: Dict[str, List[str]]
    matched_slots: int
    required_slots: int

    @property
    def confidence(self) -> float:
        if self.required_slots == 0:
            return 0.0
        return self.matched_slots / self.required_slots

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_type": self.reaction_type,
            "name": self.name,
            "category": self.category,
            "confidence": self.confidence,
            "matched_slots": self.matched_slots,
            "required_slots": self.required_slots,
            "slot_evidence": {slot: list(values) for slot, values in self.slot_evidence.items()},
        }


@dataclass(frozen=True)
class ReactionDetectionResult:
    matches: List[ReactionMatch]
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        payload = {"matches": [match.to_dict() for match in self.matches]}
        if self.error:
            payload["error"] = self.error
        return payload


BondChangeSignature = Dict[str, Set[Tuple[str, str]]]


def _has_atom_mapping(smiles: str) -> bool:
    import re

    return bool(re.search(r"\[\w+:\d+\]", smiles))


def _build_mapnum_to_element(mapped_smiles: str) -> Dict[int, str]:
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
    return label.split()[0].strip()


def _bond_change_signature(
    *,
    broken_bonds: Iterable[Tuple[Any, Any]],
    formed_bonds: Iterable[Tuple[Any, Any]],
    mapped_smiles: str,
) -> Optional[BondChangeSignature]:
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
    evidence = [f"formed:{a}-{b}" for a, b in sorted(formed)]
    evidence.extend([f"broken:{a}-{b}" for a, b in sorted(broken)])
    return evidence


def _signature_from_reaction(reaction_smiles: str) -> Tuple[Optional[BondChangeSignature], Optional[str], Optional[str]]:
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
    _load_bond_change_references.cache_clear()


@lru_cache(maxsize=1)
def _load_compound_registry() -> Dict[str, Any]:
    base = Path(__file__).resolve().parent.parent / "taxonomy" / "data"
    registry_paths = {
        "groups": base / "organic_groups.v1.3.json",
        "compounds": base / "organic_compounds.v1.3.json",
    }
    return build_compound_registry(registry_paths)


def detect_reaction_types(
    reaction_smiles: str,
    *,
    max_hits_per_compound: Optional[int] = None,
    use_bond_changes: bool = False,
    bond_change_threshold: float = 0.4,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a reaction SMILES string.
    """
    if not reaction_smiles:
        return ReactionDetectionResult(matches=[], error="empty_reaction")
    if use_bond_changes:
        bond_change_result = _detect_reaction_types_by_bond_changes(
            reaction_smiles,
            min_similarity=bond_change_threshold,
        )
        if bond_change_result.matches:
            return bond_change_result
    from chemtools.smiles import normalize_reaction
    normalized = normalize_reaction(reaction_smiles)
    reactants = normalized.get("reactants") or []
    reactant_smiles = [
        item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
        for item in reactants
    ]
    reactant_smiles = [s for s in reactant_smiles if s]

    # Only consider products if the reaction SMILES explicitly contains them (via > or >>)
    product_smiles = None
    if ">" in reaction_smiles:
        product_smiles = [
            item.get("smiles_norm") or item.get("largest_smiles") or item.get("input") or ""
            for item in (normalized.get("products") or [])
        ]
        product_smiles = [s for s in product_smiles if s]

    if not reactant_smiles:
        return ReactionDetectionResult(matches=[], error="no_reactants")
    return detect_reaction_types_from_smiles(
        reactant_smiles,
        product_smiles=product_smiles,
        max_hits_per_compound=max_hits_per_compound,
    )


def detect_reaction_types_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    product_smiles: Optional[Iterable[str]] = None,
    max_hits_per_compound: Optional[int] = None,
) -> ReactionDetectionResult:
    """
    Detect reaction types from a list of reactant SMILES strings.
    """
    if not rdkit_available():
        return ReactionDetectionResult(matches=[], error="rdkit_unavailable")

    detected_profile = _detect_motif_profile(
        reactant_smiles, max_hits_per_compound=max_hits_per_compound
    )
    if not detected_profile:
        return ReactionDetectionResult(matches=[])

    product_profile = (
        _detect_motif_profile(product_smiles, max_hits_per_compound=max_hits_per_compound)
        if product_smiles is not None
        else None
    )

    definitions, _ = load_reaction_catalog()
    matches: List[ReactionMatch] = []

    for definition in definitions.values():
        match = _match_reaction_definition(definition, detected_profile, product_profile)
        if match is not None:
            matches.append(match)

    # Priority map to break ties between overlapping definitions (e.g., Suzuki vs C-H Arylation)
    # Higher values come first.
    _PRIORITIES = {
        "suzuki_miyaura": 100,
        "buchwald_hartwig_cn": 100,
        "buchwald_hartwig_co": 100,
        "amide_formation": 100,
        "c_h_arylation": 10,
        "arylation_acidic_c_h": 10,
    }

    matches.sort(
        key=lambda m: (
            -m.matched_slots,
            -_PRIORITIES.get(m.reaction_type, 50),
            -sum(len(v) for v in m.slot_evidence.values()),
            m.reaction_type,
        )
    )
    return ReactionDetectionResult(matches=matches)


def _detect_reaction_types_by_bond_changes(
    reaction_smiles: str,
    *,
    min_similarity: float = 0.4,
) -> ReactionDetectionResult:
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


def detect_motif_ids_from_smiles(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> Set[str]:
    """Return the set of organic-compound motif IDs detected in the reactants."""
    if not rdkit_available():
        return set()
    return set(
        _detect_motif_profile(
            reactant_smiles, max_hits_per_compound=max_hits_per_compound
        ).keys()
    )


def _detect_motif_profile(
    reactant_smiles: Iterable[str],
    *,
    max_hits_per_compound: Optional[int] = None,
) -> Dict[str, Dict[str, Any]]:
    registry = _load_compound_registry()
    compiled = registry["compiled_compounds"]
    detected: Dict[str, Dict[str, Any]] = {}

    for idx, smiles in enumerate(reactant_smiles):
        mol = parse_smiles(smiles)
        if mol is None:
            continue
        
        from rdkit import Chem
        mol = Chem.AddHs(mol)
        
        hits = detect_motifs(
            mol,
            compiled,
            max_hits_per_compound=max_hits_per_compound,
        )
        for hit in hits:
            compound_id = hit.get("compound_id")
            if compound_id:
                entry = detected.setdefault(
                    compound_id,
                    {"count": 0, "molecules": set()},
                )
                entry["count"] += 1
                entry["molecules"].add(idx)

    return detected


def _match_reaction_definition(
    definition: ReactionTypeDefinition,
    detected_motifs: Dict[str, Dict[str, Any]],
    detected_products: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Optional[ReactionMatch]:
    slot_evidence: Dict[str, List[str]] = {}
    required_slots = 0
    matched_slots = 0

    if not definition.reactants and not definition.products:
        return None

    def apply_requirements(
        requirements: Dict[str, Any],
        profile: Dict[str, Dict[str, Any]],
        slot_prefix: str = "",
    ) -> bool:
        nonlocal required_slots, matched_slots
        for slot, requirement in requirements.items():
            allowed = requirement.allowed
            if not allowed:
                continue
            required_slots += 1
            hits = [
                motif
                for motif in allowed
                if profile.get(motif, {}).get("count", 0) > 0
            ]
            if not hits:
                return False
            total_hits = sum(profile.get(motif, {}).get("count", 0) for motif in allowed)
            if total_hits < requirement.min_hits:
                return False
            if requirement.min_reactants > 1:
                molecule_indices = set()
                for motif in allowed:
                    entry = profile.get(motif)
                    if entry:
                        molecule_indices.update(entry.get("molecules") or set())
                if len(molecule_indices) < requirement.min_reactants:
                    return False
            slot_name = f"{slot_prefix}{slot}" if slot_prefix else slot
            slot_evidence[slot_name] = hits
            matched_slots += 1
        return True

    if definition.reactants and not apply_requirements(definition.reactants, detected_motifs):
        return None
    
    # Only check products if they were provided in the input
    if definition.products and detected_products is not None:
        if not apply_requirements(definition.products, detected_products, slot_prefix="product:"):
            return None

    if required_slots == 0:
        return None

    # Apply top-level constraints
    constraints = definition.constraints
    if constraints:
        total_hits = 0
        all_molecule_indices = set()
        for slot_name, hits in slot_evidence.items():
            # Determine if this was a reactant or product slot
            is_product_slot = slot_name.startswith("product:")
            profile = detected_products if is_product_slot else detected_motifs
            
            # If it's a product slot but no products were provided, skip counting
            if is_product_slot and profile is None:
                continue
                
            for motif in hits:
                entry = profile.get(motif) if profile else None
                if entry:
                    total_hits += entry.get("count", 0)
                    all_molecule_indices.update(entry.get("molecules") or set())

        if "min_hits" in constraints and total_hits < int(constraints["min_hits"]):
            return None
        if "min_reactants" in constraints and len(all_molecule_indices) < int(constraints["min_reactants"]):
            return None

    return ReactionMatch(
        reaction_type=definition.id,
        name=definition.name,
        category=definition.category or None,
        slot_evidence=slot_evidence,
        matched_slots=matched_slots,
        required_slots=required_slots,
    )
