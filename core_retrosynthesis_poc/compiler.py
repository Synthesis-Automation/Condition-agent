"""Compile executable retrosynthetic SMARTS from trusted reaction cores."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple, cast

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts
from retrosynthesis_poc.chemistry import (
    canonical_smiles,
    contributing_reactants,
    digest,
    split_reaction_smiles,
)
from retrosynthesis_poc.mapping import materialize_atom_mapping

from .context import context_from_analysis
from .models import (
    AbstractionLevel,
    BondKind,
    CoreTemplate,
    CoreTemplatePrecedent,
)


ALLOWED_CORE_EVIDENCE = frozenset({"verified", "inferred"})
ALLOWED_HETEROATOMS = frozenset({"N", "O", "S"})


@dataclass(frozen=True)
class CompilationResult:
    """Core-derived templates or one deterministic rejection reason."""

    templates: Tuple[CoreTemplate, ...]
    rejection_reason: Optional[str]


def _observation(row: Dict[str, Any]) -> Dict[str, Any]:
    return dict(row.get("reaction_observation") or row.get("observation") or {})


def _source_reaction(row: Dict[str, Any], observation: Dict[str, Any]) -> str:
    return str(
        row.get("reaction_smiles")
        or observation.get("input_reaction_smiles")
        or row.get("input_reaction_smiles")
        or ""
    )


def _reference_id(row: Dict[str, Any]) -> str:
    direct = str(row.get("reference_id") or "").strip()
    if direct:
        return direct
    return str((row.get("reference_identity") or {}).get("reference_id") or "")


def _formed_center(
    observation: Dict[str, Any],
) -> Optional[tuple[BondKind, tuple[int, int]]]:
    formed = [
        edit
        for edit in observation.get("edits") or ()
        if edit.get("edit_type") == "formed" and edit.get("atom_2") is not None
    ]
    if len(formed) != 1:
        return None
    edit = formed[0]
    if str(edit.get("new_order") or "").upper() != "SINGLE":
        return None
    atoms = (edit.get("atom_1") or {}, edit.get("atom_2") or {})
    elements = tuple(str(atom.get("element") or "") for atom in atoms)
    if "C" not in elements:
        return None
    heteroatoms = set(elements).intersection(ALLOWED_HETEROATOMS)
    if len(heteroatoms) != 1:
        return None
    maps = tuple(int(atom.get("atom_map_number") or 0) for atom in atoms)
    if any(map_number <= 0 for map_number in maps):
        return None
    carbon_map = maps[elements.index("C")]
    heteroatom = next(iter(heteroatoms))
    hetero_map = maps[elements.index(heteroatom)]
    return cast(BondKind, f"C-{heteroatom}"), (carbon_map, hetero_map)


def _direct_handle_signature(molecule: Any, center_maps: set[int]) -> str:
    tokens = []
    for atom in molecule.GetAtoms():
        if int(atom.GetAtomMapNum()) not in center_maps:
            continue
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if int(neighbor.GetAtomMapNum()) > 0:
                continue
            tokens.append(
                f"{atom.GetSymbol()}-{neighbor.GetSymbol()}:{bond.GetBondType().name}"
            )
    return ";".join(sorted(tokens)) or "no_observed_precursor_handle"


def _canonical_map_lookup(
    reactants: Any,
    product: Any,
    center_maps: tuple[int, int],
) -> dict[int, int]:
    lookup = {center_maps[0]: 1, center_maps[1]: 2}
    product_atoms = {
        int(atom.GetAtomMapNum()): atom
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    reactant_atoms = {
        int(atom.GetAtomMapNum()): atom
        for atom in reactants.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    product_without_maps = Chem.Mol(product)
    for atom in product_without_maps.GetAtoms():
        atom.SetAtomMapNum(0)
    product_ranks = tuple(
        int(value)
        for value in Chem.CanonicalRankAtoms(
            product_without_maps,
            breakTies=True,
            includeChirality=True,
        )
    )
    product_rank_by_map = {
        int(atom.GetAtomMapNum()): product_ranks[int(atom.GetIdx())]
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    remaining = set(product_atoms).union(reactant_atoms).difference(lookup)

    def key(map_number: int) -> tuple[Any, ...]:
        atom = product_atoms.get(map_number) or reactant_atoms[map_number]
        return (
            product_rank_by_map.get(map_number, 10**9),
            atom.GetSymbol(),
            int(atom.GetIsAromatic()),
            int(atom.GetFormalCharge()),
            int(atom.GetDegree()),
        )

    for canonical, map_number in enumerate(sorted(remaining, key=key), start=3):
        lookup[map_number] = canonical
    return lookup


def _apply_map_lookup(molecule: Any, lookup: dict[int, int]) -> None:
    for atom in molecule.GetAtoms():
        current = int(atom.GetAtomMapNum())
        if current > 0:
            atom.SetAtomMapNum(lookup[current])


def _selected_atoms(
    molecule: Any,
    *,
    level: AbstractionLevel,
    reactant_side: bool,
) -> set[int]:
    centers = {
        int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if int(atom.GetAtomMapNum()) in {1, 2}
    }
    selected = set(centers)
    handle_frontier = []
    for atom_index in centers:
        atom = molecule.GetAtomWithIdx(atom_index)
        for neighbor in atom.GetNeighbors():
            if level == "L2" or (
                reactant_side and int(neighbor.GetAtomMapNum()) == 0
            ):
                selected.add(int(neighbor.GetIdx()))
            if reactant_side and int(neighbor.GetAtomMapNum()) == 0:
                handle_frontier.append(int(neighbor.GetIdx()))
    while handle_frontier:
        atom_index = handle_frontier.pop()
        atom = molecule.GetAtomWithIdx(atom_index)
        for neighbor in atom.GetNeighbors():
            neighbor_index = int(neighbor.GetIdx())
            if int(neighbor.GetAtomMapNum()) > 0 or neighbor_index in selected:
                continue
            selected.add(neighbor_index)
            handle_frontier.append(neighbor_index)
    return selected


def _fragment_smarts(molecule: Any, selected: set[int]) -> str:
    values = []
    for component in Chem.GetMolFrags(molecule, asMols=False, sanitizeFrags=False):
        atoms = sorted(selected.intersection(int(index) for index in component))
        if not atoms:
            continue
        atom_set = set(atoms)
        bonds = [
            int(bond.GetIdx())
            for bond in molecule.GetBonds()
            if int(bond.GetBeginAtomIdx()) in atom_set
            and int(bond.GetEndAtomIdx()) in atom_set
        ]
        values.append(
            Chem.MolFragmentToSmarts(
                molecule,
                atoms,
                bonds,
                isomericSmarts=True,
            )
        )
    return ".".join(sorted(values))


def _round_trip(
    reaction_smarts: str,
    product_smiles: str,
    expected_precursors: str,
) -> bool:
    canonical_product = canonical_smiles(product_smiles)
    if canonical_product is None:
        return False
    try:
        with redirect_stdout(io.StringIO()):
            outcomes = rdchiralRun(
                rdchiralReaction(reaction_smarts),
                rdchiralReactants(canonical_product),
            )
    except Exception:
        return False
    return expected_precursors in {
        canonical
        for outcome in outcomes
        if (canonical := canonical_smiles(str(outcome))) is not None
    }


def _compile_level(
    *,
    level: AbstractionLevel,
    reactants: Any,
    product: Any,
) -> tuple[str, str, str]:
    product_smarts = _fragment_smarts(
        product,
        _selected_atoms(product, level=level, reactant_side=False),
    )
    precursor_smarts = _fragment_smarts(
        reactants,
        _selected_atoms(reactants, level=level, reactant_side=True),
    )
    return product_smarts, precursor_smarts, f"{product_smarts}>>{precursor_smarts}"


def compile_core_templates(
    row: Dict[str, Any],
    *,
    levels: Iterable[AbstractionLevel] = ("L1", "L2"),
) -> CompilationResult:
    """Compile source-round-tripped SMARTS from one trusted reaction core."""

    requested_levels = tuple(dict.fromkeys(levels))
    if not requested_levels or any(level not in {"L1", "L2"} for level in requested_levels):
        raise ValueError("levels must contain L1 and/or L2")
    source_core = dict(row.get("reaction_core") or {})
    source_observation = _observation(row)
    if not source_core or not source_observation:
        return CompilationResult((), "missing_core_or_observation")
    if (source_core.get("quality") or {}).get("status") not in {"pass", "review"}:
        return CompilationResult((), "core_quality_blocked")
    if source_core.get("evidence_status") not in ALLOWED_CORE_EVIDENCE:
        return CompilationResult((), "core_evidence_not_allowed")
    if int(source_core.get("event_count") or 0) != 1:
        return CompilationResult((), "not_single_event")
    completeness = source_observation.get("completeness") or row.get(
        "reaction_completeness"
    ) or {}
    if completeness.get("status") != "verified":
        return CompilationResult((), "product_completeness_not_verified")
    if len(source_observation.get("products") or ()) != 1:
        return CompilationResult((), "not_single_product")
    if source_observation.get("stereo_changes"):
        return CompilationResult((), "stereo_change_out_of_scope")
    if any(
        edit.get("edit_type") == "order_changed"
        for edit in source_observation.get("edits") or ()
    ):
        return CompilationResult((), "bond_order_change_out_of_scope")

    materialized = materialize_atom_mapping(
        _source_reaction(row, source_observation)
    )
    if materialized is None:
        return CompilationResult((), "atom_mapping_unavailable")
    mapped_analysis = featurize_reaction(materialized.reaction_smiles)
    mapped_core = mapped_analysis.reaction_core
    if not mapped_analysis.valid or mapped_core is None or mapped_core.quality.status != "pass":
        return CompilationResult((), "materialized_mapping_not_verified")
    if mapped_core.center_transition_key != source_core.get("center_transition_key"):
        return CompilationResult((), "materialized_mapping_core_conflict")
    mapped_observation = mapped_analysis.to_dict().get("observation") or {}
    formed_center = _formed_center(mapped_observation)
    if formed_center is None:
        return CompilationResult((), "not_single_cx_bond_formation")
    bond_kind, center_maps = formed_center

    split = split_reaction_smiles(materialized.reaction_smiles)
    if split is None:
        return CompilationResult((), "invalid_reaction_smiles")
    reactant_smiles, product_smiles = split
    reactants = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    expected_precursors = contributing_reactants(reactant_smiles, product_smiles)
    canonical_product = canonical_smiles(product_smiles)
    if reactants is None or product is None or expected_precursors is None or canonical_product is None:
        return CompilationResult((), "participant_canonicalization_failed")
    handle_signature = _direct_handle_signature(reactants, set(center_maps))
    map_lookup = _canonical_map_lookup(reactants, product, center_maps)
    _apply_map_lookup(reactants, map_lookup)
    _apply_map_lookup(product, map_lookup)
    reference_id = _reference_id(row)
    support_identity = str(
        source_core.get("mapping_equivalence_key")
        or source_core.get("core_id")
        or row.get("reaction_id")
        or row.get("observation_id")
    )
    context = context_from_analysis(mapped_analysis, materialized.reaction_smiles)
    precedent = CoreTemplatePrecedent(
        reaction_id=str(row.get("reaction_id") or ""),
        observation_id=str(row.get("observation_id") or ""),
        reference_id=reference_id,
        support_unit_id=(
            f"reference:{reference_id}"
            if reference_id
            else f"mapping_equivalence:{support_identity}"
        ),
        core_id=str(source_core.get("core_id") or ""),
        product_smiles=canonical_product,
        precursor_smiles=expected_precursors,
        mapped_reaction_smiles=materialized.reaction_smiles,
        mapping_evidence=materialized.evidence,
        mapping_confidence=materialized.confidence,
        context=context,
    )
    templates = []
    for level in requested_levels:
        product_smarts, precursor_smarts, reaction_smarts = _compile_level(
            level=level,
            reactants=reactants,
            product=product,
        )
        if compile_smarts(product_smarts, validate=False) is None:
            continue
        if not _round_trip(reaction_smarts, product_smiles, expected_precursors):
            continue
        core_key = (
            mapped_core.shape_core_key if level == "L1" else mapped_core.typed_core_key
        )
        operator_id = digest(
            "CRO1",
            level,
            bond_kind,
            core_key,
            handle_signature,
        )
        template_id = digest("CRT1", operator_id, reaction_smarts)
        templates.append(
            CoreTemplate(
                template_id=template_id,
                operator_id=operator_id,
                abstraction_level=level,
                bond_kind=bond_kind,
                reaction_smarts=reaction_smarts,
                product_smarts=product_smarts,
                precursor_smarts=precursor_smarts,
                handle_signature=handle_signature,
                shape_core_key=mapped_core.shape_core_key,
                typed_core_key=mapped_core.typed_core_key,
                center_transition_key=mapped_core.center_transition_key,
                edit_tokens=tuple(mapped_core.edit_tokens),
                observation_support=1,
                independent_reference_support=1,
                operator_observation_support=1,
                operator_reference_support=1,
                precedents=(precedent,),
            )
        )
    if not templates:
        return CompilationResult((), "source_round_trip_failed")
    return CompilationResult(tuple(templates), None)


__all__ = ["CompilationResult", "compile_core_templates"]
