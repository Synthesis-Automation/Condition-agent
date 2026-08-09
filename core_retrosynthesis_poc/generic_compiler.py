"""Compile diverse reaction-core and RDChiral retrosynthetic templates."""

from __future__ import annotations

import io
import json
from contextlib import redirect_stdout
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, Iterable, Literal, Optional, Tuple

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

from reactive_taxonomy import featurize_reaction
from retrosynthesis_poc.chemistry import (
    canonical_smiles,
    contributing_reactants,
    digest,
    split_reaction_smiles,
)
from retrosynthesis_poc.mapping import materialize_atom_mapping

from .compiler import (
    _apply_map_lookup,
    _fragment_smarts,
    _selected_atoms,
)
from .context import context_from_analysis
from .generic_models import GenericCoreTemplate, GenericTemplatePrecedent


SUPPORTED_TRANSFORMATIONS = frozenset(
    {
        "acyl_substitution",
        "c_c_coupling",
        "carbonyl_oxidation",
        "carbonyl_reduction",
        "conjugate_addition",
        "carbonyl_condensation",
        "ring_formation",
    }
)


@dataclass(frozen=True)
class GenericCompilationResult:
    """Generic templates or one deterministic rejection reason."""

    templates: Tuple[GenericCoreTemplate, ...]
    rejection_reason: Optional[str]


@lru_cache(maxsize=20_000)
def _materialized_analysis(reaction_smiles: str) -> tuple[Any, Any] | None:
    """Cache deterministic atom mapping and featurization across compilers."""

    materialized = materialize_atom_mapping(reaction_smiles)
    if materialized is None:
        return None
    return materialized, featurize_reaction(materialized.reaction_smiles)


def _observation(row: Dict[str, Any]) -> Dict[str, Any]:
    return dict(row.get("reaction_observation") or row.get("observation") or {})


def _source_reaction(row: Dict[str, Any], observation: Dict[str, Any]) -> str:
    return str(
        row.get("reaction_smiles")
        or observation.get("input_reaction_smiles")
        or row.get("canonical_reaction_smiles")
        or ""
    )


def _reference_id(row: Dict[str, Any]) -> str:
    direct = str(row.get("reference_id") or "").strip()
    if direct:
        return direct
    return str((row.get("reference_identity") or {}).get("reference_id") or "")


def _atom_by_map(molecule: Any, map_number: int) -> Any | None:
    for atom in molecule.GetAtoms():
        if int(atom.GetAtomMapNum()) == map_number:
            return atom
    return None


def _has_carbonyl(molecule: Any, map_number: int) -> bool:
    atom = _atom_by_map(molecule, map_number)
    if atom is None or atom.GetSymbol() != "C":
        return False
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetOtherAtom(atom).GetSymbol() == "O"
        for bond in atom.GetBonds()
    )


def _has_direct_handle(molecule: Any, map_number: int, elements: set[str]) -> bool:
    atom = _atom_by_map(molecule, map_number)
    if atom is None:
        return False
    return any(
        int(neighbor.GetAtomMapNum()) == 0 and neighbor.GetSymbol() in elements
        for neighbor in atom.GetNeighbors()
    )


def _formed_edits(observation: Dict[str, Any]) -> tuple[Dict[str, Any], ...]:
    return tuple(
        edit
        for edit in observation.get("edits") or ()
        if edit.get("edit_type") == "formed" and edit.get("atom_2") is not None
    )


def _classify_transformation(
    observation: Dict[str, Any],
    reactants: Any,
    product: Any,
) -> str | None:
    edits = tuple(observation.get("edits") or ())
    formed = _formed_edits(observation)
    order_tokens = {
        (
            frozenset(
                {
                    str((edit.get("atom_1") or {}).get("element") or ""),
                    str((edit.get("atom_2") or {}).get("element") or ""),
                }
            ),
            str(edit.get("old_order") or "").upper(),
            str(edit.get("new_order") or "").upper(),
        )
        for edit in edits
        if edit.get("edit_type") == "order_changed"
    }
    if (
        len(formed) == 2
        and all(
            {
                str((edit.get("atom_1") or {}).get("element") or ""),
                str((edit.get("atom_2") or {}).get("element") or ""),
            }
            == {"C", "N"}
            for edit in formed
        )
        and any(
            elements == frozenset({"C"}) and old == "TRIPLE"
            for elements, old, _ in order_tokens
        )
    ):
        return "ring_formation"
    if not formed and len(order_tokens) == 1:
        elements, old, new = next(iter(order_tokens))
        if elements == frozenset({"C", "O"}):
            if old == "SINGLE" and new == "DOUBLE":
                return "carbonyl_oxidation"
            if old == "DOUBLE" and new == "SINGLE":
                return "carbonyl_reduction"
    if len(formed) != 1:
        return None
    endpoints = (formed[0].get("atom_1") or {}, formed[0].get("atom_2") or {})
    elements = {str(atom.get("element") or "") for atom in endpoints}
    maps = {
        str(atom.get("element") or ""): int(atom.get("atom_map_number") or 0)
        for atom in endpoints
    }
    if elements == {"C"} and any(
        token_elements == frozenset({"C"})
        and old == "DOUBLE"
        and new == "SINGLE"
        for token_elements, old, new in order_tokens
    ):
        return "conjugate_addition"
    if elements == {"C"}:
        carbon_maps = [
            int(atom.get("atom_map_number") or 0)
            for atom in endpoints
            if int(atom.get("atom_map_number") or 0) > 0
        ]
        has_boron = any(
            _has_direct_handle(reactants, value, {"B"}) for value in carbon_maps
        )
        has_leaving_group = any(
            _has_direct_handle(reactants, value, {"Cl", "Br", "I"})
            for value in carbon_maps
        )
        return "c_c_coupling" if has_boron and has_leaving_group else None
    if elements == {"C", "N"} and maps.get("C", 0) > 0:
        carbon_map = maps["C"]
        before = _has_carbonyl(reactants, carbon_map)
        after = _has_carbonyl(product, carbon_map)
        if before and after:
            return "acyl_substitution"
        if before and not after:
            return "carbonyl_condensation"
    return None


def _product_disconnection_site_key(
    observation: Dict[str, Any],
    product: Any,
) -> str:
    """Identify edited product bonds independently of precursor handle form."""

    unnumbered = Chem.Mol(product)
    for atom in unnumbered.GetAtoms():
        atom.SetAtomMapNum(0)
    ranks = tuple(
        int(rank)
        for rank in Chem.CanonicalRankAtoms(
            unnumbered,
            breakTies=False,
            includeChirality=True,
        )
    )
    rank_by_map = {
        int(atom.GetAtomMapNum()): ranks[int(atom.GetIdx())]
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    descriptors = []
    for edit in observation.get("edits") or ():
        edit_type = str(edit.get("edit_type") or "")
        if edit_type not in {"formed", "order_changed"}:
            continue
        endpoints = []
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field) or {}
            map_number = int(atom.get("atom_map_number") or 0)
            if map_number not in rank_by_map:
                endpoints = []
                break
            endpoints.append(
                (rank_by_map[map_number], str(atom.get("element") or ""))
            )
        if len(endpoints) != 2:
            continue
        descriptors.append(
            {
                "edit_type": edit_type,
                "product_atoms": sorted(endpoints),
                "product_order": str(edit.get("new_order") or ""),
            }
        )
    if not descriptors:
        return ""
    descriptors.sort(
        key=lambda descriptor: json.dumps(
            descriptor,
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    canonical_product = Chem.MolToSmiles(
        unnumbered,
        canonical=True,
        isomericSmiles=True,
    )
    return digest(
        "SITE1",
        canonical_product,
        json.dumps(descriptors, sort_keys=True, separators=(",", ":")),
    )


def _active_maps(observation: Dict[str, Any]) -> tuple[set[int], set[int]]:
    active = set()
    hydrogen = set()
    for edit in observation.get("edits") or ():
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field) or {}
            map_number = int(atom.get("atom_map_number") or 0)
            if map_number > 0:
                active.add(map_number)
        if edit.get("edit_type") == "hydrogen_change":
            map_number = int((edit.get("atom_1") or {}).get("atom_map_number") or 0)
            if map_number > 0:
                hydrogen.add(map_number)
    return active, hydrogen


def _canonical_map_lookup(
    reactants: Any,
    product: Any,
    active_maps: set[int],
) -> dict[int, int]:
    product_without_maps = Chem.Mol(product)
    for atom in product_without_maps.GetAtoms():
        atom.SetAtomMapNum(0)
    ranks = tuple(
        int(rank)
        for rank in Chem.CanonicalRankAtoms(
            product_without_maps,
            breakTies=True,
            includeChirality=True,
        )
    )
    product_rank = {
        int(atom.GetAtomMapNum()): ranks[int(atom.GetIdx())]
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    atoms = {
        int(atom.GetAtomMapNum()): atom
        for molecule in (reactants, product)
        for atom in molecule.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }

    def key(map_number: int) -> tuple[Any, ...]:
        atom = atoms[map_number]
        return (
            product_rank.get(map_number, 10**9),
            atom.GetSymbol(),
            int(atom.GetIsAromatic()),
            int(atom.GetFormalCharge()),
            map_number,
        )

    ordered_active = sorted(active_maps, key=key)
    remaining = sorted(set(atoms).difference(active_maps), key=key)
    return {
        map_number: canonical
        for canonical, map_number in enumerate(
            (*ordered_active, *remaining),
            start=1,
        )
    }


def _handle_signature(reactants: Any, active_maps: set[int]) -> str:
    values = []
    for atom in reactants.GetAtoms():
        if int(atom.GetAtomMapNum()) not in active_maps:
            continue
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if int(neighbor.GetAtomMapNum()) == 0:
                values.append(
                    f"{atom.GetSymbol()}-{neighbor.GetSymbol()}:{bond.GetBondType().name}"
                )
    return ";".join(sorted(values)) or "no_unmapped_handle"


def _without_stereo(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)


def _round_trip_policy(
    reaction_smarts: str,
    product_smiles: str,
    expected_precursors: str,
) -> Literal["exact", "relaxed"] | None:
    try:
        with redirect_stdout(io.StringIO()):
            outcomes = rdchiralRun(
                rdchiralReaction(reaction_smarts),
                rdchiralReactants(product_smiles),
            )
    except Exception:
        return None
    canonical_outcomes = {
        canonical
        for outcome in outcomes
        if (canonical := canonical_smiles(str(outcome))) is not None
    }
    if expected_precursors in canonical_outcomes:
        return "exact"
    expected_relaxed = _without_stereo(expected_precursors)
    if expected_relaxed and expected_relaxed in {
        relaxed
        for outcome in canonical_outcomes
        if (relaxed := _without_stereo(outcome)) is not None
    }:
        return "relaxed"
    return None


def compile_generic_templates(
    row: Dict[str, Any],
    *,
    engine: Literal["reaction_core", "rdchiral"] = "reaction_core",
    levels: Iterable[Literal["L1", "L2"]] = ("L1", "L2"),
) -> GenericCompilationResult:
    """Compile source-round-tripped templates for supported edit archetypes."""

    source_core = dict(row.get("reaction_core") or {})
    source_observation = _observation(row)
    if not source_core or not source_observation:
        return GenericCompilationResult((), "missing_core_or_observation")
    if (source_core.get("quality") or {}).get("status") not in {"pass", "review"}:
        return GenericCompilationResult((), "core_quality_blocked")
    completeness = source_observation.get("completeness") or row.get(
        "reaction_completeness"
    ) or {}
    if completeness.get("status") != "verified":
        return GenericCompilationResult((), "product_completeness_not_verified")
    if len(source_observation.get("products") or ()) != 1:
        return GenericCompilationResult((), "not_single_product")
    prepared = _materialized_analysis(_source_reaction(row, source_observation))
    if prepared is None:
        return GenericCompilationResult((), "atom_mapping_unavailable")
    materialized, analysis = prepared
    if (
        not analysis.valid
        or analysis.reaction_core is None
        or analysis.reaction_core.quality.status != "pass"
    ):
        return GenericCompilationResult((), "materialized_mapping_not_verified")
    value = analysis.to_dict()
    observation = value.get("observation") or {}
    split = split_reaction_smiles(materialized.reaction_smiles)
    if split is None:
        return GenericCompilationResult((), "invalid_reaction_smiles")
    reactant_smiles, product_smiles = split
    reactants = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    expected = contributing_reactants(reactant_smiles, product_smiles)
    canonical_product = canonical_smiles(product_smiles)
    if reactants is None or product is None or expected is None or canonical_product is None:
        return GenericCompilationResult((), "participant_canonicalization_failed")
    transformation = _classify_transformation(observation, reactants, product)
    if transformation not in SUPPORTED_TRANSFORMATIONS:
        return GenericCompilationResult((), "unsupported_edit_archetype")
    active_maps, hydrogen_maps = _active_maps(observation)
    if not active_maps:
        return GenericCompilationResult((), "missing_mapped_active_atoms")
    context = context_from_analysis(analysis, materialized.reaction_smiles)
    precedent = GenericTemplatePrecedent(
        reaction_id=str(row.get("reaction_id") or ""),
        reference_id=_reference_id(row),
        product_smiles=canonical_product,
        precursor_smiles=expected,
        mapped_reaction_smiles=materialized.reaction_smiles,
        context=context,
    )
    edit_tokens = tuple(analysis.reaction_core.edit_tokens)
    if engine == "rdchiral":
        try:
            with redirect_stdout(io.StringIO()):
                raw = extract_from_reaction(
                    {
                        "_id": "generic_poc",
                        "reactants": reactant_smiles,
                        "products": product_smiles,
                    }
                )
        except Exception:
            raw = None
        if not raw or not raw.get("reaction_smarts"):
            return GenericCompilationResult((), "template_extraction_failed")
        reaction_smarts = str(raw["reaction_smarts"])
        policy = _round_trip_policy(reaction_smarts, product_smiles, expected)
        if policy is None:
            return GenericCompilationResult((), "source_round_trip_failed")
        operator_id = digest(
            "GRO1",
            "rdchiral",
            transformation,
            analysis.reaction_core.typed_core_key,
        )
        template_id = digest("GRT1", operator_id, reaction_smarts)
        return GenericCompilationResult(
            (
                GenericCoreTemplate(
                    template_id=template_id,
                    operator_id=operator_id,
                    transformation_kind=transformation,
                    abstraction_level="RDCHIRAL",
                    compiler_engine="rdchiral",
                    reaction_smarts=reaction_smarts,
                    product_smarts=str(raw["products"]),
                    precursor_smarts=str(raw["reactants"]),
                    edit_tokens=edit_tokens,
                    handle_signature=_handle_signature(reactants, active_maps),
                    stereo_policy=policy,
                    observation_support=1,
                    independent_reference_support=1,
                    precedents=(precedent,),
                ),
            ),
            None,
        )

    lookup = _canonical_map_lookup(reactants, product, active_maps)
    canonical_active = {lookup[value] for value in active_maps if value in lookup}
    canonical_hydrogen = {
        lookup[value] for value in hydrogen_maps if value in lookup
    }
    handle_signature = _handle_signature(reactants, active_maps)
    _apply_map_lookup(reactants, lookup)
    _apply_map_lookup(product, lookup)
    templates = []
    for level in tuple(dict.fromkeys(levels)):
        try:
            product_smarts = _fragment_smarts(
                product,
                _selected_atoms(
                    product,
                    level=level,
                    reactant_side=False,
                    active_maps=canonical_active,
                ),
                hydrogen_carrier_maps=canonical_hydrogen,
            )
            precursor_smarts = _fragment_smarts(
                reactants,
                _selected_atoms(
                    reactants,
                    level=level,
                    reactant_side=True,
                    active_maps=canonical_active,
                ),
                hydrogen_carrier_maps=canonical_hydrogen,
            )
        except (RuntimeError, ValueError):
            continue
        reaction_smarts = f"{product_smarts}>>{precursor_smarts}"
        policy = _round_trip_policy(reaction_smarts, product_smiles, expected)
        if policy is None:
            continue
        core_key = (
            analysis.reaction_core.shape_core_key
            if level == "L1"
            else analysis.reaction_core.typed_core_key
        )
        operator_id = digest(
            "GRO1",
            "reaction_core",
            transformation,
            level,
            core_key,
            handle_signature,
        )
        template_id = digest("GRT1", operator_id, reaction_smarts)
        templates.append(
            GenericCoreTemplate(
                template_id=template_id,
                operator_id=operator_id,
                transformation_kind=transformation,
                abstraction_level=level,
                compiler_engine="reaction_core",
                reaction_smarts=reaction_smarts,
                product_smarts=product_smarts,
                precursor_smarts=precursor_smarts,
                edit_tokens=edit_tokens,
                handle_signature=handle_signature,
                stereo_policy=policy,
                observation_support=1,
                independent_reference_support=1,
                precedents=(precedent,),
            )
        )
    return GenericCompilationResult(
        tuple(templates),
        None if templates else "source_round_trip_failed",
    )


@lru_cache(maxsize=50_000)
def classify_reaction_with_site(reaction_smiles: str) -> tuple[str | None, str]:
    """Return graph-edit archetype and product-side disconnection-site key."""

    prepared = _materialized_analysis(reaction_smiles)
    if prepared is None:
        return None, ""
    materialized, analysis = prepared
    split = split_reaction_smiles(materialized.reaction_smiles)
    if not analysis.valid or split is None:
        return None, ""
    reactant_smiles, product_smiles = split
    reactants = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    observation = analysis.to_dict().get("observation") or {}
    if reactants is None or product is None:
        return None, ""
    return (
        _classify_transformation(observation, reactants, product),
        _product_disconnection_site_key(observation, product),
    )


def classify_reaction_smiles(reaction_smiles: str) -> str | None:
    """Return the supported structural archetype for a proposed reaction."""

    return classify_reaction_with_site(reaction_smiles)[0]


__all__ = [
    "GenericCompilationResult",
    "SUPPORTED_TRANSFORMATIONS",
    "classify_reaction_smiles",
    "classify_reaction_with_site",
    "compile_generic_templates",
]
