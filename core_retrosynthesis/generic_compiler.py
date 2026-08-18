"""Compile diverse reaction-core and RDChiral retrosynthetic templates."""

from __future__ import annotations

import io
import json
import re
from contextlib import redirect_stdout
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, Iterable, Literal, Optional, Tuple

from rdkit import Chem
from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

from reactive_taxonomy import featurize_reaction
from .chemistry import (
    canonical_smiles,
    contributing_reactants,
    digest,
    split_reaction_smiles,
)
from .mapping import materialize_atom_mapping

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

_SMARTS_ATOM_MAP = re.compile(r":(\d+)\]")


@dataclass(frozen=True)
class GenericCompilationResult:
    """Generic templates or one deterministic rejection reason."""

    templates: Tuple[GenericCoreTemplate, ...]
    rejection_reason: Optional[str]
    rejection_stage: str = ""
    diagnostics: Dict[str, Any] | None = None


@dataclass(frozen=True)
class GenericReactionIdentity:
    """Handle-independent and realization-level identities for one reaction."""

    named_annotation: Optional[str]
    disconnection_site_key: str
    operator_signature: str
    synthon_signature: str


@dataclass(frozen=True)
class GenericOperatorIdentity:
    """Stable operator identity derived from an observed mapped signature."""

    operator_id: str
    operator_signature: str


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


_REJECTION_STAGES = {
    "missing_source_reaction": "source",
    "atom_mapping_unavailable": "mapping",
    "materialized_analysis_invalid": "observation",
    "materialized_core_missing": "core",
    "materialized_core_not_verified": "core",
    "product_completeness_not_verified": "completeness",
    "not_single_product": "completeness",
    "invalid_reaction_smiles": "canonicalization",
    "participant_canonicalization_failed": "canonicalization",
    "unsupported_edit_archetype": "operator",
    "missing_mapped_active_atoms": "operator",
    "missing_generic_operator_signature": "operator",
    "generic_l0_compilation_failed": "template",
    "template_extraction_failed": "template",
    "source_round_trip_failed": "round_trip",
}


def generic_rejection_stage(reason: str | None) -> str:
    """Return the deterministic pipeline stage for one rejection reason."""

    return _REJECTION_STAGES.get(str(reason or ""), "unknown")


def _rejection(
    reason: str,
    diagnostics: Dict[str, Any],
) -> GenericCompilationResult:
    return GenericCompilationResult(
        (),
        reason,
        generic_rejection_stage(reason),
        diagnostics,
    )


def _normalize_reaction_atom_maps(reaction_smarts: str) -> str:
    """Relabel local SMARTS maps by deterministic product-first occurrence."""

    if reaction_smarts.count(">>") != 1:
        return reaction_smarts
    product_smarts, precursor_smarts = reaction_smarts.split(">>")
    ordered = []
    for side in (product_smarts, precursor_smarts):
        for match in _SMARTS_ATOM_MAP.finditer(side):
            value = match.group(1)
            if value not in ordered:
                ordered.append(value)
    lookup = {value: index for index, value in enumerate(ordered, start=1)}

    def replace_map(match: re.Match[str]) -> str:
        return f":{lookup[match.group(1)]}]"

    return _SMARTS_ATOM_MAP.sub(replace_map, reaction_smarts)


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


def _generic_operator_signature(
    observation: Dict[str, Any],
    product: Any,
) -> str:
    """Describe mapped product edits without precursor-only handle atoms."""

    product_maps = {
        int(atom.GetAtomMapNum())
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    tokens = []
    for edit in observation.get("edits") or ():
        edit_type = str(edit.get("edit_type") or "")
        if edit_type not in {
            "formed",
            "broken",
            "order_changed",
            "hydrogen_change",
        }:
            continue
        endpoints = []
        valid = True
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field)
            if atom is None:
                endpoints.append(("H", 0, False, "S"))
                continue
            map_number = int(atom.get("atom_map_number") or 0)
            if map_number not in product_maps:
                valid = False
                break
            endpoints.append(
                (
                    str(atom.get("element") or ""),
                    int(atom.get("formal_charge") or 0),
                    bool(atom.get("aromatic")),
                    str(atom.get("hybridization") or ""),
                )
            )
        if not valid or not endpoints:
            continue
        tokens.append(
            (
                edit_type,
                tuple(sorted(endpoints)),
                str(edit.get("old_order") or "NONE"),
                str(edit.get("new_order") or "NONE"),
            )
        )
    if not tokens:
        return ""
    return json.dumps(sorted(tokens), separators=(",", ":"))


def _synthon_signature(
    reactants: Any,
    product: Any,
    operator_signature: str,
) -> str:
    """Normalize precursor-only handles while retaining mapped skeletons."""

    product_maps = {
        int(atom.GetAtomMapNum())
        for atom in product.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }
    values = []
    for fragment in Chem.GetMolFrags(reactants, asMols=True, sanitizeFrags=False):
        retained = {
            int(atom.GetIdx())
            for atom in fragment.GetAtoms()
            if int(atom.GetAtomMapNum()) in product_maps
        }
        if not retained:
            continue
        editable = Chem.RWMol(fragment)
        for atom_index in sorted(
            set(range(fragment.GetNumAtoms())).difference(retained),
            reverse=True,
        ):
            editable.RemoveAtom(atom_index)
        skeleton = editable.GetMol()
        for atom in skeleton.GetAtoms():
            atom.SetAtomMapNum(0)
        try:
            Chem.SanitizeMol(skeleton)
            values.append(
                Chem.MolToSmiles(
                    skeleton,
                    canonical=True,
                    isomericSmiles=True,
                )
            )
        except Exception:
            return ""
    if not values:
        return ""
    return digest("SYN1", operator_signature, ".".join(sorted(values)))


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


def _generalized_product_atoms(
    molecule: Any,
    active_maps: set[int],
    hydrogen_maps: set[int],
) -> Any:
    """Generalize active product atoms for the broad L0 applicability query."""

    editable = Chem.RWMol(molecule)
    for atom in molecule.GetAtoms():
        map_number = int(atom.GetAtomMapNum())
        if map_number not in active_maps or map_number in hydrogen_maps:
            continue
        if atom.GetIsAromatic():
            symbol = atom.GetSymbol().lower()
            atom_query = symbol
        else:
            atom_query = f"#{int(atom.GetAtomicNum())}"
        charge = int(atom.GetFormalCharge())
        charge_query = "+" if charge == 1 else "-" if charge == -1 else ""
        if abs(charge) > 1:
            charge_query = f"{charge:+d}"
        query = Chem.AtomFromSmarts(
            f"[{atom_query}{charge_query}:{map_number}]"
        )
        if query is None:
            raise ValueError("generic active-atom query could not be compiled")
        editable.ReplaceAtom(int(atom.GetIdx()), query, preserveProps=False)
    return editable.GetMol()


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
    levels: Iterable[Literal["L0", "L1", "L2"]] = ("L1", "L2"),
    admission_mode: Literal["supported", "data_driven"] = "supported",
) -> GenericCompilationResult:
    """Compile source-round-tripped templates from normalized graph edits."""

    requested_levels = tuple(dict.fromkeys(levels))
    if not requested_levels or any(
        level not in {"L0", "L1", "L2"} for level in requested_levels
    ):
        raise ValueError("levels must contain L0, L1, and/or L2")
    if admission_mode not in {"supported", "data_driven"}:
        raise ValueError("unsupported generic admission mode")

    source_core = dict(row.get("reaction_core") or {})
    source_observation = _observation(row)
    source_reaction = _source_reaction(row, source_observation)
    diagnostics: Dict[str, Any] = {
        "stored_core_present": bool(source_core),
        "stored_observation_present": bool(source_observation),
        "stored_core_status": str(
            (source_core.get("quality") or {}).get("status") or "missing"
        ),
        "stored_completeness_status": str(
            (
                source_observation.get("completeness")
                or row.get("reaction_completeness")
                or {}
            ).get("status")
            or "missing"
        ),
        "recomputed_evidence": True,
    }
    if not source_reaction:
        return _rejection("missing_source_reaction", diagnostics)
    prepared = _materialized_analysis(source_reaction)
    if prepared is None:
        return _rejection("atom_mapping_unavailable", diagnostics)
    materialized, analysis = prepared
    diagnostics["mapping_evidence"] = str(materialized.evidence)
    diagnostics["mapping_confidence"] = float(materialized.confidence)
    if not analysis.valid:
        return _rejection("materialized_analysis_invalid", diagnostics)
    if analysis.reaction_core is None:
        return _rejection("materialized_core_missing", diagnostics)
    if analysis.reaction_core.quality.status != "pass":
        diagnostics["recomputed_core_status"] = str(
            analysis.reaction_core.quality.status
        )
        return _rejection("materialized_core_not_verified", diagnostics)
    value = analysis.to_dict()
    observation = value.get("observation") or {}
    diagnostics["recomputed_core_status"] = "pass"
    diagnostics["recomputed_completeness_status"] = str(
        (observation.get("completeness") or {}).get("status") or "missing"
    )
    diagnostics["recovered_stored_inputs"] = bool(
        not source_core or not source_observation
    )
    if diagnostics["recomputed_completeness_status"] != "verified":
        return _rejection("product_completeness_not_verified", diagnostics)
    if len(observation.get("products") or ()) != 1:
        return _rejection("not_single_product", diagnostics)
    split = split_reaction_smiles(materialized.reaction_smiles)
    if split is None:
        return _rejection("invalid_reaction_smiles", diagnostics)
    reactant_smiles, product_smiles = split
    reactants = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    expected = contributing_reactants(reactant_smiles, product_smiles)
    canonical_product = canonical_smiles(product_smiles)
    if reactants is None or product is None or expected is None or canonical_product is None:
        return _rejection("participant_canonicalization_failed", diagnostics)
    transformation = _classify_transformation(observation, reactants, product)
    if (
        admission_mode == "supported"
        and transformation not in SUPPORTED_TRANSFORMATIONS
    ):
        return _rejection("unsupported_edit_archetype", diagnostics)
    active_maps, hydrogen_maps = _active_maps(observation)
    if not active_maps:
        return _rejection("missing_mapped_active_atoms", diagnostics)
    operator_signature = _generic_operator_signature(observation, product)
    if not operator_signature:
        return _rejection("missing_generic_operator_signature", diagnostics)
    synthon_signature = _synthon_signature(
        reactants,
        product,
        operator_signature,
    )
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
    lookup = _canonical_map_lookup(reactants, product, active_maps)
    canonical_active = {lookup[value] for value in active_maps if value in lookup}
    canonical_hydrogen = {
        lookup[value] for value in hydrogen_maps if value in lookup
    }
    handle_signature = _handle_signature(reactants, active_maps)
    _apply_map_lookup(reactants, lookup)
    _apply_map_lookup(product, lookup)
    try:
        generalized_product = _generalized_product_atoms(
            product,
            canonical_active,
            set(),
        )
        product_l0_smarts = _fragment_smarts(
            generalized_product,
            _selected_atoms(
                generalized_product,
                level="L1",
                reactant_side=False,
                active_maps=canonical_active,
            ),
            hydrogen_carrier_maps=set(),
        )
        precursor_l0_smarts = _fragment_smarts(
            reactants,
            _selected_atoms(
                reactants,
                level="L1",
                reactant_side=True,
                active_maps=canonical_active,
            ),
            hydrogen_carrier_maps=set(),
        )
    except (RuntimeError, ValueError):
        return _rejection("generic_l0_compilation_failed", diagnostics)
    normalized_l0 = _normalize_reaction_atom_maps(
        f"{product_l0_smarts}>>{precursor_l0_smarts}"
    )
    product_l0_smarts, precursor_l0_smarts = normalized_l0.split(">>")
    operator_id = digest("OP1", operator_signature)
    realization_id = digest(
        "REAL2",
        operator_id,
        handle_signature,
        precursor_l0_smarts,
    )
    annotations = (transformation,) if transformation is not None else ()
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
            return _rejection("template_extraction_failed", diagnostics)
        reaction_smarts = _normalize_reaction_atom_maps(
            str(raw["reaction_smarts"])
        )
        normalized_product_smarts, normalized_precursor_smarts = (
            reaction_smarts.split(">>")
        )
        policy = _round_trip_policy(reaction_smarts, product_smiles, expected)
        if policy is None:
            return _rejection("source_round_trip_failed", diagnostics)
        template_id = digest("GRT3", realization_id, "RDCHIRAL", reaction_smarts)
        return GenericCompilationResult(
            (
                GenericCoreTemplate(
                    template_id=template_id,
                    operator_id=operator_id,
                    transformation_kind=transformation,
                    abstraction_level="RDCHIRAL",
                    compiler_engine="rdchiral",
                    reaction_smarts=reaction_smarts,
                    product_smarts=normalized_product_smarts,
                    precursor_smarts=normalized_precursor_smarts,
                    edit_tokens=edit_tokens,
                    handle_signature=_handle_signature(reactants, active_maps),
                    stereo_policy=policy,
                    observation_support=1,
                    independent_reference_support=1,
                    precedents=(precedent,),
                    realization_id=realization_id,
                    operator_signature=operator_signature,
                    synthon_signature=synthon_signature,
                    named_annotations=annotations,
                ),
            ),
            None,
            "accepted",
            diagnostics,
        )

    templates = []
    for level in requested_levels:
        try:
            if level == "L0":
                product_smarts = product_l0_smarts
                precursor_smarts = precursor_l0_smarts
            else:
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
        reaction_smarts = _normalize_reaction_atom_maps(
            f"{product_smarts}>>{precursor_smarts}"
        )
        product_smarts, precursor_smarts = reaction_smarts.split(">>")
        policy = _round_trip_policy(reaction_smarts, product_smiles, expected)
        if policy is None:
            continue
        template_id = digest("GRT3", realization_id, level, reaction_smarts)
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
                realization_id=realization_id,
                operator_signature=operator_signature,
                synthon_signature=synthon_signature,
                named_annotations=annotations,
            )
        )
    if not templates:
        return _rejection("source_round_trip_failed", diagnostics)
    return GenericCompilationResult(
        tuple(templates),
        None,
        "accepted",
        diagnostics,
    )


@lru_cache(maxsize=50_000)
def analyze_generic_reaction(
    reaction_smiles: str,
) -> GenericReactionIdentity | None:
    """Return data-derived identities for one structurally valid reaction."""

    prepared = _materialized_analysis(reaction_smiles)
    if prepared is None:
        return None
    materialized, analysis = prepared
    return build_generic_reaction_identity(materialized.reaction_smiles, analysis)


def build_generic_reaction_identity(
    materialized_reaction_smiles: str,
    analysis: Any,
) -> GenericReactionIdentity | None:
    """Build generic action identities from an existing mapped graph analysis."""

    split = split_reaction_smiles(materialized_reaction_smiles)
    if not analysis.valid or split is None:
        return None
    reactant_smiles, product_smiles = split
    reactants = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)
    observation = analysis.to_dict().get("observation") or {}
    if reactants is None or product is None:
        return None
    signature = _generic_operator_signature(observation, product)
    if not signature:
        return None
    return GenericReactionIdentity(
        named_annotation=_classify_transformation(observation, reactants, product),
        disconnection_site_key=_product_disconnection_site_key(
            observation,
            product,
        ),
        operator_signature=signature,
        synthon_signature=_synthon_signature(reactants, product, signature),
    )


def generic_operator_identity_from_observation(
    reaction_smiles: str,
    observation: Dict[str, Any],
) -> GenericOperatorIdentity | None:
    """Return the OP1 identity already encoded by a mapped observation.

    Route-core records already contain validated graph edits.  Reusing those
    edits avoids remapping or reinterpreting the source reaction while keeping
    operator identity identical to the generic-template compiler.
    """

    if reaction_smiles.count(">>") == 1:
        _, product_smiles = reaction_smiles.split(">>", 1)
    elif reaction_smiles.count(">") == 2:
        _, _, product_smiles = reaction_smiles.split(">", 2)
    else:
        return None
    product = Chem.MolFromSmiles(product_smiles)
    if product is None:
        return None
    signature = _generic_operator_signature(observation, product)
    if not signature:
        return None
    return GenericOperatorIdentity(
        operator_id=digest("OP1", signature),
        operator_signature=signature,
    )


def classify_reaction_with_site(reaction_smiles: str) -> tuple[str | None, str]:
    """Return optional named archetype and product disconnection-site key."""

    identity = analyze_generic_reaction(reaction_smiles)
    if identity is None:
        return None, ""
    return identity.named_annotation, identity.disconnection_site_key


def classify_reaction_smiles(reaction_smiles: str) -> str | None:
    """Return the supported structural archetype for a proposed reaction."""

    return classify_reaction_with_site(reaction_smiles)[0]


__all__ = [
    "GenericCompilationResult",
    "GenericOperatorIdentity",
    "GenericReactionIdentity",
    "SUPPORTED_TRANSFORMATIONS",
    "analyze_generic_reaction",
    "build_generic_reaction_identity",
    "classify_reaction_smiles",
    "classify_reaction_with_site",
    "compile_generic_templates",
    "generic_operator_identity_from_observation",
    "generic_rejection_stage",
]
