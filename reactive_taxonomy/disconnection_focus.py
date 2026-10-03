"""Check a selected target bond using supplied, structure-validated correspondence."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .fragment_search import indexed_product
from .reaction_edits import normalize_mapped_edits
from .reaction_models import ReactionEdit, ReactionObservation
from .reaction_parser import parse_reaction_smiles


def disconnection_focus_policy() -> dict[str, Any]:
    """Load the versioned graph constraint; arbitrary executable rules are rejected."""
    path = Path(__file__).with_name("definitions") / "disconnection_bond_focus.v1.json"
    value = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "schema_version": "disconnection_bond_focus_policy.v1",
        "definition_version": "disconnection_bond_focus.v1@1.0",
        "atom_id_scope": "zero_based_canonical_smiles",
        "required_forward_edit": "formed", "require_both_precursor_atoms": True,
        "diagnostic_example_limit": 3,
    }
    if value != expected:
        raise ValueError("Unsupported disconnection-bond focus definition")
    return value


@dataclass(frozen=True)
class DisconnectionBondFocus:
    """A bond anchored to an inspected canonical target, not to input SMILES order."""

    canonical_target_smiles: str
    target_id: str
    atom_ids: tuple[int, int]
    bond_order: str
    definition_version: str
    atom_id_scope: str = "zero_based_canonical_smiles"
    schema_version: str = "disconnection_bond_focus.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the selection and its target identity."""
        return asdict(self)


@dataclass(frozen=True)
class DisconnectionBondCheck:
    """A preliminary or confirmed selected-bond witness, including failures."""

    status: str
    reason: str
    target_atom_ids: tuple[int, int]
    target_atom_maps: tuple[int, ...]
    reaction_smiles: str
    formed_edit: ReactionEdit | None = None
    warnings: tuple[str, ...] = ()
    schema_version: str = "disconnection_bond_check.v1"


def prepare_disconnection_focus(
    target_smiles: str, required_bond: tuple[int, int] | list[int],
    focus_target_smiles: str | None,
) -> DisconnectionBondFocus:
    """Validate two canonical atom IDs against the explicitly supplied inspected target."""
    policy = disconnection_focus_policy()
    if not isinstance(focus_target_smiles, str) or not focus_target_smiles.strip():
        raise ValueError("Focused search requires focus_target_smiles from the inspected molecule")
    canonical, _ = indexed_product(target_smiles)
    reference, _ = indexed_product(focus_target_smiles)
    molecule = Chem.MolFromSmiles(canonical)
    if canonical != reference:
        raise ValueError("focus_target_smiles does not identify the requested target")
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Focused target must be one connected molecule")
    if (not isinstance(required_bond, (tuple, list)) or len(required_bond) != 2
            or any(type(index) is not int for index in required_bond)):
        raise ValueError("required_disconnection_bond must contain two canonical integer atom IDs")
    first, second = sorted(required_bond)
    if first < 0 or second >= molecule.GetNumAtoms() or first == second:
        raise ValueError("Selected canonical atom IDs must be distinct and in range")
    bond = molecule.GetBondBetweenAtoms(first, second)
    if bond is None:
        raise ValueError("Selected canonical atoms are not bonded")
    identity = hashlib.sha256(json.dumps(
        [canonical, policy["definition_version"]], separators=(",", ":"),
    ).encode()).hexdigest()
    return DisconnectionBondFocus(
        canonical, f"FOCUS_TARGET1:{identity}", (first, second),
        str(bond.GetBondType()), policy["definition_version"],
    )


def check_disconnection_bond(
    focus: DisconnectionBondFocus, reaction_smiles: str,
    canonical_target_atom_maps: tuple[int, ...], *,
    observation: ReactionObservation | None = None,
) -> DisconnectionBondCheck:
    """Check a forward formed bond without inventing correspondence.

    The caller supplies the generator's actual map numbers in canonical target
    atom order. Reconstructing that labelled target verifies the binding against
    the proposed product, including stereo. An optional final observation must
    explicitly validate the same supplied mapping without conflicting evidence.
    """
    selected_maps: tuple[int, ...] = ()

    def result(status: str, reason: str, edit: ReactionEdit | None = None,
               warnings: tuple[str, ...] = ()) -> DisconnectionBondCheck:
        return DisconnectionBondCheck(
            status, reason, focus.atom_ids, selected_maps, reaction_smiles, edit, warnings,
        )

    target = Chem.MolFromSmiles(focus.canonical_target_smiles)
    if (len(canonical_target_atom_maps) != target.GetNumAtoms()
            or len(set(canonical_target_atom_maps)) != target.GetNumAtoms()
            or any(type(number) is not int or number <= 0 for number in canonical_target_atom_maps)):
        return result("unresolved", "invalid_generator_target_binding")
    for atom, number in zip(target.GetAtoms(), canonical_target_atom_maps):
        atom.SetAtomMapNum(number)
    selected_maps = tuple(canonical_target_atom_maps[index] for index in focus.atom_ids)
    parsed = parse_reaction_smiles(reaction_smiles, include_molecular_interpretation=False)
    if not parsed.valid or len(parsed.products) != 1:
        return result("unresolved", "invalid_or_multiple_products")
    product = Chem.MolFromSmiles(parsed.products[0].input_smiles)
    if product is None or Chem.MolToSmiles(product) != Chem.MolToSmiles(target):
        return result("unresolved", "product_or_atom_binding_mismatch")
    precursor_maps = {
        atom.GetAtomMapNum()
        for component in parsed.reactants
        for atom in Chem.MolFromSmiles(component.input_smiles).GetAtoms()
    }
    if not set(selected_maps) <= precursor_maps:
        return result("unresolved", "selected_atom_missing_from_precursors")
    if observation is None:
        normalized = normalize_mapped_edits(parsed.reactants, parsed.products)
        if not normalized.valid or normalized.edit_hypotheses:
            return result("unresolved", "invalid_or_ambiguous_mapping", warnings=normalized.warnings)
        edits = normalized.edits
    else:
        if (not observation.valid or observation.input_reaction_smiles != reaction_smiles
                or observation.evidence_quality != "validated_atom_mapping"
                or observation.edit_hypotheses
                or any(candidate.status in {"ambiguous", "invalid"}
                       for candidate in observation.evidence_candidates)
                or any("CONFLICT" in warning.upper() for warning in observation.warnings)):
            return result("unresolved", "mapping_not_confirmed_by_final_observation",
                          warnings=observation.warnings)
        edits = observation.edits
    for edit in edits:
        if (edit.edit_type == "formed" and edit.atom_2 is not None
                and {edit.atom_1.atom_map_number, edit.atom_2.atom_map_number} == set(selected_maps)):
            return result("verified" if observation is not None else "matched",
                          "selected_bond_formed_forward", edit)
    return result("not_matched", "selected_bond_not_formed_forward")
