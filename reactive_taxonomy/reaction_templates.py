"""Versioned, single-event reaction-template registry.

The registry is an authoring and interpretation layer over normalized reaction
edits.  It does not execute arbitrary code and template names never replace
query-derived structural evidence.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping, Optional, Sequence, Tuple

from .chemistry.rdkit_utils import mol_to_canonical_smiles, parse_smiles
from .reaction_archetypes import infer_edit_archetype
from .reaction_edits import normalize_mapped_edits
from .reaction_events import partition_reaction_edits
from .reaction_graph_editing import bond_type, set_total_hydrogens
from .reaction_models import ReactionAtomReference, ReactionEdit
from .reaction_parser import parse_reaction_smiles


REACTION_TEMPLATE_SCHEMA_VERSION = "1.0"
REACTION_TEMPLATE_DEFINITION_VERSION = "reaction_templates.v1"
DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH = (
    Path(__file__).with_name("definitions") / "reaction_templates.v1.json"
)

_TEMPLATE_ID_RE = re.compile(r"^[a-z][a-z0-9_]{2,79}$")
_VALID_STATUS = {"draft", "active", "retired"}


class ReactionTemplateError(ValueError):
    """Raised when a template or registry violates its public contract."""


@dataclass(frozen=True)
class ReactionTemplateAtom:
    """One mapped atom participating in a reference reaction edit."""

    atom_map_number: int
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str

    @classmethod
    def from_reference(
        cls, reference: ReactionAtomReference
    ) -> "ReactionTemplateAtom":
        """Build an atom contract from a mapped reaction reference."""
        if reference.atom_map_number is None:
            raise ReactionTemplateError(
                "Template edit atoms require atom-map numbers"
            )
        return cls(
            atom_map_number=int(reference.atom_map_number),
            element=str(reference.element),
            formal_charge=int(reference.formal_charge),
            aromatic=bool(reference.aromatic),
            hybridization=str(reference.hybridization),
        )


@dataclass(frozen=True)
class ReactionTemplateEdit:
    """One normalized edit derived from the mapped reference reaction."""

    edit_type: Literal[
        "formed", "broken", "order_changed", "hydrogen_change"
    ]
    atom_1: ReactionTemplateAtom
    atom_2: Optional[ReactionTemplateAtom]
    old_order: Optional[str]
    new_order: Optional[str]

    @classmethod
    def from_reaction_edit(
        cls, edit: ReactionEdit
    ) -> "ReactionTemplateEdit":
        """Convert a normalized reaction edit into a serializable template edit."""
        return cls(
            edit_type=edit.edit_type,
            atom_1=ReactionTemplateAtom.from_reference(edit.atom_1),
            atom_2=(
                ReactionTemplateAtom.from_reference(edit.atom_2)
                if edit.atom_2 is not None
                else None
            ),
            old_order=edit.old_order,
            new_order=edit.new_order,
        )


@dataclass(frozen=True)
class ReactionTemplateParticipant:
    """One canonical species and its explicit reference multiplicity."""

    side: Literal["reactant", "agent", "product"]
    canonical_smiles: str
    explicit_count: int


@dataclass(frozen=True)
class ReactionTemplate:
    """One validated, declarative single-event reaction template."""

    template_id: str
    display_name: str
    family_id: Optional[str]
    aliases: Tuple[str, ...]
    mapped_reference_reaction: str
    participants: Tuple[ReactionTemplateParticipant, ...]
    edits: Tuple[ReactionTemplateEdit, ...]
    edit_fingerprint: str
    edit_component_count: int
    edit_archetype: str
    transformation_class: Optional[str]
    status: Literal["draft", "active", "retired"]
    provenance: str
    notes: str
    definition_hash: str
    schema_version: str = REACTION_TEMPLATE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Serialize the template with deterministic field values."""
        return asdict(self)


@dataclass(frozen=True)
class ReactionTemplateMatch:
    """One template whose edit fingerprint matches a query reaction."""

    template_id: str
    display_name: str
    family_id: Optional[str]
    status: str
    edit_fingerprint: str
    definition_hash: str
    evidence: str
    confidence: float
    provisional: bool
    predicted_product_smiles: Optional[str]
    inferred_multiplicity: bool

    def to_dict(self) -> dict[str, Any]:
        """Serialize the match."""
        return asdict(self)


@dataclass(frozen=True)
class ReactionTemplateQueryResult:
    """Query-derived signature information and matching registry templates."""

    reaction_smiles: str
    valid: bool
    evidence: str
    edit_fingerprint: Optional[str]
    signature_id: Optional[str]
    matches: Tuple[ReactionTemplateMatch, ...]
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize the query result."""
        return asdict(self)


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    )


def _digest(prefix: str, value: Any, *, length: int = 32) -> str:
    encoded = _canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:length]


def _atom_identity(atom: ReactionAtomReference) -> tuple[object, ...]:
    if atom.atom_map_number is not None:
        return ("map", int(atom.atom_map_number))
    return (
        "atom",
        str(atom.side),
        int(atom.component_index),
        int(atom.atom_index),
    )


def _atom_label(atom: ReactionAtomReference) -> str:
    return "|".join(
        (
            str(atom.element),
            str(int(atom.formal_charge)),
            "aromatic" if atom.aromatic else "aliphatic",
            str(atom.hybridization),
        )
    )


def reaction_edit_fingerprint(edits: Sequence[ReactionEdit]) -> str:
    """Return a map- and component-order-invariant edit-graph fingerprint.

    The fingerprint intentionally excludes local-environment IDs, source
    labels, template names, and atom-map numbers.  Per-atom incident edit
    multisets retain the connectivity of a compact single event more strongly
    than a flat multiset of bond changes.
    """
    if not edits:
        raise ReactionTemplateError("Cannot fingerprint an empty edit set")
    node_labels: dict[tuple[object, ...], str] = {}
    incidents: dict[tuple[object, ...], list[str]] = {}
    edit_tokens = []
    for edit in edits:
        left_key = _atom_identity(edit.atom_1)
        left_label = _atom_label(edit.atom_1)
        node_labels[left_key] = left_label
        incidents.setdefault(left_key, [])
        if edit.atom_2 is None:
            token = "|".join(
                (
                    edit.edit_type,
                    left_label,
                    "H",
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
            incidents[left_key].append(
                "|".join(
                    (
                        edit.edit_type,
                        "H",
                        edit.old_order or "NONE",
                        edit.new_order or "NONE",
                    )
                )
            )
            edit_tokens.append(token)
            continue
        right_key = _atom_identity(edit.atom_2)
        right_label = _atom_label(edit.atom_2)
        node_labels[right_key] = right_label
        incidents.setdefault(right_key, [])
        endpoint_labels = tuple(sorted((left_label, right_label)))
        token = "|".join(
            (
                edit.edit_type,
                endpoint_labels[0],
                endpoint_labels[1],
                edit.old_order or "NONE",
                edit.new_order or "NONE",
            )
        )
        edit_tokens.append(token)
        incidents[left_key].append(
            "|".join(
                (
                    edit.edit_type,
                    right_label,
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
        )
        incidents[right_key].append(
            "|".join(
                (
                    edit.edit_type,
                    left_label,
                    edit.old_order or "NONE",
                    edit.new_order or "NONE",
                )
            )
        )
    node_tokens = sorted(
        [
            {
                "label": node_labels[key],
                "incidents": sorted(values),
            }
            for key, values in incidents.items()
        ],
        key=_canonical_json,
    )
    return _digest(
        "RTE1",
        {
            "edits": sorted(edit_tokens),
            "nodes": node_tokens,
        },
    )


def _unmapped_canonical_smiles(smiles: str) -> str:
    molecule = parse_smiles(smiles)
    if molecule is None:
        raise ReactionTemplateError(f"Invalid component SMILES: {smiles}")
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = mol_to_canonical_smiles(molecule)
    if canonical is None:
        raise ReactionTemplateError(f"Cannot canonicalize component: {smiles}")
    return canonical


def _participants(parsed: Any) -> Tuple[ReactionTemplateParticipant, ...]:
    counts: dict[tuple[str, str], int] = {}
    for component in (*parsed.reactants, *parsed.agents, *parsed.products):
        canonical = _unmapped_canonical_smiles(component.input_smiles)
        key = (str(component.side), canonical)
        counts[key] = counts.get(key, 0) + 1
    return tuple(
        ReactionTemplateParticipant(
            side=side,
            canonical_smiles=canonical,
            explicit_count=count,
        )
        for (side, canonical), count in sorted(counts.items())
    )


def _neighbor_token(
    *,
    element: str,
    aromatic: bool,
    order: Optional[str],
    hydrogen_count: Optional[int] = None,
    external_heavy_count: Optional[int] = None,
) -> str:
    values = [
        str(element),
        "aromatic" if aromatic else "aliphatic",
        str(order or "NONE").upper(),
    ]
    if hydrogen_count is not None:
        values.append(f"H{int(hydrogen_count)}")
    if external_heavy_count is not None:
        values.append(f"X{int(external_heavy_count)}")
    return "|".join(values)


def _mapped_atom_states(
    components: Sequence[Any],
) -> dict[int, tuple[int, int]]:
    """Return map-number keyed H count and heavy-atom degree."""
    states: dict[int, tuple[int, int]] = {}
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            map_number = int(atom.GetAtomMapNum())
            if not map_number:
                continue
            states[map_number] = (
                int(atom.GetTotalNumHs(includeNeighbors=True)),
                sum(
                    neighbor.GetAtomicNum() > 1
                    for neighbor in atom.GetNeighbors()
                ),
            )
    return states


def _template_neighbor_token(
    *,
    atom: ReactionTemplateAtom,
    order: Optional[str],
    side_states: Mapping[int, tuple[int, int]],
) -> str:
    hydrogen_count, heavy_degree = side_states.get(
        int(atom.atom_map_number),
        (0, 1),
    )
    return _neighbor_token(
        element=atom.element,
        aromatic=atom.aromatic,
        order=order,
        hydrogen_count=hydrogen_count,
        external_heavy_count=max(0, heavy_degree - 1),
    )


@dataclass(frozen=True)
class _CenterTransitionRequirement:
    element: str
    aromatic: bool
    before_tokens: Tuple[str, ...]
    after_tokens: Tuple[str, ...]


@dataclass(frozen=True)
class _RoleAtomPattern:
    element: str
    formal_charge: int
    aromatic: bool
    minimum_hydrogen_count: int


@dataclass(frozen=True)
class _RoleBondPattern:
    atom_slot_1: int
    atom_slot_2: int
    order: str


@dataclass(frozen=True)
class _CompiledReactantRole:
    role_id: str
    required_count: int
    allow_repeated_component: bool
    atom_patterns: Tuple[_RoleAtomPattern, ...]
    bond_patterns: Tuple[_RoleBondPattern, ...]
    occurrence_map_numbers: Tuple[Tuple[int, ...], ...]


@dataclass(frozen=True)
class _RoleCandidate:
    component_index: int
    atom_indices: Tuple[int, ...]


@dataclass(frozen=True)
class _TemplateReconstruction:
    predicted_product_smiles: str
    inferred_multiplicity: bool


def _hydrogen_loss_requirements(
    edits: Sequence[ReactionTemplateEdit],
) -> dict[int, int]:
    losses: dict[int, int] = {}
    for edit in edits:
        if (
            edit.edit_type == "hydrogen_change"
            and edit.old_order is not None
            and edit.new_order is None
        ):
            map_number = int(edit.atom_1.atom_map_number)
            losses[map_number] = losses.get(map_number, 0) + 1
    return losses


def _compile_reactant_roles(
    template: ReactionTemplate,
) -> Tuple[_CompiledReactantRole, ...]:
    """Compile minimal edit-bearing reactant roles from a mapped reference."""
    parsed = parse_reaction_smiles(template.mapped_reference_reaction)
    edited_maps = {
        int(atom.atom_map_number)
        for edit in template.edits
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None
    }
    hydrogen_losses = _hydrogen_loss_requirements(template.edits)
    instances: list[
        tuple[
            str,
            Tuple[_RoleAtomPattern, ...],
            Tuple[_RoleBondPattern, ...],
            Tuple[int, ...],
            str,
        ]
    ] = []
    for component in parsed.reactants:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        mapped_atoms = {
            int(atom.GetAtomMapNum()): int(atom.GetIdx())
            for atom in molecule.GetAtoms()
            if int(atom.GetAtomMapNum()) in edited_maps
        }
        if not mapped_atoms:
            continue
        ordered_maps = tuple(
            sorted(
                mapped_atoms,
                key=lambda map_number: (
                    molecule.GetAtomWithIdx(mapped_atoms[map_number]).GetSymbol(),
                    bool(
                        molecule.GetAtomWithIdx(
                            mapped_atoms[map_number]
                        ).GetIsAromatic()
                    ),
                    int(
                        molecule.GetAtomWithIdx(
                            mapped_atoms[map_number]
                        ).GetFormalCharge()
                    ),
                    map_number,
                ),
            )
        )
        patterns = tuple(
            _RoleAtomPattern(
                element=str(
                    molecule.GetAtomWithIdx(mapped_atoms[map_number]).GetSymbol()
                ),
                formal_charge=int(
                    molecule.GetAtomWithIdx(
                        mapped_atoms[map_number]
                    ).GetFormalCharge()
                ),
                aromatic=bool(
                    molecule.GetAtomWithIdx(
                        mapped_atoms[map_number]
                    ).GetIsAromatic()
                ),
                minimum_hydrogen_count=int(
                    hydrogen_losses.get(map_number, 0)
                ),
            )
            for map_number in ordered_maps
        )
        slot_by_map = {
            map_number: slot for slot, map_number in enumerate(ordered_maps)
        }
        bonds = []
        for left_position, left_map in enumerate(ordered_maps):
            for right_map in ordered_maps[left_position + 1 :]:
                bond = molecule.GetBondBetweenAtoms(
                    mapped_atoms[left_map],
                    mapped_atoms[right_map],
                )
                if bond is not None:
                    bonds.append(
                        _RoleBondPattern(
                            atom_slot_1=slot_by_map[left_map],
                            atom_slot_2=slot_by_map[right_map],
                            order=str(bond.GetBondType()).upper(),
                        )
                    )
        pattern_key = _canonical_json(
            {
                "atoms": [asdict(pattern) for pattern in patterns],
                "bonds": [asdict(bond) for bond in bonds],
            }
        )
        instances.append(
            (
                pattern_key,
                patterns,
                tuple(bonds),
                ordered_maps,
                _unmapped_canonical_smiles(component.input_smiles),
            )
        )
    grouped: dict[
        str,
        list[
            tuple[
                Tuple[_RoleAtomPattern, ...],
                Tuple[_RoleBondPattern, ...],
                Tuple[int, ...],
                str,
            ]
        ],
    ] = {}
    for pattern_key, patterns, bonds, maps, canonical in instances:
        grouped.setdefault(pattern_key, []).append(
            (patterns, bonds, maps, canonical)
        )
    roles = []
    for role_index, pattern_key in enumerate(sorted(grouped), start=1):
        group = grouped[pattern_key]
        canonicals = {item[3] for item in group}
        roles.append(
            _CompiledReactantRole(
                role_id=f"reactant_role_{role_index}",
                required_count=len(group),
                allow_repeated_component=(
                    len(group) > 1 and len(canonicals) == 1
                ),
                atom_patterns=group[0][0],
                bond_patterns=group[0][1],
                occurrence_map_numbers=tuple(
                    item[2] for item in sorted(group, key=lambda value: value[2])
                ),
            )
        )
    return tuple(roles)


def _atom_matches_role_pattern(atom: Any, pattern: _RoleAtomPattern) -> bool:
    return bool(
        atom.GetSymbol() == pattern.element
        and int(atom.GetFormalCharge()) == pattern.formal_charge
        and bool(atom.GetIsAromatic()) == pattern.aromatic
        and int(atom.GetTotalNumHs(includeNeighbors=True))
        >= pattern.minimum_hydrogen_count
    )


def _role_candidates(
    role: _CompiledReactantRole,
    components: Sequence[Any],
) -> Tuple[_RoleCandidate, ...]:
    """Enumerate bounded atom bindings for one compiled reactant role."""
    candidates = []
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        atom_options = tuple(
            tuple(
                int(atom.GetIdx())
                for atom in molecule.GetAtoms()
                if _atom_matches_role_pattern(atom, pattern)
            )
            for pattern in role.atom_patterns
        )
        if not atom_options or any(not options for options in atom_options):
            continue
        for assignment in itertools.product(*atom_options):
            if len(set(assignment)) != len(assignment):
                continue
            if any(
                (
                    molecule.GetBondBetweenAtoms(
                        assignment[bond.atom_slot_1],
                        assignment[bond.atom_slot_2],
                    )
                    is None
                    or str(
                        molecule.GetBondBetweenAtoms(
                            assignment[bond.atom_slot_1],
                            assignment[bond.atom_slot_2],
                        ).GetBondType()
                    ).upper()
                    != bond.order
                )
                for bond in role.bond_patterns
            ):
                continue
            candidates.append(
                _RoleCandidate(
                    component_index=int(component.component_index),
                    atom_indices=tuple(int(value) for value in assignment),
                )
            )
            if len(candidates) >= 64:
                break
        if len(candidates) >= 64:
            break
    return tuple(
        sorted(
            set(candidates),
            key=lambda candidate: (
                candidate.component_index,
                candidate.atom_indices,
            ),
        )
    )


def _role_assignments(
    role: _CompiledReactantRole,
    candidates: Sequence[_RoleCandidate],
) -> Tuple[Tuple[_RoleCandidate, ...], ...]:
    if not candidates:
        return ()
    if role.required_count == 1:
        return tuple((candidate,) for candidate in candidates)
    if role.allow_repeated_component:
        assignments = itertools.combinations_with_replacement(
            candidates,
            role.required_count,
        )
    else:
        assignments = itertools.combinations(candidates, role.required_count)
    return tuple(itertools.islice(assignments, 128))


def _bond_order(molecule: Any, atom_1: int, atom_2: int) -> Optional[str]:
    bond = molecule.GetBondBetweenAtoms(int(atom_1), int(atom_2))
    return str(bond.GetBondType()).upper() if bond is not None else None


def _set_template_bond(
    molecule: Any,
    atom_1: int,
    atom_2: int,
    *,
    before: Optional[str],
    after: Optional[str],
) -> bool:
    current = _bond_order(molecule, atom_1, atom_2)
    if current != before:
        return False
    if current is not None:
        molecule.RemoveBond(int(atom_1), int(atom_2))
    if after is not None:
        molecule.AddBond(int(atom_1), int(atom_2), bond_type(after))
    return True


def _canonical_fragment_smiles(molecule: Any) -> Tuple[str, ...]:
    from rdkit import Chem

    values = []
    for fragment in Chem.GetMolFrags(
        molecule,
        asMols=True,
        sanitizeFrags=True,
    ):
        for atom in fragment.GetAtoms():
            atom.SetAtomMapNum(0)
        values.append(
            str(
                Chem.MolToSmiles(
                    fragment,
                    canonical=True,
                    isomericSmiles=True,
                )
            )
        )
    return tuple(sorted(values))


def _execute_template_assignment(
    template: ReactionTemplate,
    roles: Sequence[_CompiledReactantRole],
    assignment: Sequence[Sequence[_RoleCandidate]],
    parsed_query: Any,
    observed_products: set[str],
) -> Optional[_TemplateReconstruction]:
    """Apply stored normalized edits to bound query reactant instances."""
    from rdkit import Chem

    components = {
        int(component.component_index): component
        for component in parsed_query.reactants
    }
    source = None
    map_bindings: dict[int, int] = {}
    offset = 0
    used_components: list[int] = []
    for role, role_assignment in zip(roles, assignment):
        for occurrence_index, candidate in enumerate(role_assignment):
            component = components.get(candidate.component_index)
            molecule = (
                parse_smiles(component.input_smiles)
                if component is not None
                else None
            )
            if molecule is None:
                return None
            source = molecule if source is None else Chem.CombineMols(source, molecule)
            occurrence_maps = role.occurrence_map_numbers[occurrence_index]
            if len(occurrence_maps) != len(candidate.atom_indices):
                return None
            for slot, map_number in enumerate(occurrence_maps):
                map_bindings[int(map_number)] = (
                    offset + int(candidate.atom_indices[slot])
                )
            offset += int(molecule.GetNumAtoms())
            used_components.append(candidate.component_index)
    if source is None:
        return None
    rw = Chem.RWMol(source)
    negative = tuple(
        edit
        for edit in template.edits
        if edit.edit_type in {"broken", "order_changed"}
    )
    positive = tuple(
        edit for edit in template.edits if edit.edit_type == "formed"
    )
    for edit in negative:
        if edit.atom_2 is None:
            return None
        atom_1 = map_bindings.get(int(edit.atom_1.atom_map_number))
        atom_2 = map_bindings.get(int(edit.atom_2.atom_map_number))
        if atom_1 is None or atom_2 is None or not _set_template_bond(
            rw,
            atom_1,
            atom_2,
            before=edit.old_order,
            after=edit.new_order,
        ):
            return None
    for edit in positive:
        if edit.atom_2 is None:
            return None
        atom_1 = map_bindings.get(int(edit.atom_1.atom_map_number))
        atom_2 = map_bindings.get(int(edit.atom_2.atom_map_number))
        if atom_1 is None or atom_2 is None or not _set_template_bond(
            rw,
            atom_1,
            atom_2,
            before=None,
            after=edit.new_order,
        ):
            return None
    hydrogen_deltas: dict[int, int] = {}
    for edit in template.edits:
        if edit.edit_type != "hydrogen_change":
            continue
        atom_index = map_bindings.get(int(edit.atom_1.atom_map_number))
        if atom_index is None:
            return None
        delta = (
            -1
            if edit.old_order is not None and edit.new_order is None
            else 1
            if edit.old_order is None and edit.new_order is not None
            else 0
        )
        hydrogen_deltas[atom_index] = (
            hydrogen_deltas.get(atom_index, 0) + delta
        )
    for atom_index, delta in hydrogen_deltas.items():
        if not set_total_hydrogens(source, rw, atom_index, delta):
            return None
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        fragments = _canonical_fragment_smiles(product)
    except Exception:
        return None
    matched = sorted(observed_products.intersection(fragments))
    if not matched:
        return None
    inferred_multiplicity = len(used_components) != len(set(used_components))
    return _TemplateReconstruction(
        predicted_product_smiles=matched[0],
        inferred_multiplicity=inferred_multiplicity,
    )


def _exact_template_reconstructions(
    template: ReactionTemplate,
    parsed_query: Any,
) -> Tuple[_TemplateReconstruction, ...]:
    """Return distinct exact main-product reconstructions for one template."""
    roles = _compile_reactant_roles(template)
    if not roles:
        return ()
    role_options = tuple(
        _role_assignments(
            role,
            _role_candidates(role, parsed_query.reactants),
        )
        for role in roles
    )
    if any(not options for options in role_options):
        return ()
    observed_products = {
        _unmapped_canonical_smiles(component.input_smiles)
        for component in parsed_query.products
    }
    outcomes = []
    for assignment in itertools.islice(
        itertools.product(*role_options),
        256,
    ):
        # Different reference components must bind different supplied
        # components. Reuse within one explicitly repeatable role is the only
        # permitted multiplicity hypothesis.
        component_owners: dict[int, str] = {}
        valid = True
        for role, role_assignment in zip(roles, assignment):
            for candidate in role_assignment:
                owner = component_owners.get(candidate.component_index)
                if owner is not None and owner != role.role_id:
                    valid = False
                    break
                component_owners[candidate.component_index] = role.role_id
            if not valid:
                break
        if not valid:
            continue
        outcome = _execute_template_assignment(
            template,
            roles,
            assignment,
            parsed_query,
            observed_products,
        )
        if outcome is not None:
            outcomes.append(outcome)
    return tuple(
        sorted(
            set(outcomes),
            key=lambda outcome: (
                outcome.predicted_product_smiles,
                outcome.inferred_multiplicity,
            ),
        )
    )


def _center_transition_requirements(
    template: ReactionTemplate,
) -> Tuple[_CenterTransitionRequirement, ...]:
    """Extract minimal edited-centre states without unchanged substituents."""
    parsed = parse_reaction_smiles(template.mapped_reference_reaction)
    before_states = _mapped_atom_states(parsed.reactants)
    after_states = _mapped_atom_states(parsed.products)
    incidence: dict[int, int] = {}
    atoms: dict[int, ReactionTemplateAtom] = {}
    for edit in template.edits:
        for atom in (edit.atom_1, edit.atom_2):
            if atom is None:
                continue
            map_number = int(atom.atom_map_number)
            atoms[map_number] = atom
            incidence[map_number] = incidence.get(map_number, 0) + 1
    if not incidence:
        return ()
    maximum = max(incidence.values())
    center_maps = {
        map_number
        for map_number, count in incidence.items()
        if count == maximum
    }
    before: dict[int, list[str]] = {value: [] for value in center_maps}
    after: dict[int, list[str]] = {value: [] for value in center_maps}
    for edit in template.edits:
        left_map = int(edit.atom_1.atom_map_number)
        if edit.atom_2 is None:
            if left_map in center_maps:
                if edit.old_order is not None:
                    before[left_map].append(
                        _neighbor_token(
                            element="H",
                            aromatic=False,
                            order=edit.old_order,
                        )
                    )
                if edit.new_order is not None:
                    after[left_map].append(
                        _neighbor_token(
                            element="H",
                            aromatic=False,
                            order=edit.new_order,
                        )
                    )
            continue
        right_map = int(edit.atom_2.atom_map_number)
        if left_map in center_maps:
            if edit.old_order is not None:
                before[left_map].append(
                    _template_neighbor_token(
                        atom=edit.atom_2,
                        order=edit.old_order,
                        side_states=before_states,
                    )
                )
            if edit.new_order is not None:
                after[left_map].append(
                    _template_neighbor_token(
                        atom=edit.atom_2,
                        order=edit.new_order,
                        side_states=after_states,
                    )
                )
        if right_map in center_maps:
            if edit.old_order is not None:
                before[right_map].append(
                    _template_neighbor_token(
                        atom=edit.atom_1,
                        order=edit.old_order,
                        side_states=before_states,
                    )
                )
            if edit.new_order is not None:
                after[right_map].append(
                    _template_neighbor_token(
                        atom=edit.atom_1,
                        order=edit.new_order,
                        side_states=after_states,
                    )
                )
    return tuple(
        _CenterTransitionRequirement(
            element=atoms[map_number].element,
            aromatic=atoms[map_number].aromatic,
            before_tokens=tuple(sorted(before[map_number])),
            after_tokens=tuple(sorted(after[map_number])),
        )
        for map_number in sorted(center_maps)
        if before[map_number] and after[map_number]
    )


def _atom_neighbor_tokens(molecule: Any, atom_index: int) -> Tuple[str, ...]:
    atom = molecule.GetAtomWithIdx(int(atom_index))
    tokens = [
        _neighbor_token(
            element=neighbor.GetSymbol(),
            aromatic=bool(neighbor.GetIsAromatic()),
            order=str(
                molecule.GetBondBetweenAtoms(
                    int(atom_index), int(neighbor.GetIdx())
                ).GetBondType()
            ).upper(),
            hydrogen_count=int(
                neighbor.GetTotalNumHs(includeNeighbors=True)
            ),
            external_heavy_count=sum(
                other.GetAtomicNum() > 1
                and int(other.GetIdx()) != int(atom_index)
                for other in neighbor.GetNeighbors()
            ),
        )
        for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() > 1
    ]
    tokens.extend(
        _neighbor_token(element="H", aromatic=False, order="SINGLE")
        for _ in range(int(atom.GetTotalNumHs(includeNeighbors=True)))
    )
    return tuple(sorted(tokens))


def _contains_token_multiset(
    available: Sequence[str],
    required: Sequence[str],
) -> bool:
    counts: dict[str, int] = {}
    for token in available:
        counts[token] = counts.get(token, 0) + 1
    for token in required:
        remaining = counts.get(token, 0)
        if remaining <= 0:
            return False
        counts[token] = remaining - 1
    return True


def _center_candidates(
    components: Sequence[Any],
    requirement: _CenterTransitionRequirement,
    *,
    product_side: bool,
) -> Tuple[tuple[int, int], ...]:
    required = (
        requirement.after_tokens
        if product_side
        else requirement.before_tokens
    )
    candidates = []
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if (
                atom.GetSymbol() != requirement.element
                or bool(atom.GetIsAromatic()) != requirement.aromatic
                or not _contains_token_multiset(
                    _atom_neighbor_tokens(molecule, int(atom.GetIdx())),
                    required,
                )
            ):
                continue
            candidates.append(
                (int(component.component_index), int(atom.GetIdx()))
            )
    return tuple(candidates)


def _has_distinct_center_assignment(
    candidate_sets: Sequence[Sequence[tuple[int, int]]],
) -> bool:
    """Return whether each required center can bind to a distinct query atom."""
    ordered = sorted(candidate_sets, key=len)

    def assign(index: int, used: set[tuple[int, int]]) -> bool:
        if index == len(ordered):
            return True
        for candidate in ordered[index]:
            if candidate in used:
                continue
            used.add(candidate)
            if assign(index + 1, used):
                return True
            used.remove(candidate)
        return False

    return bool(ordered) and all(ordered) and assign(0, set())


def _matches_template_center_transition(
    template: ReactionTemplate,
    parsed_query: Any,
) -> bool:
    """Conservatively match edited before/after center-state multisets."""
    requirements = _center_transition_requirements(template)
    if not requirements:
        return False
    before_candidates = tuple(
        _center_candidates(
            parsed_query.reactants,
            requirement,
            product_side=False,
        )
        for requirement in requirements
    )
    after_candidates = tuple(
        _center_candidates(
            parsed_query.products,
            requirement,
            product_side=True,
        )
        for requirement in requirements
    )
    return _has_distinct_center_assignment(
        before_candidates
    ) and _has_distinct_center_assignment(after_candidates)


def _map_statistics(components: Iterable[Any]) -> tuple[int, int, set[int]]:
    heavy_count = 0
    mapped_heavy_count = 0
    maps: set[int] = set()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            heavy_count += 1
            map_number = int(atom.GetAtomMapNum())
            if map_number:
                mapped_heavy_count += 1
                maps.add(map_number)
    return heavy_count, mapped_heavy_count, maps


def _template_hash_payload(template: ReactionTemplate) -> dict[str, Any]:
    payload = template.to_dict()
    payload.pop("definition_hash", None)
    return payload


def _definition_hash(template: ReactionTemplate) -> str:
    return _digest("RTD1", _template_hash_payload(template), length=40)


def derive_reaction_template(
    mapped_reaction_smiles: str,
    *,
    template_id: str,
    display_name: str,
    family_id: Optional[str] = None,
    aliases: Sequence[str] = (),
    transformation_class: Optional[str] = None,
    status: Literal["draft", "active", "retired"] = "draft",
    provenance: str = "manual_mapped_reference",
    notes: str = "",
) -> ReactionTemplate:
    """Compile one fully mapped reference into a single-event template draft."""
    template_id = str(template_id or "").strip()
    if not _TEMPLATE_ID_RE.fullmatch(template_id):
        raise ReactionTemplateError(
            "template_id must use lowercase snake_case and contain 3-80 characters"
        )
    display_name = str(display_name or "").strip()
    if not display_name:
        raise ReactionTemplateError("display_name is required")
    if status not in _VALID_STATUS:
        raise ReactionTemplateError(f"Unsupported template status: {status}")
    parsed = parse_reaction_smiles(mapped_reaction_smiles)
    if not parsed.valid:
        raise ReactionTemplateError(
            f"Invalid mapped reference reaction: {parsed.error or parsed.warnings}"
        )
    reactant_heavy, reactant_mapped, reactant_maps = _map_statistics(
        parsed.reactants
    )
    product_heavy, product_mapped, product_maps = _map_statistics(parsed.products)
    if reactant_mapped != reactant_heavy:
        raise ReactionTemplateError(
            "Every reactant heavy atom in a template reference must be mapped"
        )
    if product_mapped != product_heavy:
        raise ReactionTemplateError(
            "Every product heavy atom in a template reference must be mapped"
        )
    missing_sources = product_maps - reactant_maps
    if missing_sources:
        raise ReactionTemplateError(
            "Product atom maps lack reactant provenance: "
            + ", ".join(str(value) for value in sorted(missing_sources))
        )
    normalized = normalize_mapped_edits(parsed.reactants, parsed.products)
    if not normalized.valid or not normalized.edits:
        raise ReactionTemplateError(
            "Mapped reference did not yield valid edits: "
            + ", ".join(normalized.warnings or (normalized.evidence,))
        )
    edit_groups = partition_reaction_edits(normalized.edits)
    if len(edit_groups) != 1:
        raise ReactionTemplateError(
            "Reaction templates must contain exactly one connected edit event; "
            f"found {len(edit_groups)}"
        )
    template = ReactionTemplate(
        template_id=template_id,
        display_name=display_name,
        family_id=str(family_id).strip() if family_id else None,
        aliases=tuple(
            sorted(
                {
                    str(alias).strip()
                    for alias in aliases
                    if str(alias).strip()
                }
            )
        ),
        mapped_reference_reaction=str(mapped_reaction_smiles).strip(),
        participants=_participants(parsed),
        edits=tuple(
            ReactionTemplateEdit.from_reaction_edit(edit)
            for edit in normalized.edits
        ),
        edit_fingerprint=reaction_edit_fingerprint(normalized.edits),
        edit_component_count=1,
        edit_archetype=infer_edit_archetype(normalized.edits),
        transformation_class=(
            str(transformation_class).strip()
            if transformation_class
            else None
        ),
        status=status,
        provenance=str(provenance or "manual_mapped_reference").strip(),
        notes=str(notes or "").strip(),
        definition_hash="",
    )
    return replace(template, definition_hash=_definition_hash(template))


def _atom_from_dict(payload: Mapping[str, Any]) -> ReactionTemplateAtom:
    return ReactionTemplateAtom(
        atom_map_number=int(payload["atom_map_number"]),
        element=str(payload["element"]),
        formal_charge=int(payload["formal_charge"]),
        aromatic=bool(payload["aromatic"]),
        hybridization=str(payload["hybridization"]),
    )


def reaction_template_from_dict(payload: Mapping[str, Any]) -> ReactionTemplate:
    """Load and validate one serialized reaction template."""
    try:
        template = ReactionTemplate(
            template_id=str(payload["template_id"]),
            display_name=str(payload["display_name"]),
            family_id=(
                str(payload["family_id"])
                if payload.get("family_id") is not None
                else None
            ),
            aliases=tuple(str(value) for value in payload.get("aliases") or ()),
            mapped_reference_reaction=str(payload["mapped_reference_reaction"]),
            participants=tuple(
                ReactionTemplateParticipant(
                    side=str(item["side"]),  # type: ignore[arg-type]
                    canonical_smiles=str(item["canonical_smiles"]),
                    explicit_count=int(item["explicit_count"]),
                )
                for item in payload.get("participants") or ()
            ),
            edits=tuple(
                ReactionTemplateEdit(
                    edit_type=str(item["edit_type"]),  # type: ignore[arg-type]
                    atom_1=_atom_from_dict(item["atom_1"]),
                    atom_2=(
                        _atom_from_dict(item["atom_2"])
                        if item.get("atom_2") is not None
                        else None
                    ),
                    old_order=(
                        str(item["old_order"])
                        if item.get("old_order") is not None
                        else None
                    ),
                    new_order=(
                        str(item["new_order"])
                        if item.get("new_order") is not None
                        else None
                    ),
                )
                for item in payload.get("edits") or ()
            ),
            edit_fingerprint=str(payload["edit_fingerprint"]),
            edit_component_count=int(payload["edit_component_count"]),
            edit_archetype=str(payload["edit_archetype"]),
            transformation_class=(
                str(payload["transformation_class"])
                if payload.get("transformation_class") is not None
                else None
            ),
            status=str(payload["status"]),  # type: ignore[arg-type]
            provenance=str(payload.get("provenance") or ""),
            notes=str(payload.get("notes") or ""),
            definition_hash=str(payload["definition_hash"]),
            schema_version=str(
                payload.get("schema_version")
                or REACTION_TEMPLATE_SCHEMA_VERSION
            ),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ReactionTemplateError(f"Invalid reaction template record: {exc}") from exc
    errors = validate_reaction_template(template)
    if errors:
        raise ReactionTemplateError("; ".join(errors))
    return template


def validate_reaction_template(template: ReactionTemplate) -> Tuple[str, ...]:
    """Return deterministic validation errors for one template."""
    errors = []
    if template.schema_version != REACTION_TEMPLATE_SCHEMA_VERSION:
        errors.append(
            f"{template.template_id}:unsupported_schema_version:"
            f"{template.schema_version}"
        )
    if not _TEMPLATE_ID_RE.fullmatch(template.template_id):
        errors.append(f"{template.template_id}:invalid_template_id")
    if not template.display_name.strip():
        errors.append(f"{template.template_id}:missing_display_name")
    if template.status not in _VALID_STATUS:
        errors.append(f"{template.template_id}:invalid_status:{template.status}")
    if template.edit_component_count != 1:
        errors.append(
            f"{template.template_id}:not_single_event:"
            f"{template.edit_component_count}"
        )
    if not template.edits:
        errors.append(f"{template.template_id}:missing_edits")
    if not template.edit_fingerprint.startswith("RTE1:"):
        errors.append(f"{template.template_id}:invalid_edit_fingerprint")
    if not template.definition_hash.startswith("RTD1:"):
        errors.append(f"{template.template_id}:invalid_definition_hash")
    elif _definition_hash(replace(template, definition_hash="")) != (
        template.definition_hash
    ):
        errors.append(f"{template.template_id}:definition_hash_mismatch")
    if len(set(template.aliases)) != len(template.aliases):
        errors.append(f"{template.template_id}:duplicate_aliases")
    try:
        compiled = derive_reaction_template(
            template.mapped_reference_reaction,
            template_id=template.template_id,
            display_name=template.display_name,
            family_id=template.family_id,
            aliases=template.aliases,
            transformation_class=template.transformation_class,
            status=template.status,
            provenance=template.provenance,
            notes=template.notes,
        )
    except ReactionTemplateError as exc:
        errors.append(f"{template.template_id}:invalid_mapped_reference:{exc}")
    else:
        for field_name in (
            "participants",
            "edits",
            "edit_fingerprint",
            "edit_component_count",
            "edit_archetype",
        ):
            if getattr(template, field_name) != getattr(compiled, field_name):
                errors.append(
                    f"{template.template_id}:reference_contract_mismatch:"
                    f"{field_name}"
                )
    return tuple(sorted(errors))


def empty_reaction_template_registry() -> dict[str, Any]:
    """Return the canonical empty registry document."""
    return {
        "schema_version": REACTION_TEMPLATE_SCHEMA_VERSION,
        "definition_version": REACTION_TEMPLATE_DEFINITION_VERSION,
        "templates": [],
    }


def load_reaction_template_registry(
    path: str | Path = DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
) -> Tuple[ReactionTemplate, ...]:
    """Load a versioned reaction-template registry."""
    registry_path = Path(path)
    if not registry_path.exists():
        return ()
    try:
        payload = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReactionTemplateError(
            f"Cannot read reaction-template registry {registry_path}: {exc}"
        ) from exc
    if not isinstance(payload, Mapping):
        raise ReactionTemplateError("Reaction-template registry must be an object")
    if payload.get("schema_version") != REACTION_TEMPLATE_SCHEMA_VERSION:
        raise ReactionTemplateError(
            "Unsupported reaction-template registry schema version"
        )
    if payload.get("definition_version") != REACTION_TEMPLATE_DEFINITION_VERSION:
        raise ReactionTemplateError(
            "Unsupported reaction-template registry definition version"
        )
    records = payload.get("templates")
    if not isinstance(records, list):
        raise ReactionTemplateError(
            "Reaction-template registry templates must be a list"
        )
    templates = tuple(reaction_template_from_dict(record) for record in records)
    ids = [template.template_id for template in templates]
    if len(ids) != len(set(ids)):
        raise ReactionTemplateError("Duplicate reaction-template IDs")
    return tuple(sorted(templates, key=lambda item: item.template_id))


def save_reaction_template_registry(
    templates: Sequence[ReactionTemplate],
    path: str | Path = DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
) -> Path:
    """Atomically save validated templates to a versioned registry file."""
    registry_path = Path(path)
    normalized = tuple(sorted(templates, key=lambda item: item.template_id))
    ids = [template.template_id for template in normalized]
    if len(ids) != len(set(ids)):
        raise ReactionTemplateError("Duplicate reaction-template IDs")
    errors = tuple(
        error
        for template in normalized
        for error in validate_reaction_template(template)
    )
    if errors:
        raise ReactionTemplateError("; ".join(errors))
    payload = empty_reaction_template_registry()
    payload["templates"] = [template.to_dict() for template in normalized]
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, indent=2
    ) + "\n"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{registry_path.name}.",
        suffix=".tmp",
        dir=str(registry_path.parent),
        text=True,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, registry_path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    return registry_path


def upsert_reaction_template(
    template: ReactionTemplate,
    path: str | Path = DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
    *,
    replace_existing: bool = False,
) -> Path:
    """Add one template, requiring explicit permission to replace an ID."""
    templates = list(load_reaction_template_registry(path))
    existing_index = next(
        (
            index
            for index, current in enumerate(templates)
            if current.template_id == template.template_id
        ),
        None,
    )
    if existing_index is not None and not replace_existing:
        raise ReactionTemplateError(
            f"Template already exists: {template.template_id}"
        )
    if existing_index is None:
        templates.append(template)
    else:
        templates[existing_index] = template
    return save_reaction_template_registry(templates, path)


def validate_reaction_template_registry(
    path: str | Path = DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
) -> Tuple[str, ...]:
    """Return registry validation errors without raising."""
    try:
        templates = load_reaction_template_registry(path)
    except ReactionTemplateError as exc:
        return (str(exc),)
    errors = [
        error
        for template in templates
        for error in validate_reaction_template(template)
    ]
    aliases: dict[str, str] = {}
    for template in templates:
        for alias in (template.display_name, *template.aliases):
            normalized = alias.casefold().strip()
            owner = aliases.get(normalized)
            if owner is not None and owner != template.template_id:
                errors.append(
                    f"ambiguous_template_alias:{alias}:{owner}:"
                    f"{template.template_id}"
                )
            aliases[normalized] = template.template_id
    return tuple(sorted(set(errors)))


def match_reaction_templates(
    reaction_smiles: str,
    *,
    path: str | Path = DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH,
    include_drafts: bool = False,
) -> ReactionTemplateQueryResult:
    """Match a query-derived edit graph against registered templates.

    This function never copies a stored signature.  The query signature is
    generated by :func:`featurize_reaction`; the template registry contributes
    only structural interpretation candidates.
    """
    from .reaction_api import featurize_reaction

    analysis = featurize_reaction(reaction_smiles)
    signature = analysis.reaction_signature
    edits: Tuple[ReactionEdit, ...] = (
        tuple(signature.edits) if signature is not None else ()
    )
    evidence = analysis.evidence_quality
    warnings = list(analysis.warnings)
    if not edits:
        parsed = parse_reaction_smiles(reaction_smiles)
        if parsed.valid:
            mapped = normalize_mapped_edits(parsed.reactants, parsed.products)
            if mapped.valid:
                edits = mapped.edits
                evidence = mapped.evidence
                warnings.extend(mapped.warnings)
    fingerprint = reaction_edit_fingerprint(edits) if edits else None
    templates = load_reaction_template_registry(path)
    exact_matches = tuple(
        ReactionTemplateMatch(
            template_id=template.template_id,
            display_name=template.display_name,
            family_id=template.family_id,
            status=template.status,
            edit_fingerprint=template.edit_fingerprint,
            definition_hash=template.definition_hash,
            evidence="query_derived_edit_fingerprint",
            confidence=1.0,
            provisional=False,
            predicted_product_smiles=None,
            inferred_multiplicity=False,
        )
        for template in templates
        if fingerprint is not None
        and template.edit_fingerprint == fingerprint
        and (include_drafts or template.status == "active")
    )
    matches = exact_matches
    if fingerprint is None and analysis.valid:
        parsed_query = parse_reaction_smiles(reaction_smiles)
        reconstruction_matches = []
        for template in templates:
            if not (include_drafts or template.status == "active"):
                continue
            outcomes = _exact_template_reconstructions(template, parsed_query)
            if not outcomes:
                continue
            inferred_multiplicity = all(
                outcome.inferred_multiplicity for outcome in outcomes
            )
            preferred_outcome = next(
                (
                    outcome
                    for outcome in outcomes
                    if not outcome.inferred_multiplicity
                ),
                outcomes[0],
            )
            reconstruction_matches.append(
                ReactionTemplateMatch(
                    template_id=template.template_id,
                    display_name=template.display_name,
                    family_id=template.family_id,
                    status=template.status,
                    edit_fingerprint=template.edit_fingerprint,
                    definition_hash=template.definition_hash,
                    evidence=(
                        "exact_template_reconstruction_with_"
                        "inferred_multiplicity"
                        if inferred_multiplicity
                        else "exact_template_reconstruction"
                    ),
                    confidence=0.85 if inferred_multiplicity else 0.95,
                    provisional=inferred_multiplicity,
                    predicted_product_smiles=(
                        preferred_outcome.predicted_product_smiles
                    ),
                    inferred_multiplicity=inferred_multiplicity,
                )
            )
        if reconstruction_matches:
            matches = tuple(reconstruction_matches)
            candidate_fingerprints = {
                match.edit_fingerprint for match in matches
            }
            candidate_evidence = {match.evidence for match in matches}
            if len(candidate_fingerprints) == 1:
                fingerprint = next(iter(candidate_fingerprints))
                evidence = (
                    next(iter(candidate_evidence))
                    if len(candidate_evidence) == 1
                    else "exact_template_reconstruction"
                )
                warnings.append(
                    "EXACT_MAIN_PRODUCT_RECONSTRUCTION_FROM_TEMPLATE"
                )
                if any(match.provisional for match in matches):
                    warnings.append("INFERRED_REACTANT_MULTIPLICITY")
            else:
                warnings.append("AMBIGUOUS_TEMPLATE_RECONSTRUCTIONS")
        else:
            provisional_matches = tuple(
                ReactionTemplateMatch(
                    template_id=template.template_id,
                    display_name=template.display_name,
                    family_id=template.family_id,
                    status=template.status,
                    edit_fingerprint=template.edit_fingerprint,
                    definition_hash=template.definition_hash,
                    evidence="template_center_transition_hypothesis",
                    confidence=0.7,
                    provisional=True,
                    predicted_product_smiles=None,
                    inferred_multiplicity=False,
                )
                for template in templates
                if (include_drafts or template.status == "active")
                and _matches_template_center_transition(template, parsed_query)
            )
            matches = provisional_matches
            candidate_fingerprints = {
                match.edit_fingerprint for match in provisional_matches
            }
            if len(candidate_fingerprints) == 1:
                fingerprint = next(iter(candidate_fingerprints))
                evidence = "template_center_transition_hypothesis"
                warnings.append(
                    "PROVISIONAL_TEMPLATE_MATCH_WITHOUT_ATOM_PROVENANCE"
                )
            elif len(candidate_fingerprints) > 1:
                warnings.append("AMBIGUOUS_TEMPLATE_CENTER_TRANSITIONS")
    return ReactionTemplateQueryResult(
        reaction_smiles=reaction_smiles,
        valid=bool(analysis.valid),
        evidence=evidence,
        edit_fingerprint=fingerprint,
        signature_id=signature.signature_id if signature else None,
        matches=matches,
        warnings=tuple(sorted(set(warnings))),
        error=analysis.error,
    )


__all__ = [
    "DEFAULT_REACTION_TEMPLATE_REGISTRY_PATH",
    "REACTION_TEMPLATE_DEFINITION_VERSION",
    "REACTION_TEMPLATE_SCHEMA_VERSION",
    "ReactionTemplate",
    "ReactionTemplateAtom",
    "ReactionTemplateEdit",
    "ReactionTemplateError",
    "ReactionTemplateMatch",
    "ReactionTemplateParticipant",
    "ReactionTemplateQueryResult",
    "derive_reaction_template",
    "empty_reaction_template_registry",
    "load_reaction_template_registry",
    "match_reaction_templates",
    "reaction_edit_fingerprint",
    "reaction_template_from_dict",
    "save_reaction_template_registry",
    "upsert_reaction_template",
    "validate_reaction_template",
    "validate_reaction_template_registry",
]
