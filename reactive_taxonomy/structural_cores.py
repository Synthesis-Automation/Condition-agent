"""Deterministic structural-core observations for normalized target graphs.

The module observes bounded target subgraphs.  It does not infer a preferred
synthetic strategy, retrieve precedents, or modify retrosynthesis ranking.
Target atom references are scoped to one normalized molecule and are not
presented as observed reaction atom mapping.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from .chemistry.rdkit_utils import (
    mol_to_canonical_smiles,
    parse_smiles,
    prepare_fragment_serialization_copy,
)


STRUCTURAL_CORE_SCHEMA_VERSION = "structural_core_observation.v1"
STRUCTURAL_CORE_OBSERVATION_DEFINITION_ID = "structural_core_observations.v1"
STRUCTURAL_CORE_MATCHING_DEFINITION_ID = "structural_core_matching.v1"

_DEFINITION_DIRECTORY = Path(__file__).with_name("definitions")
_OBSERVATION_DEFINITION_PATH = (
    _DEFINITION_DIRECTORY / "structural_core_observations.v1.json"
)
_MATCHING_DEFINITION_PATH = (
    _DEFINITION_DIRECTORY / "structural_core_matching.v1.json"
)

StructuralCoreKind = Literal[
    "scaffold_backbone",
    "bridge_interface",
    "linker_region",
    "carbon_framework",
    "stereo_backbone_region",
    "ring_system",
]

_CORE_KINDS = frozenset(
    {
        "scaffold_backbone",
        "bridge_interface",
        "linker_region",
        "carbon_framework",
        "stereo_backbone_region",
        "ring_system",
    }
)
_OBSERVATION_DEFINITION_FIELDS = {
    "definition_id",
    "definition_version",
    "schema_version",
    "maximum_observations",
    "component_policy",
    "minimum_candidate_heavy_atoms",
    "minimum_carbon_framework_atoms",
    "minimum_bridge_side_heavy_atoms",
    "minimum_bridge_balance_fraction",
    "bridge_neighborhood_radius",
    "stereo_neighborhood_radius",
    "maximum_coverage_by_kind",
    "maximum_per_kind",
    "simple_ring_policy",
    "kind_priority",
}
_MATCHING_DEFINITION_FIELDS = {
    "definition_id",
    "definition_version",
    "schema_version",
    "weisfeiler_lehman_iterations",
    "exact_atom_fields",
    "typed_atom_fields",
    "shape_atom_fields",
    "exact_attachment_fields",
    "typed_attachment_fields",
    "shape_attachment_fields",
}
_EXPECTED_MATCHING_FIELDS = {
    "exact_atom_fields": (
        "element",
        "isotope",
        "formal_charge",
        "aromatic",
        "hybridization",
        "chiral_tag",
        "total_hydrogens",
    ),
    "typed_atom_fields": (
        "element",
        "formal_charge",
        "aromatic",
        "hybridization",
    ),
    "shape_atom_fields": ("internal_degree", "ring_membership"),
    "exact_attachment_fields": (
        "anchor_label",
        "bond_order",
        "external_element",
        "external_formal_charge",
        "external_aromatic",
    ),
    "typed_attachment_fields": (
        "anchor_label",
        "bond_order",
        "external_element",
        "external_aromatic",
    ),
    "shape_attachment_fields": (
        "anchor_label",
        "bond_order",
    ),
}


def _stable_id(prefix: str, value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return f"{prefix}:{hashlib.sha256(encoded).hexdigest()[:20]}"


def _positive_integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def validate_structural_core_observation_definition(
    value: Mapping[str, Any],
) -> list[str]:
    """Return stable validation errors for the core-observation policy."""

    errors: list[str] = []
    if set(value) != _OBSERVATION_DEFINITION_FIELDS:
        errors.append("structural_core_observations:invalid_fields")
    if value.get("definition_id") != STRUCTURAL_CORE_OBSERVATION_DEFINITION_ID:
        errors.append("structural_core_observations:invalid_definition_id")
    if value.get("schema_version") != "1.0":
        errors.append("structural_core_observations:invalid_schema_version")
    if not str(value.get("definition_version") or "").strip():
        errors.append("structural_core_observations:missing_definition_version")
    if value.get("component_policy") != "largest_heavy_atom_component":
        errors.append("structural_core_observations:invalid_component_policy")
    for field in (
        "maximum_observations",
        "minimum_candidate_heavy_atoms",
        "minimum_carbon_framework_atoms",
        "minimum_bridge_side_heavy_atoms",
    ):
        if not _positive_integer(value.get(field)):
            errors.append(f"structural_core_observations:invalid_{field}")
    for field in (
        "bridge_neighborhood_radius",
        "stereo_neighborhood_radius",
    ):
        raw = value.get(field)
        if (
            not isinstance(raw, int)
            or isinstance(raw, bool)
            or raw < 0
            or raw > 5
        ):
            errors.append(f"structural_core_observations:invalid_{field}")
    balance = value.get("minimum_bridge_balance_fraction")
    if (
        not isinstance(balance, (int, float))
        or isinstance(balance, bool)
        or not 0.0 < float(balance) <= 0.5
    ):
        errors.append(
            "structural_core_observations:invalid_minimum_bridge_balance_fraction"
        )
    coverage = value.get("maximum_coverage_by_kind")
    if (
        not isinstance(coverage, dict)
        or set(coverage) != _CORE_KINDS
        or any(
            not isinstance(limit, (int, float))
            or isinstance(limit, bool)
            or not 0.0 < float(limit) <= 1.0
            for limit in coverage.values()
        )
    ):
        errors.append("structural_core_observations:invalid_coverage_limits")
    kind_limits = value.get("maximum_per_kind")
    if (
        not isinstance(kind_limits, dict)
        or set(kind_limits) != _CORE_KINDS
        or any(not _positive_integer(limit) for limit in kind_limits.values())
    ):
        errors.append("structural_core_observations:invalid_kind_limits")
    if value.get("simple_ring_policy") != "distinctive_only":
        errors.append("structural_core_observations:invalid_simple_ring_policy")
    priority = value.get("kind_priority")
    if (
        not isinstance(priority, list)
        or len(priority) != len(set(priority))
        or set(priority) != _CORE_KINDS
    ):
        errors.append("structural_core_observations:invalid_kind_priority")
    return errors


def validate_structural_core_matching_definition(
    value: Mapping[str, Any],
) -> list[str]:
    """Return stable validation errors for structural-core key generation."""

    errors: list[str] = []
    if set(value) != _MATCHING_DEFINITION_FIELDS:
        errors.append("structural_core_matching:invalid_fields")
    if value.get("definition_id") != STRUCTURAL_CORE_MATCHING_DEFINITION_ID:
        errors.append("structural_core_matching:invalid_definition_id")
    if value.get("schema_version") != "1.0":
        errors.append("structural_core_matching:invalid_schema_version")
    if not str(value.get("definition_version") or "").strip():
        errors.append("structural_core_matching:missing_definition_version")
    iterations = value.get("weisfeiler_lehman_iterations")
    if (
        not isinstance(iterations, int)
        or isinstance(iterations, bool)
        or iterations < 1
        or iterations > 16
    ):
        errors.append("structural_core_matching:invalid_iterations")
    for field, expected in _EXPECTED_MATCHING_FIELDS.items():
        raw = value.get(field)
        if not isinstance(raw, list) or tuple(raw) != expected:
            errors.append(f"structural_core_matching:invalid_{field}")
    return errors


@lru_cache(maxsize=1)
def load_structural_core_observation_definition() -> Mapping[str, Any]:
    """Load and validate the canonical structural-core observation policy."""

    value = json.loads(_OBSERVATION_DEFINITION_PATH.read_text(encoding="utf-8"))
    errors = validate_structural_core_observation_definition(value)
    if errors:
        raise ValueError("; ".join(errors))
    return value


@lru_cache(maxsize=1)
def load_structural_core_matching_definition() -> Mapping[str, Any]:
    """Load and validate the canonical structural-core matching policy."""

    value = json.loads(_MATCHING_DEFINITION_PATH.read_text(encoding="utf-8"))
    errors = validate_structural_core_matching_definition(value)
    if errors:
        raise ValueError("; ".join(errors))
    return value


@dataclass(frozen=True)
class MoleculeAtomReference:
    """One atom reference scoped to a normalized target molecule."""

    molecule_id: str
    component_index: int
    atom_index: int
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    environment_id: str
    schema_version: str = STRUCTURAL_CORE_SCHEMA_VERSION


@dataclass(frozen=True)
class StructuralCoreObservation:
    """One bounded graph-derived target subgraph observation."""

    core_observation_id: str
    rank: int
    molecule_id: str
    component_index: int
    kind: StructuralCoreKind
    atom_references: tuple[MoleculeAtomReference, ...]
    attachment_bond_keys: tuple[str, ...]
    focus_atom_indices: tuple[int, ...]
    focus_bond_keys: tuple[str, ...]
    structural_exact_key: str
    structural_typed_key: str
    structural_shape_key: str
    topology_tokens: tuple[str, ...]
    descriptor_values: tuple[tuple[str, float], ...]
    evidence_codes: tuple[str, ...]
    warnings: tuple[str, ...]
    definition_id: str
    definition_version: str
    matching_definition_id: str
    matching_definition_version: str
    schema_version: str = STRUCTURAL_CORE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.rank < 1:
            raise ValueError("structural-core rank must be positive")
        if self.kind not in _CORE_KINDS:
            raise ValueError(f"unsupported structural-core kind: {self.kind}")
        if not self.atom_references:
            raise ValueError("structural-core observation requires atoms")
        if any(
            reference.molecule_id != self.molecule_id
            or reference.component_index != self.component_index
            for reference in self.atom_references
        ):
            raise ValueError("structural-core atom reference scope is inconsistent")
        if not set(self.focus_atom_indices) <= set(self.atom_indices):
            raise ValueError("structural-core focus atoms must belong to the core")

    @property
    def atom_indices(self) -> tuple[int, ...]:
        """Return target-scoped atom indices for display and verification."""

        return tuple(reference.atom_index for reference in self.atom_references)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible observation."""

        return {
            "core_observation_id": self.core_observation_id,
            "rank": self.rank,
            "molecule_id": self.molecule_id,
            "component_index": self.component_index,
            "kind": self.kind,
            "atom_references": [asdict(item) for item in self.atom_references],
            "attachment_bond_keys": list(self.attachment_bond_keys),
            "focus_atom_indices": list(self.focus_atom_indices),
            "focus_bond_keys": list(self.focus_bond_keys),
            "structural_exact_key": self.structural_exact_key,
            "structural_typed_key": self.structural_typed_key,
            "structural_shape_key": self.structural_shape_key,
            "topology_tokens": list(self.topology_tokens),
            "descriptor_values": dict(self.descriptor_values),
            "evidence_codes": list(self.evidence_codes),
            "warnings": list(self.warnings),
            "definition_id": self.definition_id,
            "definition_version": self.definition_version,
            "matching_definition_id": self.matching_definition_id,
            "matching_definition_version": self.matching_definition_version,
            "schema_version": self.schema_version,
        }


@dataclass(frozen=True)
class StructuralCoreAnalysis:
    """Deterministic structural-core observations for one target."""

    input_smiles: str
    canonical_smiles: str | None
    molecule_id: str | None
    valid: bool
    observations: tuple[StructuralCoreObservation, ...] = ()
    warnings: tuple[str, ...] = ()
    error: str | None = None
    definition_id: str = STRUCTURAL_CORE_OBSERVATION_DEFINITION_ID
    definition_version: str = ""
    matching_definition_id: str = STRUCTURAL_CORE_MATCHING_DEFINITION_ID
    matching_definition_version: str = ""
    schema_version: str = STRUCTURAL_CORE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible target analysis."""

        return {
            "input_smiles": self.input_smiles,
            "canonical_smiles": self.canonical_smiles,
            "molecule_id": self.molecule_id,
            "valid": self.valid,
            "observations": [item.to_dict() for item in self.observations],
            "warnings": list(self.warnings),
            "error": self.error,
            "definition_id": self.definition_id,
            "definition_version": self.definition_version,
            "matching_definition_id": self.matching_definition_id,
            "matching_definition_version": self.matching_definition_version,
            "schema_version": self.schema_version,
        }


@dataclass(frozen=True)
class _CandidateCore:
    kind: StructuralCoreKind
    atom_indices: frozenset[int]
    evidence_codes: tuple[str, ...]
    selection_score: float = 0.0
    focus_atom_indices: tuple[int, ...] = ()
    focus_bond_pairs: tuple[tuple[int, int], ...] = ()


def _fragment_smiles(molecule: Any, atom_indices: Sequence[int]) -> str:
    from rdkit import Chem  # type: ignore

    copied = prepare_fragment_serialization_copy(molecule, atom_indices)
    if copied is None:
        return ""
    try:
        return str(
            Chem.MolFragmentToSmiles(
                copied,
                atomsToUse=sorted(int(value) for value in atom_indices),
                canonical=True,
                isomericSmiles=True,
            )
        )
    except Exception:
        return ""


def _component_records(molecule: Any) -> tuple[tuple[int, tuple[int, ...], str], ...]:
    from rdkit import Chem  # type: ignore

    records = []
    for indices in Chem.GetMolFrags(molecule, asMols=False, sanitizeFrags=False):
        atom_indices = tuple(sorted(int(value) for value in indices))
        heavy_count = sum(
            molecule.GetAtomWithIdx(index).GetAtomicNum() > 1
            for index in atom_indices
        )
        records.append((heavy_count, atom_indices, _fragment_smiles(molecule, atom_indices)))
    ordered = sorted(records, key=lambda item: (-item[0], item[2], item[1]))
    return tuple(
        (component_index, atom_indices, smiles)
        for component_index, (_, atom_indices, smiles) in enumerate(ordered)
    )


def _ring_systems(molecule: Any, allowed: frozenset[int]) -> tuple[frozenset[int], ...]:
    systems = [
        set(int(index) for index in ring if int(index) in allowed)
        for ring in molecule.GetRingInfo().AtomRings()
    ]
    systems = [value for value in systems if value]
    merged: list[set[int]] = []
    while systems:
        current = systems.pop(0)
        changed = True
        while changed:
            changed = False
            remaining = []
            for candidate in systems:
                if current.intersection(candidate):
                    current.update(candidate)
                    changed = True
                else:
                    remaining.append(candidate)
            systems = remaining
        merged.append(current)
    return tuple(
        frozenset(value)
        for value in sorted(merged, key=lambda item: (-len(item), tuple(sorted(item))))
    )


def _murcko_scaffold_atoms(
    molecule: Any,
    allowed: frozenset[int],
) -> frozenset[int]:
    """Return target atom indices retained by the Bemis-Murcko scaffold."""

    from rdkit import Chem  # type: ignore
    from rdkit.Chem.Scaffolds import MurckoScaffold  # type: ignore

    mapped = Chem.Mol(molecule)
    for atom in mapped.GetAtoms():
        atom.SetAtomMapNum(int(atom.GetIdx()) + 1)
    scaffold = MurckoScaffold.GetScaffoldForMol(mapped)
    return frozenset(
        int(atom.GetAtomMapNum()) - 1
        for atom in scaffold.GetAtoms()
        if atom.GetAtomMapNum() and int(atom.GetAtomMapNum()) - 1 in allowed
    )


def _acyclic_main_skeleton(
    molecule: Any,
    allowed: frozenset[int],
) -> frozenset[int]:
    """Return the deterministic heavy-atom graph diameter for an acyclic target."""

    from rdkit import Chem  # type: ignore

    paths = [
        tuple(int(value) for value in Chem.GetShortestPath(molecule, left, right))
        for position, left in enumerate(sorted(allowed))
        for right in sorted(allowed)[position + 1 :]
    ]
    if not paths:
        return allowed
    return frozenset(min(paths, key=lambda path: (-len(path), path)))


def _direct_ring_linkers(
    molecule: Any,
    ring_systems: Sequence[frozenset[int]],
) -> tuple[
    tuple[frozenset[int], float, tuple[str, ...], tuple[tuple[int, int], ...]],
    ...,
]:
    """Return minimal paths directly connecting pairs of ring systems."""

    from rdkit import Chem  # type: ignore

    all_ring_atoms = frozenset().union(*ring_systems) if ring_systems else frozenset()
    records: dict[
        frozenset[int],
        tuple[float, tuple[str, ...], tuple[tuple[int, int], ...]],
    ] = {}
    for left_position, left in enumerate(ring_systems):
        for right in ring_systems[left_position + 1 :]:
            paths = [
                tuple(
                    int(value)
                    for value in Chem.GetShortestPath(molecule, left_atom, right_atom)
                )
                for left_atom in sorted(left)
                for right_atom in sorted(right)
            ]
            direct = [
                path
                for path in paths
                if path
                and not set(path[1:-1]).intersection(all_ring_atoms)
            ]
            if not direct:
                continue
            path = min(direct, key=lambda value: (len(value), value))
            selected = frozenset(path)
            internal = path[1:-1]
            evidence = {"DIRECT_RING_LINKER"}
            if len(path) == 2:
                evidence.add("DIRECT_RING_TO_RING_BOND")
            if any(molecule.GetAtomWithIdx(index).GetAtomicNum() != 6 for index in internal):
                evidence.add("HETEROATOM_LINKER")
            score = (
                min(len(left), len(right)) / max(1, max(len(left), len(right)))
                + 1.0 / len(path)
            )
            focus_bonds = tuple(
                (min(first, second), max(first, second))
                for first, second in zip(path, path[1:])
            )
            previous = records.get(selected)
            if previous is None or score > previous[0]:
                records[selected] = (
                    score,
                    tuple(sorted(evidence)),
                    focus_bonds,
                )
    return tuple(
        (indices, score, evidence, focus_bonds)
        for indices, (score, evidence, focus_bonds) in sorted(
            records.items(),
            key=lambda item: (-item[1][0], tuple(sorted(item[0]))),
        )
    )


def _induced_components(
    molecule: Any,
    atom_indices: frozenset[int],
) -> tuple[frozenset[int], ...]:
    pending = set(atom_indices)
    components = []
    while pending:
        start = min(pending)
        selected = {start}
        frontier = {start}
        pending.remove(start)
        while frontier:
            next_frontier = {
                int(neighbor.GetIdx())
                for index in frontier
                for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
                if int(neighbor.GetIdx()) in pending
            }
            selected.update(next_frontier)
            pending.difference_update(next_frontier)
            frontier = next_frontier
        components.append(frozenset(selected))
    return tuple(
        sorted(components, key=lambda item: (-len(item), tuple(sorted(item))))
    )


def _bond_partition(
    molecule: Any,
    allowed: frozenset[int],
    left: int,
    right: int,
) -> tuple[frozenset[int], frozenset[int]]:
    selected = {left}
    frontier = {left}
    blocked = frozenset({left, right})
    while frontier:
        next_frontier = {
            int(neighbor.GetIdx())
            for index in frontier
            for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
            if int(neighbor.GetIdx()) in allowed
            and int(neighbor.GetIdx()) not in selected
            and frozenset({index, int(neighbor.GetIdx())}) != blocked
        }
        selected.update(next_frontier)
        frontier = next_frontier
    left_side = frozenset(selected)
    return left_side, frozenset(allowed.difference(left_side))


def _neighborhood(
    molecule: Any,
    seeds: Sequence[int],
    allowed: frozenset[int],
    radius: int,
) -> frozenset[int]:
    selected = {int(value) for value in seeds if int(value) in allowed}
    frontier = set(selected)
    for _ in range(radius):
        next_frontier = set()
        for index in frontier:
            atom = molecule.GetAtomWithIdx(index)
            next_frontier.update(
                int(neighbor.GetIdx())
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() > 1
                and int(neighbor.GetIdx()) in allowed
                and int(neighbor.GetIdx()) not in selected
            )
        if not next_frontier:
            break
        selected.update(next_frontier)
        frontier = next_frontier
    return frozenset(selected)


def _merge_overlapping(values: Sequence[frozenset[int]]) -> tuple[frozenset[int], ...]:
    pending = [set(value) for value in values if value]
    merged = []
    while pending:
        current = pending.pop(0)
        changed = True
        while changed:
            changed = False
            remaining = []
            for candidate in pending:
                if current.intersection(candidate):
                    current.update(candidate)
                    changed = True
                else:
                    remaining.append(candidate)
            pending = remaining
        merged.append(frozenset(current))
    return tuple(sorted(merged, key=lambda item: (-len(item), tuple(sorted(item)))))


def _attachment_bonds(
    molecule: Any,
    selected: frozenset[int],
) -> tuple[Any, ...]:
    return tuple(
        sorted(
            (
                bond
                for bond in molecule.GetBonds()
                if (
                    int(bond.GetBeginAtomIdx()) in selected
                    and int(bond.GetEndAtomIdx()) not in selected
                )
                or (
                    int(bond.GetEndAtomIdx()) in selected
                    and int(bond.GetBeginAtomIdx()) not in selected
                )
            ),
            key=lambda bond: (
                min(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
                max(int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())),
            ),
        )
    )


def _bond_token(bond: Any) -> str:
    return (
        f"{float(bond.GetBondTypeAsDouble()):.1f}:"
        f"{'aromatic' if bond.GetIsAromatic() else 'aliphatic'}"
    )


def _atom_environment_id(atom: Any) -> str:
    neighbors = sorted(
        (
            neighbor.GetSymbol(),
            int(neighbor.GetFormalCharge()),
            bool(neighbor.GetIsAromatic()),
            _bond_token(atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())),
        )
        for neighbor in atom.GetNeighbors()
        if neighbor.GetAtomicNum() > 1
    )
    return _stable_id(
        "ENV1",
        {
            "element": atom.GetSymbol(),
            "formal_charge": int(atom.GetFormalCharge()),
            "aromatic": bool(atom.GetIsAromatic()),
            "hybridization": str(atom.GetHybridization()),
            "neighbors": neighbors,
        },
    )


def _initial_atom_label(
    atom: Any,
    selected: frozenset[int],
    ring_memberships: Mapping[int, int],
    level: str,
) -> str:
    internal_degree = sum(
        int(neighbor.GetIdx()) in selected and neighbor.GetAtomicNum() > 1
        for neighbor in atom.GetNeighbors()
    )
    if level == "exact":
        value = (
            atom.GetSymbol(),
            int(atom.GetIsotope()),
            int(atom.GetFormalCharge()),
            bool(atom.GetIsAromatic()),
            str(atom.GetHybridization()),
            str(atom.GetChiralTag()),
            int(atom.GetTotalNumHs()),
        )
    elif level == "typed":
        value = (
            atom.GetSymbol(),
            int(atom.GetFormalCharge()),
            bool(atom.GetIsAromatic()),
            str(atom.GetHybridization()),
        )
    else:
        value = (internal_degree, int(ring_memberships.get(int(atom.GetIdx()), 0) > 0))
    return _stable_id("AL1", value)


def _refined_labels(
    molecule: Any,
    selected: frozenset[int],
    ring_memberships: Mapping[int, int],
    level: str,
    iterations: int,
) -> dict[int, str]:
    labels = {
        index: _initial_atom_label(
            molecule.GetAtomWithIdx(index),
            selected,
            ring_memberships,
            level,
        )
        for index in selected
    }
    for _ in range(iterations):
        labels = {
            index: _stable_id(
                "WL1",
                {
                    "self": labels[index],
                    "neighbors": sorted(
                        (
                            _bond_token(bond),
                            labels[
                                int(
                                    bond.GetOtherAtomIdx(index)
                                )
                            ],
                        )
                        for bond in molecule.GetAtomWithIdx(index).GetBonds()
                        if int(bond.GetOtherAtomIdx(index)) in selected
                    ),
                },
            )
            for index in selected
        }
    return labels


def _structural_key(
    molecule: Any,
    selected: frozenset[int],
    ring_memberships: Mapping[int, int],
    level: str,
    matching_definition: Mapping[str, Any],
) -> str:
    labels = _refined_labels(
        molecule,
        selected,
        ring_memberships,
        level,
        int(matching_definition["weisfeiler_lehman_iterations"]),
    )
    edges = sorted(
        (
            min(labels[int(bond.GetBeginAtomIdx())], labels[int(bond.GetEndAtomIdx())]),
            max(labels[int(bond.GetBeginAtomIdx())], labels[int(bond.GetEndAtomIdx())]),
            _bond_token(bond),
        )
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in selected
        and int(bond.GetEndAtomIdx()) in selected
    )
    attachments = []
    for bond in _attachment_bonds(molecule, selected):
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        anchor = left if left in selected else right
        external_index = right if anchor == left else left
        external = molecule.GetAtomWithIdx(external_index)
        if level == "exact":
            token = (
                labels[anchor],
                _bond_token(bond),
                external.GetSymbol(),
                int(external.GetFormalCharge()),
                bool(external.GetIsAromatic()),
            )
        elif level == "typed":
            token = (
                labels[anchor],
                _bond_token(bond),
                external.GetSymbol(),
                bool(external.GetIsAromatic()),
            )
        else:
            token = (labels[anchor], _bond_token(bond))
        attachments.append(token)
    payload = {
        "level": level,
        "nodes": sorted(Counter(labels.values()).items()),
        "edges": edges,
        "attachments": sorted(attachments),
        "fragment_smiles": (
            _fragment_smiles(molecule, tuple(selected)) if level == "exact" else None
        ),
        "matching_definition_id": matching_definition["definition_id"],
        "matching_definition_version": matching_definition["definition_version"],
    }
    prefixes = {"exact": "SCEX1", "typed": "SCTY1", "shape": "SCSH1"}
    return _stable_id(prefixes[level], payload)


def _candidate_topology_tokens(
    molecule: Any,
    selected: frozenset[int],
    ring_memberships: Mapping[int, int],
) -> tuple[str, ...]:
    internal_bonds = tuple(
        bond
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in selected
        and int(bond.GetEndAtomIdx()) in selected
    )
    tokens = {
        f"HEAVY_ATOMS:{len(selected)}",
        f"INTERNAL_BONDS:{len(internal_bonds)}",
        f"ATTACHMENTS:{len(_attachment_bonds(molecule, selected))}",
    }
    ring_atoms = sum(ring_memberships.get(index, 0) > 0 for index in selected)
    if ring_atoms:
        tokens.add("CONTAINS_RING_ATOMS")
    if any(ring_memberships.get(index, 0) > 1 for index in selected):
        tokens.add("CONTAINS_MULTI_RING_ATOM")
    if any(
        sum(
            neighbor.GetAtomicNum() > 1 and int(neighbor.GetIdx()) in selected
            for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
        )
        >= 3
        for index in selected
    ):
        tokens.add("CONTAINS_BRANCH_POINT")
    return tuple(sorted(tokens))


def _descriptor_values(
    molecule: Any,
    selected: frozenset[int],
    primary_heavy_count: int,
    ring_systems: Sequence[frozenset[int]],
    stereocenter_indices: frozenset[int],
) -> tuple[tuple[str, float], ...]:
    internal_bonds = [
        bond
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in selected
        and int(bond.GetEndAtomIdx()) in selected
    ]
    components = 1 if selected else 0
    cycle_rank = max(0, len(internal_bonds) - len(selected) + components)
    branch_atoms = sum(
        sum(
            neighbor.GetAtomicNum() > 1 and int(neighbor.GetIdx()) in selected
            for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
        )
        >= 3
        for index in selected
    )
    values = {
        "attachment_count": float(len(_attachment_bonds(molecule, selected))),
        "branch_atom_count": float(branch_atoms),
        "carbon_atom_count": float(
            sum(molecule.GetAtomWithIdx(index).GetAtomicNum() == 6 for index in selected)
        ),
        "cycle_rank": float(cycle_rank),
        "heavy_atom_count": float(len(selected)),
        "heteroatom_count": float(
            sum(molecule.GetAtomWithIdx(index).GetAtomicNum() != 6 for index in selected)
        ),
        "primary_component_fraction": round(
            len(selected) / max(1, primary_heavy_count), 8
        ),
        "ring_atom_count": float(
            sum(any(index in system for system in ring_systems) for index in selected)
        ),
        "ring_system_count": float(
            sum(bool(selected.intersection(system)) for system in ring_systems)
        ),
        "stereocenter_count": float(len(selected.intersection(stereocenter_indices))),
    }
    return tuple(sorted(values.items()))


def _observation(
    molecule: Any,
    molecule_id: str,
    component_index: int,
    candidate: _CandidateCore,
    primary_heavy_count: int,
    ring_systems: Sequence[frozenset[int]],
    ring_memberships: Mapping[int, int],
    stereocenter_indices: frozenset[int],
    observation_definition: Mapping[str, Any],
    matching_definition: Mapping[str, Any],
) -> StructuralCoreObservation:
    selected = candidate.atom_indices
    references = tuple(
        MoleculeAtomReference(
            molecule_id=molecule_id,
            component_index=component_index,
            atom_index=index,
            element=str(molecule.GetAtomWithIdx(index).GetSymbol()),
            formal_charge=int(molecule.GetAtomWithIdx(index).GetFormalCharge()),
            aromatic=bool(molecule.GetAtomWithIdx(index).GetIsAromatic()),
            hybridization=str(molecule.GetAtomWithIdx(index).GetHybridization()),
            environment_id=_atom_environment_id(molecule.GetAtomWithIdx(index)),
        )
        for index in sorted(selected)
    )
    attachment_keys = []
    for bond in _attachment_bonds(molecule, selected):
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        anchor = left if left in selected else right
        external = right if anchor == left else left
        attachment_keys.append(
            f"{component_index}:{anchor}:{external}:{_bond_token(bond)}"
        )
    focus_atom_indices = tuple(
        sorted(index for index in candidate.focus_atom_indices if index in selected)
    )
    focus_bond_keys = []
    for left, right in candidate.focus_bond_pairs:
        bond = molecule.GetBondBetweenAtoms(int(left), int(right))
        if bond is None or left not in selected or right not in selected:
            continue
        focus_bond_keys.append(
            f"{component_index}:{min(left, right)}:{max(left, right)}:{_bond_token(bond)}"
        )
    exact_key = _structural_key(
        molecule, selected, ring_memberships, "exact", matching_definition
    )
    typed_key = _structural_key(
        molecule, selected, ring_memberships, "typed", matching_definition
    )
    shape_key = _structural_key(
        molecule, selected, ring_memberships, "shape", matching_definition
    )
    identity = {
        "molecule_id": molecule_id,
        "component_index": component_index,
        "kind": candidate.kind,
        "atom_indices": sorted(selected),
        "attachment_bond_keys": sorted(attachment_keys),
        "focus_atom_indices": list(focus_atom_indices),
        "focus_bond_keys": sorted(focus_bond_keys),
        "observation_definition_id": observation_definition["definition_id"],
        "observation_definition_version": observation_definition["definition_version"],
        "matching_definition_id": matching_definition["definition_id"],
        "matching_definition_version": matching_definition["definition_version"],
        "schema_version": STRUCTURAL_CORE_SCHEMA_VERSION,
    }
    return StructuralCoreObservation(
        core_observation_id=_stable_id("SCOBS1", identity),
        rank=1,
        molecule_id=molecule_id,
        component_index=component_index,
        kind=candidate.kind,
        atom_references=references,
        attachment_bond_keys=tuple(sorted(attachment_keys)),
        focus_atom_indices=focus_atom_indices,
        focus_bond_keys=tuple(sorted(focus_bond_keys)),
        structural_exact_key=exact_key,
        structural_typed_key=typed_key,
        structural_shape_key=shape_key,
        topology_tokens=_candidate_topology_tokens(
            molecule, selected, ring_memberships
        ),
        descriptor_values=_descriptor_values(
            molecule,
            selected,
            primary_heavy_count,
            ring_systems,
            stereocenter_indices,
        ),
        evidence_codes=tuple(sorted(set(candidate.evidence_codes))),
        warnings=(),
        definition_id=str(observation_definition["definition_id"]),
        definition_version=str(observation_definition["definition_version"]),
        matching_definition_id=str(matching_definition["definition_id"]),
        matching_definition_version=str(matching_definition["definition_version"]),
    )


def observe_structural_cores(smiles: str) -> StructuralCoreAnalysis:
    """Return bounded deterministic structural-core observations for a target."""

    observation_definition = load_structural_core_observation_definition()
    matching_definition = load_structural_core_matching_definition()
    raw = str(smiles or "")
    parsed = parse_smiles(raw)
    canonical = mol_to_canonical_smiles(parsed)
    if parsed is None or not canonical:
        return StructuralCoreAnalysis(
            input_smiles=raw,
            canonical_smiles=canonical,
            molecule_id=None,
            valid=False,
            error="INVALID_TARGET_SMILES",
            definition_version=str(observation_definition["definition_version"]),
            matching_definition_version=str(matching_definition["definition_version"]),
        )
    molecule = parse_smiles(canonical)
    if molecule is None:
        return StructuralCoreAnalysis(
            input_smiles=raw,
            canonical_smiles=canonical,
            molecule_id=None,
            valid=False,
            error="CANONICAL_SMILES_REPARSE_FAILED",
            definition_version=str(observation_definition["definition_version"]),
            matching_definition_version=str(matching_definition["definition_version"]),
        )

    from rdkit import Chem  # type: ignore

    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    molecule_id = _stable_id(
        "MOL1",
        {
            "canonical_smiles": canonical,
            "schema_version": STRUCTURAL_CORE_SCHEMA_VERSION,
        },
    )
    components = _component_records(molecule)
    if not components:
        return StructuralCoreAnalysis(
            input_smiles=raw,
            canonical_smiles=canonical,
            molecule_id=molecule_id,
            valid=False,
            error="TARGET_HAS_NO_COMPONENTS",
            definition_version=str(observation_definition["definition_version"]),
            matching_definition_version=str(matching_definition["definition_version"]),
        )
    primary_component_index, primary_indices, _ = components[0]
    primary = frozenset(
        index
        for index in primary_indices
        if molecule.GetAtomWithIdx(index).GetAtomicNum() > 1
    )
    primary_heavy_count = len(primary)
    warnings: list[str] = []
    if len(components) > 1:
        warnings.append("MULTICOMPONENT_TARGET_PRIMARY_COMPONENT_SELECTED")
        heavy_counts = [
            sum(molecule.GetAtomWithIdx(index).GetAtomicNum() > 1 for index in indices)
            for _, indices, _ in components
        ]
        if len(heavy_counts) > 1 and heavy_counts[0] == heavy_counts[1]:
            warnings.append("AMBIGUOUS_PRIMARY_COMPONENT")
    if primary_heavy_count < int(observation_definition["minimum_candidate_heavy_atoms"]):
        warnings.append("NO_STRUCTURAL_CORE_OBSERVATION")
        return StructuralCoreAnalysis(
            input_smiles=raw,
            canonical_smiles=canonical,
            molecule_id=molecule_id,
            valid=True,
            warnings=tuple(warnings),
            definition_version=str(observation_definition["definition_version"]),
            matching_definition_version=str(matching_definition["definition_version"]),
        )

    atom_rings = tuple(
        frozenset(int(index) for index in ring)
        for ring in molecule.GetRingInfo().AtomRings()
    )
    ring_memberships = {
        index: sum(index in ring for ring in atom_rings)
        for index in primary
    }
    ring_systems = _ring_systems(molecule, primary)
    stereocenters = tuple(
        (int(index), str(assignment))
        for index, assignment in Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
            includeCIP=True,
        )
        if int(index) in primary
    )
    stereocenter_indices = frozenset(index for index, _ in stereocenters)
    candidates: list[_CandidateCore] = []
    scaffold = _murcko_scaffold_atoms(molecule, primary)
    if scaffold:
        candidates.append(
            _CandidateCore(
                "scaffold_backbone",
                scaffold,
                ("BEMIS_MURCKO_RING_LINKER_SCAFFOLD",),
                1.0,
            )
        )
    elif not ring_systems:
        main_skeleton = _acyclic_main_skeleton(molecule, primary)
        candidates.append(
            _CandidateCore(
                "scaffold_backbone",
                main_skeleton,
                ("ACYCLIC_HEAVY_ATOM_MAIN_SKELETON",),
                len(main_skeleton) / max(1, primary_heavy_count),
            )
        )

    for linker, score, evidence, focus_bonds in _direct_ring_linkers(
        molecule,
        ring_systems,
    ):
        candidates.append(
            _CandidateCore(
                "linker_region",
                linker,
                evidence,
                score,
                tuple(sorted(linker)),
                focus_bonds,
            )
        )

    carbon_atoms = frozenset(
        index
        for index in primary
        if molecule.GetAtomWithIdx(index).GetAtomicNum() == 6
    )
    minimum_carbon_atoms = int(
        observation_definition["minimum_carbon_framework_atoms"]
    )
    for framework in _induced_components(molecule, carbon_atoms):
        if len(framework) < minimum_carbon_atoms or framework in ring_systems:
            continue
        branch_count = sum(
            sum(
                int(neighbor.GetIdx()) in framework
                for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()
            )
            >= 3
            for index in framework
        )
        evidence = {"CONTINUOUS_CARBON_FRAMEWORK"}
        if branch_count:
            evidence.add("BRANCHED_CARBON_FRAMEWORK")
        if framework.intersection(stereocenter_indices):
            evidence.add("STEREOCENTER_ON_CARBON_FRAMEWORK")
        candidates.append(
            _CandidateCore(
                "carbon_framework",
                framework,
                tuple(sorted(evidence)),
                len(framework) / max(1, primary_heavy_count) + 0.05 * branch_count,
            )
        )

    all_ring_atoms = frozenset().union(*ring_systems) if ring_systems else frozenset()
    minimum_bridge_side = int(
        observation_definition["minimum_bridge_side_heavy_atoms"]
    )
    minimum_balance = float(
        observation_definition["minimum_bridge_balance_fraction"]
    )
    bridge_radius = int(observation_definition["bridge_neighborhood_radius"])
    for bond in molecule.GetBonds():
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        if (
            left not in primary
            or right not in primary
            or bond.GetIsAromatic()
            or bond.IsInRing()
            or float(bond.GetBondTypeAsDouble()) != 1.0
        ):
            continue
        left_side, right_side = _bond_partition(
            molecule,
            primary,
            left,
            right,
        )
        smaller = min(len(left_side), len(right_side))
        if smaller < minimum_bridge_side:
            continue
        balance = smaller / max(1, primary_heavy_count)
        if balance < minimum_balance:
            continue
        interface = _neighborhood(
            molecule,
            (left, right),
            primary,
            bridge_radius,
        )
        evidence = {"BALANCED_GRAPH_BRIDGE_INTERFACE"}
        if left_side.intersection(all_ring_atoms) and right_side.intersection(
            all_ring_atoms
        ):
            evidence.add("CONNECTS_RING_BEARING_REGIONS")
        endpoint_atomic_numbers = {
            molecule.GetAtomWithIdx(left).GetAtomicNum(),
            molecule.GetAtomWithIdx(right).GetAtomicNum(),
        }
        if endpoint_atomic_numbers == {6}:
            evidence.add("CARBON_CARBON_BRIDGE")
        else:
            evidence.add("HETEROATOM_BRIDGE")
        candidates.append(
            _CandidateCore(
                "bridge_interface",
                interface,
                tuple(sorted(evidence)),
                balance,
                tuple(sorted((left, right))),
                ((min(left, right), max(left, right)),),
            )
        )

    stereo_evidence = {"LOCAL_STEREOCHEMICAL_BACKBONE"}
    if any(assignment == "?" for _, assignment in stereocenters):
        stereo_evidence.add("UNASSIGNED_STEREOCENTER_PRESENT")
    candidates.extend(
        _CandidateCore(
            "stereo_backbone_region",
            _neighborhood(
                molecule,
                (index,),
                primary,
                int(observation_definition["stereo_neighborhood_radius"]),
            ),
            tuple(sorted(stereo_evidence)),
            1.0 if assignment != "?" else 0.5,
            (index,),
        )
        for index, assignment in stereocenters
    )

    for system in ring_systems:
        internal_bonds = [
            bond
            for bond in molecule.GetBonds()
            if int(bond.GetBeginAtomIdx()) in system
            and int(bond.GetEndAtomIdx()) in system
        ]
        cycle_rank = max(0, len(internal_bonds) - len(system) + 1)
        heteroatoms = sum(
            molecule.GetAtomWithIdx(index).GetAtomicNum() != 6 for index in system
        )
        saturated_atoms = sum(
            not molecule.GetAtomWithIdx(index).GetIsAromatic() for index in system
        )
        if not (cycle_rank > 1 or heteroatoms or saturated_atoms):
            continue
        evidence = {"DISTINCTIVE_RING_SYSTEM"}
        if cycle_rank > 1:
            evidence.add("POLYCYCLIC_RING_SYSTEM")
        if heteroatoms:
            evidence.add("HETEROCYCLIC_RING_SYSTEM")
        if saturated_atoms:
            evidence.add("NONAROMATIC_RING_SYSTEM")
        candidates.append(
            _CandidateCore(
                "ring_system",
                system,
                tuple(sorted(evidence)),
                cycle_rank + 0.2 * heteroatoms + 0.05 * saturated_atoms,
            )
        )

    minimum_atoms = int(observation_definition["minimum_candidate_heavy_atoms"])
    priority = {
        str(kind): index
        for index, kind in enumerate(observation_definition["kind_priority"])
    }
    deduplicated: dict[
        tuple[
            frozenset[int],
            tuple[int, ...],
            tuple[tuple[int, int], ...],
        ],
        _CandidateCore,
    ] = {}
    for candidate in candidates:
        selected = frozenset(
            index
            for index in candidate.atom_indices
            if index in primary and molecule.GetAtomWithIdx(index).GetAtomicNum() > 1
        )
        if len(selected) < minimum_atoms:
            continue
        coverage = len(selected) / max(1, primary_heavy_count)
        coverage_limit = float(
            observation_definition["maximum_coverage_by_kind"][candidate.kind]
        )
        if coverage > coverage_limit:
            continue
        normalized = replace(candidate, atom_indices=selected)
        identity = (
            selected,
            tuple(sorted(normalized.focus_atom_indices)),
            tuple(sorted(normalized.focus_bond_pairs)),
        )
        prior = deduplicated.get(identity)
        if prior is None:
            deduplicated[identity] = normalized
            continue
        chosen_kind = min((prior.kind, normalized.kind), key=lambda kind: priority[kind])
        deduplicated[identity] = _CandidateCore(
            chosen_kind,
            selected,
            tuple(sorted({*prior.evidence_codes, *normalized.evidence_codes})),
            max(prior.selection_score, normalized.selection_score),
            tuple(sorted({*prior.focus_atom_indices, *normalized.focus_atom_indices})),
            tuple(sorted({*prior.focus_bond_pairs, *normalized.focus_bond_pairs})),
        )

    candidate_observations = [
        (
            candidate,
            _observation(
            molecule,
            molecule_id,
            primary_component_index,
            candidate,
            primary_heavy_count,
            ring_systems,
            ring_memberships,
            stereocenter_indices,
            observation_definition,
            matching_definition,
            ),
        )
        for candidate in deduplicated.values()
    ]
    candidate_observations.sort(
        key=lambda item: (
            priority[item[1].kind],
            -item[0].selection_score,
            -len(item[1].atom_references),
            item[1].structural_exact_key,
            item[1].atom_indices,
        )
    )
    by_kind: dict[str, list[StructuralCoreObservation]] = {
        kind: [] for kind in observation_definition["kind_priority"]
    }
    for _, observation in candidate_observations:
        by_kind[observation.kind].append(observation)
    selected_observations = [
        values[0]
        for kind in observation_definition["kind_priority"]
        if (values := by_kind[str(kind)])
    ]
    maximum_observations = int(observation_definition["maximum_observations"])
    if len(selected_observations) < maximum_observations:
        already_selected = {item.core_observation_id for item in selected_observations}
        for _, observation in candidate_observations:
            if observation.core_observation_id in already_selected:
                continue
            kind_limit = int(
                observation_definition["maximum_per_kind"][observation.kind]
            )
            if sum(item.kind == observation.kind for item in selected_observations) >= kind_limit:
                continue
            selected_observations.append(observation)
            already_selected.add(observation.core_observation_id)
            if len(selected_observations) >= maximum_observations:
                break
    selected_observations.sort(
        key=lambda item: (
            priority[item.kind],
            next(
                -candidate.selection_score
                for candidate, observation in candidate_observations
                if observation.core_observation_id == item.core_observation_id
            ),
            item.structural_exact_key,
        )
    )
    observations = [
        replace(item, rank=index)
        for index, item in enumerate(
            selected_observations[:maximum_observations],
            start=1,
        )
    ]
    if not observations:
        warnings.append("NO_STRUCTURAL_CORE_OBSERVATION")
    return StructuralCoreAnalysis(
        input_smiles=raw,
        canonical_smiles=canonical,
        molecule_id=molecule_id,
        valid=True,
        observations=tuple(observations),
        warnings=tuple(dict.fromkeys(warnings)),
        definition_version=str(observation_definition["definition_version"]),
        matching_definition_version=str(matching_definition["definition_version"]),
    )


__all__ = [
    "STRUCTURAL_CORE_MATCHING_DEFINITION_ID",
    "STRUCTURAL_CORE_OBSERVATION_DEFINITION_ID",
    "STRUCTURAL_CORE_SCHEMA_VERSION",
    "MoleculeAtomReference",
    "StructuralCoreAnalysis",
    "StructuralCoreKind",
    "StructuralCoreObservation",
    "load_structural_core_matching_definition",
    "load_structural_core_observation_definition",
    "observe_structural_cores",
    "validate_structural_core_matching_definition",
    "validate_structural_core_observation_definition",
]
