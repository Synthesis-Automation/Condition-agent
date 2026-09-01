"""Role-neutral contracts for target partitions and their evidence.

The contracts in this module describe target-derived atom ownership.  They do
not claim that a partition is an executable reaction or a feasible route.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable, Mapping, Sequence

from rdkit import Chem

from .chemistry import canonical_smiles, digest


SYNTHETIC_PARTITION_SCHEMA_VERSION = "1.0"
SYNTHETIC_PARTITION_DEFINITION_ID = "synthetic_partition_policy.v1"
PARTITION_EVIDENCE_LEVELS = frozenset({"E0", "E1", "E2", "E3", "E4"})
PARTITION_SOURCE_KINDS = frozenset(
    {
        "structural_baseline",
        "validated_operator_projection",
        "operator_combination_unrealized",
        "observed_route_frontier",
        "planned_route_frontier",
    }
)
PARTITION_REALIZATION_STATUSES = frozenset(
    {
        "fully_realized",
        "partially_realized",
        "unrealized_but_plausible",
        "contradicted",
        "not_attempted",
    }
)


def _positive_sorted_unique(values: Iterable[int], field: str) -> tuple[int, ...]:
    normalized = tuple(sorted({int(value) for value in values}))
    if not normalized or normalized[0] <= 0:
        raise ValueError(f"{field} must contain positive atom-map numbers")
    return normalized


def _sorted_unique_strings(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({str(value) for value in values if str(value)}))


@dataclass(frozen=True)
class TargetAtomReference:
    """Stable target-side atom reference used by every partition."""

    atom_map: int
    atom_index: int
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    local_environment_id: str

    def __post_init__(self) -> None:
        if self.atom_map <= 0 or self.atom_index < 0 or not self.element:
            raise ValueError("invalid target atom reference")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible atom reference."""

        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "TargetAtomReference":
        """Reconstruct a target atom reference."""

        return cls(
            atom_map=int(value["atom_map"]),
            atom_index=int(value["atom_index"]),
            element=str(value["element"]),
            formal_charge=int(value["formal_charge"]),
            aromatic=bool(value["aromatic"]),
            hybridization=str(value["hybridization"]),
            local_environment_id=str(value["local_environment_id"]),
        )


@dataclass(frozen=True)
class TargetModule:
    """One role-neutral, nonempty set of target-derived atoms."""

    module_id: str
    target_id: str
    target_atom_maps: tuple[int, ...]
    attachment_atom_maps: tuple[int, ...] = ()
    graph_complexity: float | None = None

    def __post_init__(self) -> None:
        normalized = _positive_sorted_unique(
            self.target_atom_maps,
            "target_atom_maps",
        )
        if normalized != self.target_atom_maps:
            raise ValueError("target module atom maps must be sorted and unique")
        attachments = tuple(sorted(set(self.attachment_atom_maps)))
        if attachments != self.attachment_atom_maps:
            raise ValueError("attachment atom maps must be sorted and unique")
        if not set(attachments).issubset(normalized):
            raise ValueError("attachment atoms must belong to the module")
        if not self.module_id.startswith("SPMOD1:") or not self.target_id:
            raise ValueError("target module requires content-derived identities")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible target module."""

        value = asdict(self)
        value["target_atom_maps"] = list(self.target_atom_maps)
        value["attachment_atom_maps"] = list(self.attachment_atom_maps)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "TargetModule":
        """Reconstruct a target module."""

        return cls(
            module_id=str(value["module_id"]),
            target_id=str(value["target_id"]),
            target_atom_maps=tuple(
                int(item) for item in value.get("target_atom_maps") or ()
            ),
            attachment_atom_maps=tuple(
                int(item) for item in value.get("attachment_atom_maps") or ()
            ),
            graph_complexity=(
                float(value["graph_complexity"])
                if value.get("graph_complexity") is not None
                else None
            ),
        )


@dataclass(frozen=True)
class ModuleAnnotation:
    """Optional interpretation that never participates in module identity."""

    module_id: str
    label: str
    proposed_role: str
    confidence: float | None = None
    evidence: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.confidence is not None and not 0.0 <= self.confidence <= 1.0:
            raise ValueError("module annotation confidence must be in [0, 1]")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible annotation."""

        value = asdict(self)
        value["evidence"] = list(self.evidence)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ModuleAnnotation":
        """Reconstruct an optional module annotation."""

        return cls(
            module_id=str(value["module_id"]),
            label=str(value.get("label") or ""),
            proposed_role=str(value.get("proposed_role") or "unresolved"),
            confidence=(
                float(value["confidence"])
                if value.get("confidence") is not None
                else None
            ),
            evidence=tuple(str(item) for item in value.get("evidence") or ()),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class NonTargetAtom:
    """One precursor atom that is not owned by the target partition."""

    element: str
    atom_map: int | None
    classification: str


@dataclass(frozen=True)
class LatentModuleState:
    """A concrete precursor state contributing one or more target modules."""

    latent_state_id: str
    module_ids: tuple[str, ...]
    mapped_smiles: str
    target_atom_maps: tuple[int, ...]
    non_target_atoms: tuple[NonTargetAtom, ...] = ()
    state_annotations: tuple[str, ...] = ()
    evidence_status: str = "unresolved"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible latent state."""

        return {
            "latent_state_id": self.latent_state_id,
            "module_ids": list(self.module_ids),
            "mapped_smiles": self.mapped_smiles,
            "target_atom_maps": list(self.target_atom_maps),
            "non_target_atoms": [asdict(atom) for atom in self.non_target_atoms],
            "state_annotations": list(self.state_annotations),
            "evidence_status": self.evidence_status,
        }


@dataclass(frozen=True)
class InterfaceHypothesis:
    """Unmaterialized interface evidence used while constructing a partition."""

    interface_kind: str
    target_bonds: tuple[tuple[int, int, str], ...]
    candidate_operator_ids: tuple[str, ...] = ()
    candidate_strategy_ids: tuple[str, ...] = ()
    precedent_reaction_ids: tuple[str, ...] = ()
    evidence_level: str = "E0"
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.evidence_level not in PARTITION_EVIDENCE_LEVELS:
            raise ValueError("unsupported interface evidence level")


@dataclass(frozen=True)
class StrategicInterface:
    """Evidence for a target-side relationship between partition modules."""

    interface_id: str
    interface_kind: str
    module_ids: tuple[str, ...]
    target_bonds: tuple[tuple[int, int, str], ...]
    candidate_operator_ids: tuple[str, ...]
    candidate_strategy_ids: tuple[str, ...]
    precedent_reaction_ids: tuple[str, ...]
    evidence_level: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.interface_id.startswith("SPINT1:"):
            raise ValueError("strategic interface requires a content-derived ID")
        if not self.module_ids:
            raise ValueError("strategic interface requires participating modules")
        if self.evidence_level not in PARTITION_EVIDENCE_LEVELS:
            raise ValueError("unsupported interface evidence level")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible strategic interface."""

        return {
            "interface_id": self.interface_id,
            "interface_kind": self.interface_kind,
            "module_ids": list(self.module_ids),
            "target_bonds": [list(bond) for bond in self.target_bonds],
            "candidate_operator_ids": list(self.candidate_operator_ids),
            "candidate_strategy_ids": list(self.candidate_strategy_ids),
            "precedent_reaction_ids": list(self.precedent_reaction_ids),
            "evidence_level": self.evidence_level,
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StrategicInterface":
        """Reconstruct a strategic interface."""

        return cls(
            interface_id=str(value["interface_id"]),
            interface_kind=str(value["interface_kind"]),
            module_ids=tuple(str(item) for item in value.get("module_ids") or ()),
            target_bonds=tuple(
                (int(item[0]), int(item[1]), str(item[2]))
                for item in value.get("target_bonds") or ()
            ),
            candidate_operator_ids=tuple(
                str(item) for item in value.get("candidate_operator_ids") or ()
            ),
            candidate_strategy_ids=tuple(
                str(item) for item in value.get("candidate_strategy_ids") or ()
            ),
            precedent_reaction_ids=tuple(
                str(item) for item in value.get("precedent_reaction_ids") or ()
            ),
            evidence_level=str(value["evidence_level"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class SyntheticPartition:
    """One unordered partition of all target-derived heavy atoms."""

    partition_id: str
    target_id: str
    target_smiles: str
    target_atom_maps: tuple[int, ...]
    modules: tuple[TargetModule, ...]
    interfaces: tuple[StrategicInterface, ...]
    annotations: tuple[ModuleAnnotation, ...]
    generation_evidence: tuple[str, ...]
    evidence_level: str
    source_kind: str
    heuristic_score: float = 0.0
    realization_status: str = "not_attempted"
    warnings: tuple[str, ...] = ()
    definition_id: str = SYNTHETIC_PARTITION_DEFINITION_ID
    schema_version: str = SYNTHETIC_PARTITION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        report = validate_synthetic_partition(self)
        if not report.valid:
            raise ValueError(
                "invalid synthetic partition: " + ", ".join(report.issues)
            )

    @property
    def k(self) -> int:
        """Return the number of role-neutral target modules."""

        return len(self.modules)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible partition."""

        return {
            "partition_id": self.partition_id,
            "target_id": self.target_id,
            "target_smiles": self.target_smiles,
            "target_atom_maps": list(self.target_atom_maps),
            "k": self.k,
            "modules": [module.to_dict() for module in self.modules],
            "interfaces": [item.to_dict() for item in self.interfaces],
            "annotations": [item.to_dict() for item in self.annotations],
            "generation_evidence": list(self.generation_evidence),
            "evidence_level": self.evidence_level,
            "source_kind": self.source_kind,
            "heuristic_score": self.heuristic_score,
            "realization_status": self.realization_status,
            "warnings": list(self.warnings),
            "definition_id": self.definition_id,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SyntheticPartition":
        """Reconstruct and validate a serialized partition."""

        return cls(
            partition_id=str(value["partition_id"]),
            target_id=str(value["target_id"]),
            target_smiles=str(value["target_smiles"]),
            target_atom_maps=tuple(
                int(item) for item in value.get("target_atom_maps") or ()
            ),
            modules=tuple(
                TargetModule.from_dict(item) for item in value.get("modules") or ()
            ),
            interfaces=tuple(
                StrategicInterface.from_dict(item)
                for item in value.get("interfaces") or ()
            ),
            annotations=tuple(
                ModuleAnnotation.from_dict(item)
                for item in value.get("annotations") or ()
            ),
            generation_evidence=tuple(
                str(item) for item in value.get("generation_evidence") or ()
            ),
            evidence_level=str(value["evidence_level"]),
            source_kind=str(value["source_kind"]),
            heuristic_score=float(value.get("heuristic_score") or 0.0),
            realization_status=str(
                value.get("realization_status") or "not_attempted"
            ),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            definition_id=str(
                value.get("definition_id") or SYNTHETIC_PARTITION_DEFINITION_ID
            ),
            schema_version=str(
                value.get("schema_version") or SYNTHETIC_PARTITION_SCHEMA_VERSION
            ),
        )


@dataclass(frozen=True)
class PartitionValidationReport:
    """Deterministic validation result for one target partition."""

    valid: bool
    issues: tuple[str, ...]


@dataclass(frozen=True)
class PartitionSearchDiagnostics:
    """Review-only counters for partition generation."""

    supplied_candidate_count: int = 0
    accepted_seed_count: int = 0
    rejected_candidate_count: int = 0
    generated_combination_count: int = 0
    duplicate_partition_count: int = 0
    returned_partition_count: int = 0
    rejection_counts: tuple[tuple[str, int], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible search diagnostics."""

        value = asdict(self)
        value["rejection_counts"] = dict(self.rejection_counts)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionSearchDiagnostics":
        """Reconstruct partition search diagnostics."""

        raw_rejections = value.get("rejection_counts") or {}
        if not isinstance(raw_rejections, Mapping):
            raise ValueError("partition rejection counts must be an object")
        return cls(
            supplied_candidate_count=int(
                value.get("supplied_candidate_count") or 0
            ),
            accepted_seed_count=int(value.get("accepted_seed_count") or 0),
            rejected_candidate_count=int(
                value.get("rejected_candidate_count") or 0
            ),
            generated_combination_count=int(
                value.get("generated_combination_count") or 0
            ),
            duplicate_partition_count=int(
                value.get("duplicate_partition_count") or 0
            ),
            returned_partition_count=int(
                value.get("returned_partition_count") or 0
            ),
            rejection_counts=tuple(
                sorted(
                    (str(key), int(count))
                    for key, count in raw_rejections.items()
                )
            ),
        )


@dataclass(frozen=True)
class SyntheticPartitionLandscape:
    """A bounded, role-neutral collection of target partition views."""

    target_id: str
    target_smiles: str
    target_atoms: tuple[TargetAtomReference, ...]
    partitions: tuple[SyntheticPartition, ...]
    searched_k_values: tuple[int, ...]
    generation_coverage: tuple[str, ...]
    unresolved_motifs: tuple[str, ...]
    abstained: bool
    abstention_reasons: tuple[str, ...]
    diagnostics: PartitionSearchDiagnostics
    warnings: tuple[str, ...] = ()
    definition_id: str = SYNTHETIC_PARTITION_DEFINITION_ID
    schema_version: str = SYNTHETIC_PARTITION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.target_atoms or not self.partitions:
            raise ValueError("partition landscape requires target atoms and views")
        if tuple(sorted(set(self.searched_k_values))) != self.searched_k_values:
            raise ValueError("searched k values must be sorted and unique")
        if len({item.partition_id for item in self.partitions}) != len(
            self.partitions
        ):
            raise ValueError("partition landscape contains duplicate identities")
        if any(item.target_id != self.target_id for item in self.partitions):
            raise ValueError("partition landscape mixes targets")
        if build_target_id(self.target_smiles) != self.target_id:
            raise ValueError("partition landscape target identity mismatch")
        if tuple(atom.atom_map for atom in self.target_atoms) != tuple(
            range(1, len(self.target_atoms) + 1)
        ):
            raise ValueError("partition landscape target maps are not canonical")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible partition landscape."""

        return {
            "target_id": self.target_id,
            "target_smiles": self.target_smiles,
            "target_atoms": [atom.to_dict() for atom in self.target_atoms],
            "partitions": [partition.to_dict() for partition in self.partitions],
            "searched_k_values": list(self.searched_k_values),
            "generation_coverage": list(self.generation_coverage),
            "unresolved_motifs": list(self.unresolved_motifs),
            "abstained": self.abstained,
            "abstention_reasons": list(self.abstention_reasons),
            "diagnostics": self.diagnostics.to_dict(),
            "warnings": list(self.warnings),
            "definition_id": self.definition_id,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SyntheticPartitionLandscape":
        """Reconstruct and validate a serialized partition landscape."""

        diagnostics = value.get("diagnostics")
        if not isinstance(diagnostics, Mapping):
            raise ValueError("partition landscape diagnostics must be an object")
        return cls(
            target_id=str(value["target_id"]),
            target_smiles=str(value["target_smiles"]),
            target_atoms=tuple(
                TargetAtomReference.from_dict(item)
                for item in value.get("target_atoms") or ()
            ),
            partitions=tuple(
                SyntheticPartition.from_dict(item)
                for item in value.get("partitions") or ()
            ),
            searched_k_values=tuple(
                int(item) for item in value.get("searched_k_values") or ()
            ),
            generation_coverage=tuple(
                str(item) for item in value.get("generation_coverage") or ()
            ),
            unresolved_motifs=tuple(
                str(item) for item in value.get("unresolved_motifs") or ()
            ),
            abstained=bool(value.get("abstained")),
            abstention_reasons=tuple(
                str(item) for item in value.get("abstention_reasons") or ()
            ),
            diagnostics=PartitionSearchDiagnostics.from_dict(diagnostics),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            definition_id=str(
                value.get("definition_id") or SYNTHETIC_PARTITION_DEFINITION_ID
            ),
            schema_version=str(
                value.get("schema_version") or SYNTHETIC_PARTITION_SCHEMA_VERSION
            ),
        )


def build_target_id(target_smiles: str) -> str:
    """Return the content-derived identity of one canonical target."""

    canonical = canonical_smiles(target_smiles)
    if canonical is None or "." in canonical:
        raise ValueError("target must be one valid molecule")
    return digest("SPTGT1", canonical, SYNTHETIC_PARTITION_SCHEMA_VERSION)


def build_module_id(target_id: str, atom_maps: Iterable[int]) -> str:
    """Return an order-independent target-module identity."""

    normalized = _positive_sorted_unique(atom_maps, "module atom maps")
    return digest("SPMOD1", target_id, ",".join(map(str, normalized)))


def build_partition_id(
    target_id: str,
    module_atom_sets: Iterable[Iterable[int]],
) -> str:
    """Return a partition identity invariant to module permutation."""

    blocks = tuple(
        sorted(
            ",".join(map(str, _positive_sorted_unique(block, "module atom maps")))
            for block in module_atom_sets
        )
    )
    if len(set(blocks)) != len(blocks):
        raise ValueError("partition contains duplicate modules")
    return digest(
        "SPART1",
        target_id,
        *blocks,
        SYNTHETIC_PARTITION_SCHEMA_VERSION,
        SYNTHETIC_PARTITION_DEFINITION_ID,
    )


def analyze_partition_target(
    target_smiles: str,
) -> tuple[str, str, tuple[TargetAtomReference, ...], dict[int, int]]:
    """Canonicalize a target and assign deterministic heavy-atom references."""

    canonical = canonical_smiles(target_smiles)
    if canonical is None or "." in canonical:
        raise ValueError("target must be one valid molecule")
    molecule = Chem.MolFromSmiles(canonical)
    if molecule is None:
        raise ValueError("target could not be parsed")
    target_id = build_target_id(canonical)
    atom_index_to_map: dict[int, int] = {}
    references: list[TargetAtomReference] = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        atom_map = len(references) + 1
        atom_index_to_map[atom.GetIdx()] = atom_map
        neighbors = tuple(
            sorted(
                (
                    neighbor.GetSymbol(),
                    str(molecule.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()),
                )
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() > 1
            )
        )
        local_environment_id = digest(
            "SPENV1",
            atom.GetSymbol(),
            str(atom.GetFormalCharge()),
            str(atom.GetIsAromatic()),
            str(atom.GetHybridization()),
            repr(neighbors),
        )
        references.append(
            TargetAtomReference(
                atom_map=atom_map,
                atom_index=atom.GetIdx(),
                element=atom.GetSymbol(),
                formal_charge=atom.GetFormalCharge(),
                aromatic=atom.GetIsAromatic(),
                hybridization=str(atom.GetHybridization()),
                local_environment_id=local_environment_id,
            )
        )
    if not references:
        raise ValueError("target must contain at least one heavy atom")
    return canonical, target_id, tuple(references), atom_index_to_map


def _target_bonds(
    target_smiles: str,
    atom_index_to_map: Mapping[int, int],
) -> tuple[tuple[int, int, str], ...]:
    molecule = Chem.MolFromSmiles(target_smiles)
    if molecule is None:
        raise ValueError("target could not be parsed")
    values = []
    for bond in molecule.GetBonds():
        left = atom_index_to_map.get(bond.GetBeginAtomIdx())
        right = atom_index_to_map.get(bond.GetEndAtomIdx())
        if left is None or right is None:
            continue
        values.append((min(left, right), max(left, right), str(bond.GetBondType())))
    return tuple(sorted(values))


def create_synthetic_partition(
    target_smiles: str,
    module_atom_sets: Iterable[Iterable[int]],
    *,
    source_kind: str,
    evidence_level: str,
    interface_hypotheses: Sequence[InterfaceHypothesis] = (),
    generation_evidence: Iterable[str] = (),
    annotations: Sequence[ModuleAnnotation] = (),
    heuristic_score: float = 0.0,
    warnings: Iterable[str] = (),
) -> SyntheticPartition:
    """Build and validate one canonical role-neutral partition."""

    canonical, target_id, references, index_to_map = analyze_partition_target(
        target_smiles
    )
    if source_kind not in PARTITION_SOURCE_KINDS:
        raise ValueError("unsupported partition source kind")
    if evidence_level not in PARTITION_EVIDENCE_LEVELS:
        raise ValueError("unsupported partition evidence level")
    blocks = tuple(
        sorted(
            (
                _positive_sorted_unique(block, "module atom maps")
                for block in module_atom_sets
            ),
            key=lambda item: (item[0], len(item), item),
        )
    )
    target_maps = tuple(reference.atom_map for reference in references)
    owner: dict[int, int] = {}
    for module_index, block in enumerate(blocks):
        for atom_map in block:
            if atom_map in owner:
                raise ValueError("partition duplicates target atom ownership")
            owner[atom_map] = module_index
    if set(owner) != set(target_maps):
        raise ValueError("partition must cover every target heavy atom exactly once")
    all_target_bonds = _target_bonds(canonical, index_to_map)
    attachment_maps: dict[int, set[int]] = {
        index: set() for index in range(len(blocks))
    }
    for left, right, _ in all_target_bonds:
        left_owner = owner[left]
        right_owner = owner[right]
        if left_owner != right_owner:
            attachment_maps[left_owner].add(left)
            attachment_maps[right_owner].add(right)
    modules = tuple(
        TargetModule(
            module_id=build_module_id(target_id, block),
            target_id=target_id,
            target_atom_maps=block,
            attachment_atom_maps=tuple(sorted(attachment_maps[index])),
        )
        for index, block in enumerate(blocks)
    )
    partition_id = build_partition_id(
        target_id,
        (module.target_atom_maps for module in modules),
    )
    module_by_map = {
        atom_map: module.module_id
        for module in modules
        for atom_map in module.target_atom_maps
    }
    interfaces = []
    for hypothesis in interface_hypotheses:
        bonds = tuple(
            sorted(
                {
                    (min(int(left), int(right)), max(int(left), int(right)), str(kind))
                    for left, right, kind in hypothesis.target_bonds
                }
            )
        )
        participating = {
            module_by_map[atom_map]
            for left, right, _ in bonds
            for atom_map in (left, right)
            if atom_map in module_by_map
        }
        if not participating and len(modules) == 1:
            participating.add(modules[0].module_id)
        if not participating:
            participating.update(module.module_id for module in modules)
        module_ids = tuple(sorted(participating))
        encoding = tuple(f"{left}-{right}-{kind}" for left, right, kind in bonds)
        interface_id = digest(
            "SPINT1",
            partition_id,
            hypothesis.interface_kind,
            *module_ids,
            *encoding,
        )
        interfaces.append(
            StrategicInterface(
                interface_id=interface_id,
                interface_kind=hypothesis.interface_kind,
                module_ids=module_ids,
                target_bonds=bonds,
                candidate_operator_ids=_sorted_unique_strings(
                    hypothesis.candidate_operator_ids
                ),
                candidate_strategy_ids=_sorted_unique_strings(
                    hypothesis.candidate_strategy_ids
                ),
                precedent_reaction_ids=_sorted_unique_strings(
                    hypothesis.precedent_reaction_ids
                ),
                evidence_level=hypothesis.evidence_level,
                warnings=_sorted_unique_strings(hypothesis.warnings),
            )
        )
    ordered_interfaces = tuple(
        sorted(
            {item.interface_id: item for item in interfaces}.values(),
            key=lambda item: item.interface_id,
        )
    )
    module_ids = {module.module_id for module in modules}
    ordered_annotations = tuple(
        sorted(
            (item for item in annotations if item.module_id in module_ids),
            key=lambda item: (item.module_id, item.proposed_role, item.label),
        )
    )
    return SyntheticPartition(
        partition_id=partition_id,
        target_id=target_id,
        target_smiles=canonical,
        target_atom_maps=target_maps,
        modules=modules,
        interfaces=ordered_interfaces,
        annotations=ordered_annotations,
        generation_evidence=_sorted_unique_strings(generation_evidence),
        evidence_level=evidence_level,
        source_kind=source_kind,
        heuristic_score=round(float(heuristic_score), 8),
        warnings=_sorted_unique_strings(warnings),
    )


def validate_synthetic_partition(
    partition: SyntheticPartition,
) -> PartitionValidationReport:
    """Return stable contract violations without modifying the partition."""

    issues: list[str] = []
    if partition.schema_version != SYNTHETIC_PARTITION_SCHEMA_VERSION:
        issues.append("invalid_schema_version")
    if partition.definition_id != SYNTHETIC_PARTITION_DEFINITION_ID:
        issues.append("invalid_definition_id")
    if partition.evidence_level not in PARTITION_EVIDENCE_LEVELS:
        issues.append("invalid_evidence_level")
    if partition.source_kind not in PARTITION_SOURCE_KINDS:
        issues.append("invalid_source_kind")
    if partition.realization_status not in PARTITION_REALIZATION_STATUSES:
        issues.append("invalid_realization_status")
    if not partition.partition_id.startswith("SPART1:"):
        issues.append("invalid_partition_id")
    canonical_target = canonical_smiles(partition.target_smiles)
    if canonical_target != partition.target_smiles:
        issues.append("noncanonical_target_smiles")
    try:
        if build_target_id(partition.target_smiles) != partition.target_id:
            issues.append("target_identity_mismatch")
    except ValueError:
        issues.append("invalid_target")
    if tuple(sorted(set(partition.target_atom_maps))) != partition.target_atom_maps:
        issues.append("invalid_target_atom_maps")
    module_ids = tuple(module.module_id for module in partition.modules)
    if len(set(module_ids)) != len(module_ids):
        issues.append("duplicate_module_id")
    if any(module.target_id != partition.target_id for module in partition.modules):
        issues.append("module_target_identity_mismatch")
    if any(
        module.module_id
        != build_module_id(partition.target_id, module.target_atom_maps)
        for module in partition.modules
    ):
        issues.append("module_identity_mismatch")
    if tuple(
        sorted(
            partition.modules,
            key=lambda module: (
                module.target_atom_maps[0],
                len(module.target_atom_maps),
                module.target_atom_maps,
            ),
        )
    ) != partition.modules:
        issues.append("noncanonical_module_order")
    owned = [
        atom_map
        for module in partition.modules
        for atom_map in module.target_atom_maps
    ]
    if len(owned) != len(set(owned)):
        issues.append("duplicate_atom_ownership")
    if set(owned) != set(partition.target_atom_maps):
        issues.append("incomplete_atom_coverage")
    expected_partition_id = None
    try:
        expected_partition_id = build_partition_id(
            partition.target_id,
            (module.target_atom_maps for module in partition.modules),
        )
    except ValueError:
        issues.append("invalid_module_partition")
    if expected_partition_id and expected_partition_id != partition.partition_id:
        issues.append("partition_identity_mismatch")
    known_modules = set(module_ids)
    if any(
        not set(interface.module_ids).issubset(known_modules)
        for interface in partition.interfaces
    ):
        issues.append("interface_unknown_module")
    if any(annotation.module_id not in known_modules for annotation in partition.annotations):
        issues.append("annotation_unknown_module")
    return PartitionValidationReport(valid=not issues, issues=tuple(issues))


__all__ = [
    "InterfaceHypothesis",
    "LatentModuleState",
    "ModuleAnnotation",
    "NonTargetAtom",
    "PARTITION_EVIDENCE_LEVELS",
    "PARTITION_REALIZATION_STATUSES",
    "PARTITION_SOURCE_KINDS",
    "PartitionSearchDiagnostics",
    "PartitionValidationReport",
    "SYNTHETIC_PARTITION_DEFINITION_ID",
    "SYNTHETIC_PARTITION_SCHEMA_VERSION",
    "StrategicInterface",
    "SyntheticPartition",
    "SyntheticPartitionLandscape",
    "TargetAtomReference",
    "TargetModule",
    "analyze_partition_target",
    "build_module_id",
    "build_partition_id",
    "build_target_id",
    "create_synthetic_partition",
    "validate_synthetic_partition",
]
