"""Immutable public contracts for context-aware reactivity descriptors."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Literal, Tuple, Union


DescriptorStatus = Literal[
    "observed",
    "derived",
    "unresolved",
    "not_applicable",
    "not_computed",
]

AccessibilityClass = Literal["open", "moderate", "hindered", "severe"]
BurdenClass = Literal["none", "low", "medium", "high"]


def _validate_confidence(value: float) -> None:
    if not 0.0 <= float(value) <= 1.0:
        raise ValueError("confidence must be between 0 and 1")


def _validate_score(value: float) -> None:
    if not -1.0 <= float(value) <= 1.0:
        raise ValueError("normalized score must be between -1 and 1")


@dataclass(frozen=True)
class DescriptorEvidence:
    """Method and atom provenance for one derived descriptor."""

    source: str
    method: str
    confidence: float
    contributing_atom_indices: Tuple[int, ...] = ()
    warnings: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        _validate_confidence(self.confidence)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""
        return asdict(self)


@dataclass(frozen=True)
class StericContribution:
    """One bounded graph branch contributing to site shielding."""

    origin_atom_index: int
    relation: str
    heavy_atom_count: int
    branch_count: int
    score: float

    def __post_init__(self) -> None:
        if self.heavy_atom_count < 0 or self.branch_count < 0:
            raise ValueError("steric counts must be non-negative")
        if not 0.0 <= float(self.score) <= 1.0:
            raise ValueError("steric contribution score must be between 0 and 1")


@dataclass(frozen=True)
class ElectronicContribution:
    """One auditable contribution to a local electronic activation axis."""

    source_id: str
    effect: Literal["withdrawing", "donating", "mixed"]
    pathway: Literal[
        "inductive",
        "resonance",
        "charge",
        "aromatic_intrinsic",
    ]
    positional_relation: str
    contribution: float
    atom_indices: Tuple[int, ...] = ()

    def __post_init__(self) -> None:
        _validate_score(self.contribution)


@dataclass(frozen=True)
class ReactiveCenterProfile:
    """Atom state and reactivity-relevant identity of a selected center."""

    atom_index: int
    element: str
    hybridization: str
    formal_charge: int
    radical_electrons: int
    hydrogen_count: int
    heavy_atom_attachment_count: int
    aromatic: bool
    in_ring: bool
    ring_sizes: Tuple[int, ...]
    substitution_class: str | None
    conjugation_class: str | None
    lone_pair_class: str | None
    lone_pair_availability: str | None
    acidity_class: str | None
    evidence: DescriptorEvidence

    def __post_init__(self) -> None:
        if (
            self.atom_index < 0
            or self.radical_electrons < 0
            or self.hydrogen_count < 0
            or self.heavy_atom_attachment_count < 0
        ):
            raise ValueError("reactive-center indices and counts must be non-negative")


@dataclass(frozen=True)
class StericProfile:
    """Context-relative graph accessibility and its branch evidence."""

    accessibility_class: AccessibilityClass
    accessibility_score: float
    approach_burden_class: BurdenClass
    branch_contributions: Tuple[StericContribution, ...]
    context_metrics: Tuple[Tuple[str, str], ...]
    evidence: DescriptorEvidence

    def __post_init__(self) -> None:
        if not 0.0 <= float(self.accessibility_score) <= 1.0:
            raise ValueError("accessibility score must be between 0 and 1")
        if tuple(sorted(self.context_metrics)) != self.context_metrics:
            raise ValueError("context_metrics must be sorted")
        if len({key for key, _ in self.context_metrics}) != len(
            self.context_metrics
        ):
            raise ValueError("context_metrics keys must be unique")

    def metric(self, key: str) -> str | None:
        """Return one normalized context metric by key."""
        return dict(self.context_metrics).get(key)


@dataclass(frozen=True)
class ElectronicProfile:
    """Context-specific activation axis with normalized contribution evidence."""

    activation_axis: str
    activation_class: str
    activation_score: float
    contributions: Tuple[ElectronicContribution, ...]
    evidence: DescriptorEvidence

    def __post_init__(self) -> None:
        _validate_score(self.activation_score)


@dataclass(frozen=True)
class AromaticHeteroatom:
    """One heteroatom positioned relative to an aromatic reactive anchor."""

    atom_index: int
    element: str
    formal_charge: int
    aromatic_role: Literal["pyridine_like", "pyrrole_like", "other"]
    ring_distance_from_anchor: int
    positional_relation: Literal[
        "anchor",
        "ortho",
        "meta",
        "para",
        "remote",
        "fused_other_ring",
    ]
    same_anchor_ring: bool

    def __post_init__(self) -> None:
        if self.atom_index < 0 or self.ring_distance_from_anchor < 0:
            raise ValueError("aromatic heteroatom indices must be non-negative")


@dataclass(frozen=True)
class AromaticContextDescriptor:
    """Ring-system topology and position-normalized aromatic identity."""

    context_kind: Literal["aromatic"]
    system_class: Literal["carbocyclic", "heteroaromatic", "mixed_fused"]
    ring_family: str
    ring_sizes: Tuple[int, ...]
    aromatic_ring_count: int
    fused: bool
    anchor_in_ring: bool
    heteroatoms: Tuple[AromaticHeteroatom, ...]
    ortho_occupancy_count: int
    ortho_capacity: int
    ortho_burden_class: BurdenClass
    ortho_burden_score: float

    def __post_init__(self) -> None:
        if (
            self.aromatic_ring_count < 1
            or self.ortho_occupancy_count < 0
            or self.ortho_capacity < 0
        ):
            raise ValueError("aromatic counts must be non-negative")
        if self.ortho_occupancy_count > self.ortho_capacity:
            raise ValueError("ortho occupancy cannot exceed capacity")
        if not 0.0 <= float(self.ortho_burden_score) <= 1.0:
            raise ValueError("ortho burden score must be between 0 and 1")


@dataclass(frozen=True)
class AlkylContextDescriptor:
    """Substitution, branching, activation, and cyclicity of an sp3 carbon."""

    context_kind: Literal["alkyl"]
    carbon_substitution: Literal["methyl", "primary", "secondary", "tertiary"]
    alpha_carbon_neighbor_count: int
    alpha_branched: bool
    beta_branch_count: int
    beta_hydrogen_count: int
    cyclic: bool
    ring_sizes: Tuple[int, ...]
    benzylic: bool
    allylic: bool
    propargylic: bool
    adjacent_heteroatoms: Tuple[str, ...]


@dataclass(frozen=True)
class AlkenylContextDescriptor:
    """Local topology of a non-aromatic carbon-carbon double bond."""

    context_kind: Literal["alkenyl"]
    endpoint_substitution: Tuple[int, int]
    alkene_class: str
    stereochemistry: str | None
    cyclic: bool
    ring_size: int | None
    conjugation_class: str
    allylic_branch_count: int


@dataclass(frozen=True)
class AlkynylContextDescriptor:
    """Local topology of a non-aromatic carbon-carbon triple bond."""

    context_kind: Literal["alkynyl"]
    terminal: bool
    endpoint_substitution: Tuple[int, int]
    conjugation_class: str
    propargylic_branch_count: int


@dataclass(frozen=True)
class ActivatedCenterContextDescriptor:
    """Acyl, sulfonyl, or phosphoryl center topology."""

    context_kind: Literal["acyl", "sulfonyl", "phosphoryl"]
    center_class: str
    attached_group_classes: Tuple[str, ...]
    conjugation_class: str
    heteroatom_substitution: Tuple[str, ...]
    enolizable: bool | None


@dataclass(frozen=True)
class AttachedGroupProfile:
    """One heavy group directly attached to a reactive heteroatom."""

    atom_index: int
    context: str
    element: str
    attachment_carbon_class: str | None
    alpha_branched: bool
    beta_branch_count: int


@dataclass(frozen=True)
class HeteroatomContextDescriptor:
    """Attachment and lone-pair context for a reactive heteroatom."""

    context_kind: Literal["heteroatom"]
    element: str
    substitution_class: str
    attached_contexts: Tuple[str, ...]
    resonance_class: str
    lone_pair_class: str
    proton_count: int
    alpha_branched_group_count: int
    attached_groups: Tuple[AttachedGroupProfile, ...]


@dataclass(frozen=True)
class OtherContextDescriptor:
    """Conservative fallback for a center outside supported context variants."""

    context_kind: Literal["other"]
    center_element: str
    reason: str


ContextDescriptor = Union[
    AromaticContextDescriptor,
    AlkylContextDescriptor,
    AlkenylContextDescriptor,
    AlkynylContextDescriptor,
    ActivatedCenterContextDescriptor,
    HeteroatomContextDescriptor,
    OtherContextDescriptor,
]


@dataclass(frozen=True)
class ReactivityModifier:
    """Orthogonal handle or liability attached to a scaffold context."""

    modifier_type: Literal[
        "leaving_group",
        "transfer_carrier",
        "removable_hydrogen",
        "coordination",
        "redox",
        "elimination",
        "strain",
        "other",
    ]
    modifier_id: str
    class_name: str
    attributes: Tuple[Tuple[str, str], ...]
    evidence: DescriptorEvidence

    def __post_init__(self) -> None:
        if tuple(sorted(self.attributes)) != self.attributes:
            raise ValueError("modifier attributes must be sorted")


@dataclass(frozen=True)
class SiteReactivityProfile:
    """Complete typed reactivity observation for one molecular site."""

    site_id: str
    center_atom_index: int
    context_kind: str
    context: ContextDescriptor
    reactive_center: ReactiveCenterProfile
    steric: StericProfile
    electronic: ElectronicProfile
    modifiers: Tuple[ReactivityModifier, ...]
    flags: Tuple[str, ...]
    status: DescriptorStatus
    definition_versions: Tuple[Tuple[str, str], ...]
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if self.center_atom_index != self.reactive_center.atom_index:
            raise ValueError("profile center must match reactive-center atom")
        if self.context_kind != self.context.context_kind:
            raise ValueError("profile context_kind must match context descriptor")
        if tuple(sorted(self.flags)) != self.flags:
            raise ValueError("profile flags must be sorted")
        if tuple(sorted(self.definition_versions)) != self.definition_versions:
            raise ValueError("definition_versions must be sorted")

    def to_dict(self) -> Dict[str, Any]:
        """Return a canonical JSON-serializable representation."""
        return asdict(self)


__all__ = [
    "ActivatedCenterContextDescriptor",
    "AlkylContextDescriptor",
    "AlkenylContextDescriptor",
    "AlkynylContextDescriptor",
    "AromaticContextDescriptor",
    "AromaticHeteroatom",
    "AttachedGroupProfile",
    "ContextDescriptor",
    "DescriptorEvidence",
    "DescriptorStatus",
    "ElectronicContribution",
    "ElectronicProfile",
    "HeteroatomContextDescriptor",
    "OtherContextDescriptor",
    "ReactiveCenterProfile",
    "ReactivityModifier",
    "SiteReactivityProfile",
    "StericContribution",
    "StericProfile",
]
