"""Typed molecular structure and optional interpretation contracts."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple

from .descriptors.models import SiteReactivityProfile


SiteType = Literal[
    "leaving_group",
    "pronucleophile_XH",
    "nucleophile_anion",
    "transfer_group",
    "addition_donor",
    "eliminable_pair",
    "electrophilic_center",
    "aromatic_CH",
    "unsaturated_bond",
    "dipolar_group",
    "heteroatom_bond",
]
SiteTopology = Literal["edge", "atom", "center", "bond"]


@dataclass(frozen=True)
class MolecularAtomObservation:
    """One atom as represented in the parsed molecular graph."""

    atom_index: int
    element: str
    isotope: int
    formal_charge: int
    aromatic: bool
    hybridization: str
    total_hydrogens: int
    degree: int
    in_ring: bool
    atom_map_number: Optional[int]


@dataclass(frozen=True)
class MolecularBondObservation:
    """One bond as represented in the parsed molecular graph."""

    bond_index: int
    atom_1_index: int
    atom_2_index: int
    order: str
    aromatic: bool
    conjugated: bool
    in_ring: bool
    stereo: str


@dataclass(frozen=True)
class MolecularComponentStructure:
    """Structure-only observation for one disconnected component."""

    component_index: int
    input_smiles: str
    canonical_smiles: str
    atom_offset: int
    atoms: Tuple[MolecularAtomObservation, ...]
    bonds: Tuple[MolecularBondObservation, ...]
    schema_version: str = "1.0"


@dataclass(frozen=True)
class MolecularStructureObservation:
    """Parsed molecular graph facts without reactivity interpretation."""

    input_smiles: str
    canonical_smiles: Optional[str]
    valid: bool
    components: Tuple[MolecularComponentStructure, ...] = ()
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ContextClassification:
    """One definition-derived context facet around a molecular locus."""

    token: str
    attachment_atom_index: int
    fragment_atom_indices: Tuple[int, ...]
    classification_method: str
    facet: str = "fallback"
    semantic_id: str = "other"
    display_token: Optional[str] = None
    subtype: Optional[str] = None
    matched_pattern: Optional[str] = None
    features: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactiveSiteCandidate:
    """Internal molecular-reactivity candidate before overlap resolution."""

    site_type: SiteType
    topology: SiteTopology
    atom_roles: Dict[str, Tuple[int, ...]]
    atom_indices: Tuple[int, ...]
    bond_indices: Tuple[int, ...]
    canonical_signature: str
    render_kind: str
    render_data: Dict[str, Any]
    matched_patterns: Tuple[str, ...]
    details: Dict[str, Any] = field(default_factory=dict)
    context_records: Tuple[ContextClassification, ...] = ()
    availability: str = "available"
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactiveSiteHypothesis:
    """Optional hypothesis describing how a molecular locus may react."""

    hypothesis_id: str
    site_type: SiteType
    topology: SiteTopology
    component_index: int
    atom_indices: Tuple[int, ...]
    bond_indices: Tuple[int, ...]
    canonical_signature: str
    chemist_label: str
    availability: str = "available"
    details: Dict[str, Any] = field(default_factory=dict)
    context_features: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 1.0
    evidence: Tuple[str, ...] = ("molecular_graph_pattern",)
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("reactive-site confidence must be between 0 and 1")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MolecularMotifMatch:
    """Non-exclusive, definition-derived molecular structural motif."""

    motif_id: str
    chemist_label: str
    component_index: int
    atom_indices: Tuple[int, ...]
    anchor_atom_index: int
    tags: Tuple[str, ...] = ()
    matched_pattern: Optional[str] = None
    confidence: float = 1.0
    evidence: Tuple[str, ...] = ("molecular_graph_pattern",)
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactiveSiteEnvironment:
    """Optional local reactivity interpretation for one site hypothesis."""

    hypothesis_id: str
    center_atom_index: int
    reactivity_profile: SiteReactivityProfile
    first_shell: Tuple[str, ...] = ()
    nearby_motifs: Tuple[Dict[str, Any], ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MolecularComponentInterpretation:
    """Optional motifs and reactivity hypotheses for one component."""

    component_index: int
    motifs: Tuple[MolecularMotifMatch, ...] = ()
    reactive_site_hypotheses: Tuple[ReactiveSiteHypothesis, ...] = ()
    reactive_site_environments: Tuple[ReactiveSiteEnvironment, ...] = ()
    connectivity_hypotheses: Tuple[Any, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MolecularInterpretation:
    """Optional molecular annotations excluded from graph identity."""

    components: Tuple[MolecularComponentInterpretation, ...] = ()
    motifs: Tuple[MolecularMotifMatch, ...] = ()
    reactive_site_hypotheses: Tuple[ReactiveSiteHypothesis, ...] = ()
    reactive_site_environments: Tuple[ReactiveSiteEnvironment, ...] = ()
    connectivity_hypotheses: Tuple[Any, ...] = ()
    schema_version: str = "1.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MoleculeAnalysis:
    """Composed molecular structure and optional interpretation."""

    structure: MolecularStructureObservation
    interpretation: MolecularInterpretation = field(
        default_factory=MolecularInterpretation
    )
    schema_version: str = "3.0"

    @property
    def valid(self) -> bool:
        return self.structure.valid

    @property
    def input_smiles(self) -> str:
        return self.structure.input_smiles

    @property
    def canonical_smiles(self) -> Optional[str]:
        return self.structure.canonical_smiles

    @property
    def components(self) -> Tuple[MolecularComponentStructure, ...]:
        return self.structure.components

    @property
    def motifs(self) -> Tuple[MolecularMotifMatch, ...]:
        return self.interpretation.motifs

    @property
    def reactive_site_hypotheses(self) -> Tuple[ReactiveSiteHypothesis, ...]:
        return self.interpretation.reactive_site_hypotheses

    @property
    def reactive_site_environments(self) -> Tuple[ReactiveSiteEnvironment, ...]:
        return self.interpretation.reactive_site_environments

    @property
    def connectivity_hypotheses(self) -> Tuple[Any, ...]:
        return self.interpretation.connectivity_hypotheses

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "structure": self.structure.to_dict(),
            "interpretation": self.interpretation.to_dict(),
        }


__all__ = [
    "ContextClassification",
    "MolecularAtomObservation",
    "MolecularBondObservation",
    "MolecularComponentInterpretation",
    "MolecularComponentStructure",
    "MolecularInterpretation",
    "MolecularMotifMatch",
    "MolecularStructureObservation",
    "MoleculeAnalysis",
    "ReactiveSiteCandidate",
    "ReactiveSiteEnvironment",
    "ReactiveSiteHypothesis",
    "SiteTopology",
    "SiteType",
]
