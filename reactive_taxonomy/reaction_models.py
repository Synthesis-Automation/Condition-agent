"""Typed models for reactive-taxonomy reaction analysis."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

from .models import CompoundAnalysis


REACTION_SIGNATURE_SCHEMA_VERSION = "1.5"


@dataclass(frozen=True)
class ReactionComponent:
    side: Literal["reactant", "agent", "product"]
    component_index: int
    input_smiles: str
    canonical_smiles: str
    atom_mapped: bool
    compound_analysis: CompoundAnalysis


@dataclass(frozen=True)
class ReactionSiteReference:
    side: str
    component_index: int
    site_id: str
    site_type: str
    canonical_signature: str
    chemist_label: str
    availability: str
    atom_roles: Dict[str, Tuple[int, ...]]
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class BondChange:
    change_type: Literal["formed", "broken", "order_changed", "hydrogen_change"]
    atom_1_role: str
    atom_2_role: Optional[str]
    old_order: Optional[str]
    new_order: Optional[str]
    evidence: str


@dataclass(frozen=True)
class ReactionAtomReference:
    """One atom participating in an observed or reconstructed reaction edit."""

    side: str
    component_index: int
    atom_index: int
    atom_map_number: Optional[int]
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    local_environment_id: str
    chiral_tag: Optional[str] = None
    cip_code: Optional[str] = None


@dataclass(frozen=True)
class ReactionEdit:
    """Normalized graph change with atom provenance and evidence."""

    edit_type: Literal["formed", "broken", "order_changed", "hydrogen_change"]
    atom_1: ReactionAtomReference
    atom_2: Optional[ReactionAtomReference]
    old_order: Optional[str]
    new_order: Optional[str]
    evidence: str
    confidence: float


@dataclass(frozen=True)
class ReactionStereoChange:
    """Observed atom or bond stereochemical descriptor correspondence."""

    stereo_type: Literal["atom", "bond"]
    atom_1: ReactionAtomReference
    atom_2: Optional[ReactionAtomReference]
    old_descriptor: Optional[str]
    new_descriptor: Optional[str]
    change_type: Literal[
        "retained", "created", "destroyed", "descriptor_changed"
    ]
    evidence: str
    confidence: float


@dataclass(frozen=True)
class ReactionCompletenessAssessment:
    """Product-atom accounting without requiring reported byproducts."""

    status: Literal["verified", "incomplete", "unresolved"]
    reactant_heavy_atom_count: int
    product_heavy_atom_count: int
    reactant_element_counts: Dict[str, int]
    product_element_counts: Dict[str, int]
    product_element_excess: Dict[str, int]
    reactant_element_excess: Dict[str, int]
    reactant_mapped_heavy_atom_count: int
    product_mapped_heavy_atom_count: int
    shared_mapped_heavy_atom_count: int
    reactant_mapping_coverage: float
    product_mapping_coverage: float
    product_heavy_atom_coverage: Optional[float]
    suspected_missing_reactant: bool
    suspected_insufficient_reactant_multiplicity: bool
    evidence: str
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"


@dataclass(frozen=True)
class ReactionLabelClause:
    """One deterministic, evidence-backed human-readable edit clause."""

    edit_type: Literal["formed", "broken", "order_changed", "hydrogen_change"]
    concise: str
    detailed: str
    elements: Tuple[str, ...]
    atom_map_numbers: Tuple[int, ...]
    old_order: Optional[str]
    new_order: Optional[str]
    evidence: str
    confidence: float


@dataclass(frozen=True)
class ReactionDisplayLabel:
    """Structured display label that never participates in signature identity."""

    concise: str
    detailed: str
    status: Literal[
        "observed_edits",
        "exact_reconstruction",
        "family_overlay",
        "generic_pattern",
        "conflicting_evidence",
        "multi_event",
        "reactant_only",
        "ambiguous_reactants",
        "unavailable",
    ]
    clauses: Tuple[ReactionLabelClause, ...]
    evidence: str
    confidence: float
    warnings: Tuple[str, ...]
    style: str
    definition_version: str
    structural_label: Optional[str] = None
    transformation_label: Optional[str] = None
    grammar_label: Optional[str] = None
    family_label: Optional[str] = None
    pattern_id: Optional[str] = None
    pattern_definition_version: Optional[str] = None
    grammar_id: Optional[str] = None
    contextual_label: Optional[str] = None
    reactant_context_label: Optional[str] = None
    product_context_label: Optional[str] = None
    event_labels: Tuple[str, ...] = ()
    event_count: int = 0
    schema_version: str = "1.3"


@dataclass(frozen=True)
class ReactionSpectatorGroup:
    """Functional group retained outside the selected reaction event."""

    group_id: str
    chemist_label: str
    component_index: int
    atom_indices: Tuple[int, ...]
    nearest_reactive_site_id: Optional[str]
    graph_distance: Optional[int]
    tags: Tuple[str, ...] = ()
    unchanged_evidence: str = "selected_event_exclusion"


@dataclass(frozen=True)
class ReactionPartnerEnvironment:
    """Reaction-family interpretation of one selected substrate role."""

    role: str
    component_index: int
    site_id: str
    handle_token: Optional[str]
    anchor_context: Optional[str]
    chemist_label: str
    steric: Dict[str, Any] = field(default_factory=dict)
    electronic: Dict[str, Any] = field(default_factory=dict)
    nearby_groups: Tuple[Dict[str, Any], ...] = ()
    spectator_group_ids: Tuple[str, ...] = ()
    competing_site_labels: Tuple[str, ...] = ()
    coordination_group_ids: Tuple[str, ...] = ()
    unprotected_xh_group_ids: Tuple[str, ...] = ()
    flags: Tuple[str, ...] = ()
    features: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ReactionFamilyEnvironment:
    """Family overlay assembled from mechanism-neutral observations."""

    family_id: str
    partners: Tuple[ReactionPartnerEnvironment, ...]
    flags: Tuple[str, ...] = ()
    evidence: str = "selected_candidate"
    feature_version: str = "1.0"


@dataclass(frozen=True)
class ReactionPartner:
    """Mechanism-neutral reaction partner with an optional interpreted role."""

    partner_id: str
    component_index: int
    role: Optional[str]
    role_confidence: float
    reactive_site_ids: Tuple[str, ...]
    handle_tokens: Tuple[str, ...]
    anchor_contexts: Tuple[str, ...]
    chemist_label: str
    steric: Dict[str, Any] = field(default_factory=dict)
    electronic: Dict[str, Any] = field(default_factory=dict)
    nearby_groups: Tuple[Dict[str, Any], ...] = ()
    spectator_group_ids: Tuple[str, ...] = ()
    flags: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ProductConnectionEndpoint:
    """One reactant-derived endpoint of a newly formed product bond."""

    role: str
    component_index: int
    atom_index: int
    context: str
    source_site_id: str


@dataclass(frozen=True)
class ProductConnection:
    """Verified product bond with reactant-role and atom provenance."""

    endpoint_1: ProductConnectionEndpoint
    endpoint_2: ProductConnectionEndpoint
    bond_order: str
    connection_type: str
    concise_label: str
    evidence: str
    schema_version: str = "1.0"


@dataclass(frozen=True)
class ProductTransformation:
    """General product transformation supporting one or more graph edits."""

    edits: Tuple[ReactionEdit, ...]
    stereo_changes: Tuple[ReactionStereoChange, ...]
    formed_connection_labels: Tuple[str, ...]
    concise_label: Optional[str]
    exact_product_verified: bool
    evidence: str


@dataclass(frozen=True)
class ReactionTopology:
    """Component and ring topology of an observed or reconstructed event."""

    reaction_scope: Literal[
        "intramolecular", "intermolecular", "mixed", "unimolecular", "unresolved"
    ]
    participating_component_indices: Tuple[int, ...]
    role_component_indices: Dict[str, int]
    same_component_role_groups: Tuple[Tuple[str, ...], ...]
    formed_bond_scopes: Tuple[Literal["intramolecular", "intermolecular"], ...]
    reactant_tether_distances: Tuple[int, ...]
    formed_ring_sizes: Tuple[int, ...]
    ring_count_delta: Optional[int]
    evidence: str
    confidence: float
    schema_version: str = "1.0"


@dataclass(frozen=True)
class ReactionEvent:
    """One connected, atom-provenanced transformation within a reaction."""

    event_id: str
    event_signature_key: str
    edits: Tuple[ReactionEdit, ...]
    partner_ids: Tuple[str, ...]
    reactive_site_ids: Tuple[str, ...]
    formed_bond_types: Tuple[str, ...]
    broken_bond_types: Tuple[str, ...]
    order_changes: Tuple[str, ...]
    hydrogen_changes: Tuple[str, ...]
    formed_connection_labels: Tuple[str, ...]
    topology: ReactionTopology
    transformation_class: Optional[str]
    transformation_confidence: float
    named_family: Optional[str]
    family_confidence: float
    compatible_named_families: Tuple[str, ...]
    evidence: str
    confidence: float
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.1"


@dataclass(frozen=True)
class ReactionEventRelation:
    """Observed structural relationship between two reaction events."""

    event_id_1: str
    event_id_2: str
    relation_type: Literal[
        "independent_sites", "shared_component", "shared_atom", "unresolved"
    ]
    shared_component_indices: Tuple[int, ...]
    shared_atom_map_numbers: Tuple[int, ...]
    evidence: str
    schema_version: str = "1.0"


@dataclass(frozen=True)
class ReactionSignature:
    """Versioned chemistry identity used by conversion and recommendation."""

    signature_id: str
    exact_signature_key: str
    handle_signature_key: str
    transformation_signature_key: str
    bond_edit_signature_key: str
    environment_signature_key: str
    edit_signature: str
    formed_bond_types: Tuple[str, ...]
    broken_bond_types: Tuple[str, ...]
    order_changes: Tuple[str, ...]
    hydrogen_changes: Tuple[str, ...]
    stereo_changes: Tuple[ReactionStereoChange, ...]
    edits: Tuple[ReactionEdit, ...]
    events: Tuple[ReactionEvent, ...]
    event_count: int
    event_scope: Literal["single_event", "multi_event", "unresolved"]
    event_relations: Tuple[ReactionEventRelation, ...]
    partners: Tuple[ReactionPartner, ...]
    product_transformation: Optional[ProductTransformation]
    topology: ReactionTopology
    transformation_class: Optional[str]
    transformation_confidence: float
    named_family: Optional[str]
    family_confidence: float
    compatible_named_families: Tuple[str, ...]
    spectator_groups: Tuple[ReactionSpectatorGroup, ...]
    completeness: ReactionCompletenessAssessment
    global_descriptors: Dict[str, Any]
    warnings: Tuple[str, ...]
    evidence_quality: str
    definition_versions: Dict[str, str]
    schema_version: str = REACTION_SIGNATURE_SCHEMA_VERSION


@dataclass(frozen=True)
class ReactionCandidate:
    grammar_id: str
    transformation_class: str
    role_assignments: Dict[str, ReactionSiteReference]
    predicted_bond_changes: Tuple[BondChange, ...]
    predicted_product_smiles: Optional[str]
    verification: str
    reaction_label: Optional[str]
    compatible_named_families: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class ReactionAnalysis:
    input_reaction_smiles: str
    valid: bool
    reactants: Tuple[ReactionComponent, ...] = ()
    agents: Tuple[ReactionComponent, ...] = ()
    products: Tuple[ReactionComponent, ...] = ()
    candidates: Tuple[ReactionCandidate, ...] = ()
    selected_candidate: Optional[ReactionCandidate] = None
    selected_events: Tuple[ReactionCandidate, ...] = ()
    transformation_class: Optional[str] = None
    compatible_named_families: Tuple[str, ...] = ()
    named_family: Optional[str] = None
    reaction_label: Optional[str] = None
    reaction_label_status: str = "unavailable"
    display_label: Optional[ReactionDisplayLabel] = None
    evidence_quality: str = "unresolved"
    mapped_bond_changes: Tuple[Dict[str, Any], ...] = ()
    spectator_groups: Tuple[ReactionSpectatorGroup, ...] = ()
    family_environment: Optional[ReactionFamilyEnvironment] = None
    product_connection: Optional[ProductConnection] = None
    reaction_topology: Optional[ReactionTopology] = None
    reaction_signature: Optional[ReactionSignature] = None
    reaction_completeness: Optional[ReactionCompletenessAssessment] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "1.8"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


__all__ = [
    "BondChange",
    "ProductConnection",
    "ProductConnectionEndpoint",
    "ProductTransformation",
    "REACTION_SIGNATURE_SCHEMA_VERSION",
    "ReactionAnalysis",
    "ReactionAtomReference",
    "ReactionCandidate",
    "ReactionCompletenessAssessment",
    "ReactionComponent",
    "ReactionDisplayLabel",
    "ReactionEdit",
    "ReactionEvent",
    "ReactionEventRelation",
    "ReactionFamilyEnvironment",
    "ReactionLabelClause",
    "ReactionPartner",
    "ReactionPartnerEnvironment",
    "ReactionSignature",
    "ReactionSiteReference",
    "ReactionSpectatorGroup",
    "ReactionStereoChange",
    "ReactionTopology",
]
