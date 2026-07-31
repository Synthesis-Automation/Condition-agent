"""Typed models for reactive-taxonomy reaction analysis."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

from .descriptors.models import SiteReactivityProfile
from .models import CompoundAnalysis


REACTION_SIGNATURE_SCHEMA_VERSION = "3.0"
REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION = "1.3"
REACTION_CORE_PROJECTION_SCHEMA_VERSION = "2.1"
REACTION_CORE_PROJECTION_ALGORITHM_VERSION = "reaction_core_projection.v6"

EditArchetype = Literal[
    "substitution",
    "addition",
    "elimination",
    "bond_order_change",
    "composite",
    "unresolved",
]

BondStateKind = Literal[
    "bond",
    "no_bond",
    "endpoint_absent",
    "unknown",
]

ConnectivityObservationScope = Literal[
    "observed_product",
    "main_product_projection",
    "exact_reconstruction",
    "correspondence_inference",
    "unresolved",
]


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
class PredictedStereoChange:
    """Role-addressed stereochemical outcome emitted by a graph rewrite."""

    stereo_type: Literal["atom", "bond"]
    atom_1_role: str
    atom_2_role: Optional[str]
    old_descriptor: Optional[str]
    new_descriptor: Optional[str]
    change_type: Literal[
        "retained", "inverted", "created", "destroyed", "descriptor_changed"
    ]
    evidence: str


@dataclass(frozen=True)
class RewriteOutcome:
    """One deterministic constitutional outcome of a connectivity rewrite."""

    outcome_id: str
    predicted_product_smiles: Optional[str]
    predicted_bond_changes: Tuple[BondChange, ...]
    predicted_stereo_changes: Tuple[PredictedStereoChange, ...] = ()
    warnings: Tuple[str, ...] = ()


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
class BondState:
    """One typed bond state before or after a connectivity transition."""

    state_kind: BondStateKind
    order: Optional[str]

    def __post_init__(self) -> None:
        if self.state_kind == "bond":
            if not self.order:
                raise ValueError("bond state requires an order")
            object.__setattr__(self, "order", str(self.order).upper())
        elif self.order is not None:
            raise ValueError(f"{self.state_kind} state cannot carry a bond order")


@dataclass(frozen=True)
class BondTransition:
    """Evidence-scoped before/after state for one reactant-origin atom pair."""

    atom_1: ReactionAtomReference
    atom_2: ReactionAtomReference
    before_state: BondState
    after_state: BondState
    delta_units: Optional[int]
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        localized_units = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3}
        if self.delta_units is not None and (
            self.before_state.state_kind not in {"bond", "no_bond"}
            or self.after_state.state_kind not in {"bond", "no_bond"}
        ):
            raise ValueError(
                "delta_units requires definite bond or no-bond endpoint states"
            )
        if self.delta_units is not None:
            before_units = (
                0
                if self.before_state.state_kind == "no_bond"
                else localized_units.get(str(self.before_state.order or ""))
            )
            after_units = (
                0
                if self.after_state.state_kind == "no_bond"
                else localized_units.get(str(self.after_state.order or ""))
            )
            if (
                before_units is None
                or after_units is None
                or after_units - before_units != self.delta_units
            ):
                raise ValueError(
                    "delta_units must match localized before/after bond states"
                )
        if (
            self.after_state.state_kind == "endpoint_absent"
            and self.observation_scope != "main_product_projection"
        ):
            raise ValueError("endpoint_absent requires main_product_projection scope")
        if (
            "unknown"
            in {
                self.before_state.state_kind,
                self.after_state.state_kind,
            }
            and self.observation_scope != "unresolved"
        ):
            raise ValueError("unknown bond state requires unresolved scope")


@dataclass(frozen=True)
class HydrogenDelta:
    """Schema-level hydrogen-count change without invented H correspondence."""

    atom: ReactionAtomReference
    before_count: int
    after_count: int
    delta_count: int
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float

    def __post_init__(self) -> None:
        if self.before_count < 0 or self.after_count < 0:
            raise ValueError("hydrogen counts must be non-negative")
        if self.after_count - self.before_count != self.delta_count:
            raise ValueError("hydrogen delta does not match before/after counts")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")


@dataclass(frozen=True)
class AtomStateTransition:
    """Observed atom-state change kept orthogonal to covalent connectivity."""

    reactant_atom: ReactionAtomReference
    product_atom: ReactionAtomReference
    before_formal_charge: int
    after_formal_charge: int
    before_radical_electrons: Optional[int]
    after_radical_electrons: Optional[int]
    before_isotope: Optional[int]
    after_isotope: Optional[int]
    observation_scope: ConnectivityObservationScope
    evidence: str
    confidence: float

    def __post_init__(self) -> None:
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")


@dataclass(frozen=True)
class ConnectivityEditGraph:
    """Internal canonical connectivity observation used for shadow evaluation."""

    bond_transitions: Tuple[BondTransition, ...]
    hydrogen_deltas: Tuple[HydrogenDelta, ...]
    atom_state_transitions: Tuple[AtomStateTransition, ...]
    edit_component_keys: Tuple[str, ...]
    shadow_key: str
    evidence: str
    confidence: float
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        if not self.shadow_key.startswith("CEG1:"):
            raise ValueError("shadow_key must use the CEG1 namespace")


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
        "retained", "inverted", "created", "destroyed", "descriptor_changed"
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
class ProductFragmentAttachment:
    """One observed bond from a conserved product scaffold to an origin gap."""

    scaffold_atom: ReactionAtomReference
    fragment_atom: ReactionAtomReference
    bond_order: str


@dataclass(frozen=True)
class ProductFragmentSourceCandidate:
    """One structurally compatible supplied source for a product-origin gap."""

    side: Literal["reactant", "agent"]
    component_index: int
    canonical_smiles: str
    evidence: str
    confidence: float


@dataclass(frozen=True)
class ProductAtomProvenance:
    """One product-heavy-atom source assignment or explicit provenance gap."""

    product_atom: ReactionAtomReference
    source_kind: Literal[
        "reactant_correspondence",
        "agent_supported",
        "unresolved",
    ]
    source_atom: Optional[ReactionAtomReference]
    evidence: str
    confidence: float


@dataclass(frozen=True)
class ProductOriginGap:
    """Connected product fragment lacking supplied reactant atom mapping."""

    product_component_index: int
    atom_references: Tuple[ReactionAtomReference, ...]
    internal_bond_types: Tuple[str, ...]
    attachments: Tuple[ProductFragmentAttachment, ...]
    canonical_fragment_smiles: str
    rooted_fragment_smiles: str
    fragment_key: str
    element_counts: Dict[str, int]
    formal_charge: int
    source_status: Literal[
        "unresolved",
        "reactant_supported",
        "agent_supported",
        "ambiguous",
    ]
    source_candidates: Tuple[ProductFragmentSourceCandidate, ...]
    evidence: str
    confidence: float
    schema_version: str = "1.0"


@dataclass(frozen=True)
class PartialProductTransformation:
    """Conservative transformation observed despite incomplete atom provenance."""

    transformation_type: Literal[
        "attachment_replacement",
        "attachment_fragment_replacement",
    ]
    transformation_class: str
    reactant_center: ReactionAtomReference
    product_center: ReactionAtomReference
    removed_attachment: ReactionAtomReference
    added_attachment: ReactionAtomReference
    removed_fragment_atom_indices: Tuple[int, ...]
    removed_fragment_smiles: str
    removed_fragment_key: str
    installed_fragment: ProductOriginGap
    product_atom_provenance: Tuple[ProductAtomProvenance, ...]
    old_order: str
    new_order: str
    conserved_atom_count: int
    product_heavy_atom_coverage: float
    missing_product_atom_elements: Tuple[str, ...]
    evidence: Literal["partial_product_correspondence"]
    confidence: float
    warnings: Tuple[str, ...] = ()
    schema_version: str = "2.0"


@dataclass(frozen=True)
class FragmentSourceRequirement:
    """Product-only fragment that a reported condition component must supply."""

    requirement_id: str
    fragment_key: str
    canonical_fragment_smiles: str
    rooted_fragment_smiles: str
    element_counts: Dict[str, int]
    atom_count: int
    center_element: str
    attachment_element: str
    attachment_bond_order: str
    removed_attachment_element: str
    evidence: Literal["partial_product_correspondence"]
    confidence: float
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
        "partial_product_correspondence",
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
    schema_version: str = "1.4"


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
    nearby_groups: Tuple[Dict[str, Any], ...] = ()
    spectator_group_ids: Tuple[str, ...] = ()
    competing_site_labels: Tuple[str, ...] = ()
    coordination_group_ids: Tuple[str, ...] = ()
    unprotected_xh_group_ids: Tuple[str, ...] = ()
    flags: Tuple[str, ...] = ()
    features: Dict[str, Any] = field(default_factory=dict)
    reactivity_profile: Optional[SiteReactivityProfile] = None


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
    nearby_groups: Tuple[Dict[str, Any], ...] = ()
    spectator_group_ids: Tuple[str, ...] = ()
    flags: Tuple[str, ...] = ()
    reactivity_profile: Optional[SiteReactivityProfile] = None


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
class ReactionEditHypothesis:
    """One unresolved but chemically explicit alternative edit interpretation."""

    hypothesis_id: str
    provider: str
    evidence: str
    confidence: float
    edits: Tuple[ReactionEdit, ...]
    stereo_changes: Tuple[ReactionStereoChange, ...]
    correspondence_count: int
    edit_cost: Tuple[int, int, int]
    topology: Optional[ReactionTopology] = None
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not self.hypothesis_id.startswith("REH1:"):
            raise ValueError("hypothesis_id must use the REH1 namespace")
        if not self.provider or not self.evidence:
            raise ValueError("an edit hypothesis requires provider and evidence")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        if not self.edits:
            raise ValueError("an edit hypothesis requires at least one edit")
        if self.correspondence_count < 1:
            raise ValueError("correspondence_count must be positive")
        if len(self.edit_cost) != 3 or any(value < 0 for value in self.edit_cost):
            raise ValueError("edit_cost must contain three non-negative values")


@dataclass(frozen=True)
class ReactionEvidenceCandidate:
    """Provider-scoped structural evidence before authority reconciliation."""

    provider: str
    status: Literal["verified", "ambiguous", "unresolved", "invalid"]
    evidence: str
    confidence: float
    edits: Tuple[ReactionEdit, ...] = ()
    stereo_changes: Tuple[ReactionStereoChange, ...] = ()
    edit_hypotheses: Tuple[ReactionEditHypothesis, ...] = ()
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        if not self.provider:
            raise ValueError("evidence provider must be present")
        if not self.evidence:
            raise ValueError("evidence description must be present")
        if self.status not in {"verified", "ambiguous", "unresolved", "invalid"}:
            raise ValueError("unsupported evidence-provider status")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        if self.status == "verified" and not self.edits:
            raise ValueError("verified evidence requires edits")
        if self.status == "ambiguous" and not self.edit_hypotheses:
            raise ValueError("ambiguous evidence requires edit hypotheses")


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
    edit_archetype: EditArchetype
    transformation_class: Optional[str]
    transformation_confidence: float
    named_family: Optional[str]
    family_confidence: float
    compatible_named_families: Tuple[str, ...]
    evidence: str
    confidence: float
    warnings: Tuple[str, ...] = ()
    schema_version: str = "1.2"


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
class ReactionFallbackDescriptor:
    """Structure-derived retrieval features that do not assert atom correspondence."""

    descriptor_id: str
    evidence_mode: Literal[
        "verified_signature",
        "partial_product_correspondence",
        "candidate_hypotheses",
        "structure_inventory_only",
    ]
    confidence: float
    retrieval_eligible: bool
    ineligibility_reasons: Tuple[str, ...]
    reactant_component_tokens: Tuple[str, ...]
    product_component_tokens: Tuple[str, ...]
    reactant_site_tokens: Tuple[str, ...]
    product_site_tokens: Tuple[str, ...]
    reactant_group_tokens: Tuple[str, ...]
    product_group_tokens: Tuple[str, ...]
    context_tokens: Tuple[str, ...]
    candidate_grammar_tokens: Tuple[str, ...]
    candidate_transformation_tokens: Tuple[str, ...]
    candidate_handle_tokens: Tuple[str, ...]
    candidate_edit_tokens: Tuple[str, ...]
    candidate_hypothesis_tokens: Tuple[str, ...]
    verified_edit_tokens: Tuple[str, ...]
    reaction_center_core_tokens: Tuple[str, ...]
    reaction_center_radius_1_tokens: Tuple[str, ...]
    reaction_center_radius_2_tokens: Tuple[str, ...]
    reaction_center_radius_3_tokens: Tuple[str, ...]
    bond_inventory_delta_tokens: Tuple[str, ...]
    element_delta_tokens: Tuple[str, ...]
    compatibility_tags: Tuple[str, ...]
    partial_transformation_key: Optional[str]
    source_requirements: Tuple[FragmentSourceRequirement, ...]
    required_condition_source_elements: Tuple[str, ...]
    condition_source_requirement_id: Optional[str]
    definition_versions: Dict[str, str]
    warnings: Tuple[str, ...] = ()
    schema_version: str = REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION


ReactionCoreRemoteClass = Literal[
    "aryl",
    "heteroaryl",
    "alkyl",
    "alkenyl",
    "alkynyl",
    "acyl",
    "ring_aliphatic",
    "heteroatom",
    "generic_R",
]


@dataclass(frozen=True)
class ReactionCoreAtomState:
    """One observed atom state on one side of a minimized reaction center."""

    side: Literal["reactant", "product"]
    component_index: int
    atom_index: int
    atom_map_number: Optional[int]
    element: str
    formal_charge: int
    aromatic: bool
    hybridization: str
    total_hydrogens: int
    heavy_atom_degree: int
    neighbor_tokens: Tuple[str, ...]
    functional_group_ids: Tuple[str, ...]
    concise_label: str
    state_key: str


@dataclass(frozen=True)
class ReactionCoreAtomTransition:
    """Before/after state of one atom participating in the minimized graph."""

    transition_id: str
    atom_map_number: Optional[int]
    before_state: Optional[ReactionCoreAtomState]
    after_state: Optional[ReactionCoreAtomState]
    incident_edit_count: int
    stable_remote_subgraph_count: int
    role: Literal["primary_center", "participant"]


@dataclass(frozen=True)
class ReactionCoreAttachmentPort:
    """One cut connection between an active atom and a remote subgraph."""

    side: Literal["reactant", "product"]
    core_component_index: int
    core_atom_index: int
    core_atom_map_number: Optional[int]
    attachment_atom_index: int
    attachment_atom_map_number: Optional[int]
    attachment_element: str
    bond_order: str


@dataclass(frozen=True)
class ReactionCoreRemoteSubgraph:
    """One connected graph removed from the active minimized reaction graph."""

    subgraph_id: str
    side: Literal["reactant", "product"]
    component_index: int
    atom_indices: Tuple[int, ...]
    atom_map_numbers: Tuple[int, ...]
    remote_class: ReactionCoreRemoteClass
    continuity: Literal[
        "retained", "departing", "appearing", "changed", "unresolved"
    ]
    attachment_ports: Tuple[ReactionCoreAttachmentPort, ...]
    fragment_smiles: str
    fragment_heavy_atom_count: int
    fragment_heteroatom_count: int
    fragment_aromatic_atom_count: int
    functional_group_ids: Tuple[str, ...]


@dataclass(frozen=True)
class ReactionCoreEvent:
    """One connected minimized edit event."""

    event_id: str
    transition_ids: Tuple[str, ...]
    edit_tokens: Tuple[str, ...]


@dataclass(frozen=True)
class ReactionCoreAbstraction:
    """Broad graph-derived motif plus narrower environment limiters."""

    motif_id: str
    motif_key: str
    general_label: str
    limiter_label: str
    motif_tokens: Tuple[str, ...]
    limiter_tokens: Tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.motif_id:
            raise ValueError("reaction-core abstraction requires a motif ID")
        if not self.motif_key.startswith("RCM1:"):
            raise ValueError("motif_key must use the RCM1 namespace")
        if not self.motif_tokens:
            raise ValueError("reaction-core abstraction requires motif tokens")


@dataclass(frozen=True)
class ReactionCoreProjection:
    """Template-free, scaffold-abstracted view of observed reaction edits."""

    core_id: str
    exact_core_key: str
    typed_core_key: str
    shape_core_key: str
    center_transition_key: str
    atom_transitions: Tuple[ReactionCoreAtomTransition, ...]
    events: Tuple[ReactionCoreEvent, ...]
    remote_subgraphs: Tuple[ReactionCoreRemoteSubgraph, ...]
    edit_tokens: Tuple[str, ...]
    participant_tokens: Tuple[str, ...]
    generic_label: str
    abstraction: Optional[ReactionCoreAbstraction]
    active_atom_count: int
    event_count: int
    evidence: str
    evidence_status: Literal["verified", "inferred", "external", "hypothesis"]
    confidence: float
    warnings: Tuple[str, ...] = ()
    algorithm_version: str = REACTION_CORE_PROJECTION_ALGORITHM_VERSION
    schema_version: str = REACTION_CORE_PROJECTION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.core_id.startswith("RCP2:"):
            raise ValueError("core_id must use the RCP2 namespace")
        if not self.exact_core_key.startswith("RCX2:"):
            raise ValueError("exact_core_key must use the RCX2 namespace")
        if not self.typed_core_key.startswith("RCT2:"):
            raise ValueError("typed_core_key must use the RCT2 namespace")
        if not self.shape_core_key.startswith("RSH2:"):
            raise ValueError("shape_core_key must use the RSH2 namespace")
        if not self.center_transition_key.startswith("RCS2:"):
            raise ValueError(
                "center_transition_key must use the RCS2 namespace"
            )
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("confidence must be between 0 and 1")
        if not self.atom_transitions:
            raise ValueError("a reaction-core projection requires active atoms")
        if not any(
            transition.role == "primary_center"
            for transition in self.atom_transitions
        ):
            raise ValueError("a reaction-core projection requires a primary center")


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
    edit_archetype: EditArchetype
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
    rewrite_outcome_id: str
    edit_archetype: EditArchetype
    transformation_class: str
    role_assignments: Dict[str, ReactionSiteReference]
    predicted_bond_changes: Tuple[BondChange, ...]
    predicted_product_smiles: Optional[str]
    verification: str
    reaction_label: Optional[str]
    predicted_stereo_changes: Tuple[PredictedStereoChange, ...] = ()
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
    edit_archetype: EditArchetype = "unresolved"
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
    evidence_candidates: Tuple[ReactionEvidenceCandidate, ...] = ()
    edit_hypotheses: Tuple[ReactionEditHypothesis, ...] = ()
    reaction_signature: Optional[ReactionSignature] = None
    fallback_descriptor: Optional[ReactionFallbackDescriptor] = None
    partial_product_transformation: Optional[PartialProductTransformation] = None
    reaction_completeness: Optional[ReactionCompletenessAssessment] = None
    reaction_core: Optional[ReactionCoreProjection] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "3.5"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


__all__ = [
    "AtomStateTransition",
    "BondState",
    "BondStateKind",
    "BondTransition",
    "BondChange",
    "ConnectivityEditGraph",
    "ConnectivityObservationScope",
    "EditArchetype",
    "FragmentSourceRequirement",
    "HydrogenDelta",
    "RewriteOutcome",
    "PartialProductTransformation",
    "ProductAtomProvenance",
    "ProductFragmentAttachment",
    "ProductFragmentSourceCandidate",
    "ProductOriginGap",
    "ProductConnection",
    "ProductConnectionEndpoint",
    "ProductTransformation",
    "REACTION_CORE_PROJECTION_ALGORITHM_VERSION",
    "REACTION_CORE_PROJECTION_SCHEMA_VERSION",
    "REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION",
    "REACTION_SIGNATURE_SCHEMA_VERSION",
    "ReactionAnalysis",
    "ReactionAtomReference",
    "ReactionCandidate",
    "ReactionCompletenessAssessment",
    "ReactionComponent",
    "ReactionCoreAtomState",
    "ReactionCoreAtomTransition",
    "ReactionCoreAttachmentPort",
    "ReactionCoreAbstraction",
    "ReactionCoreEvent",
    "ReactionCoreProjection",
    "ReactionCoreRemoteClass",
    "ReactionCoreRemoteSubgraph",
    "ReactionDisplayLabel",
    "ReactionEdit",
    "ReactionEditHypothesis",
    "ReactionEvidenceCandidate",
    "ReactionEvent",
    "ReactionEventRelation",
    "ReactionFallbackDescriptor",
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
