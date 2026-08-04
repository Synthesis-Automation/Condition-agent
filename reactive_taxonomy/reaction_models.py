"""Typed models for reactive-taxonomy reaction analysis."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

from .descriptors.models import SiteReactivityProfile
from .models import MoleculeAnalysis, MolecularStructureObservation


REACTION_SIGNATURE_SCHEMA_VERSION = "3.4"
REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION = "3.0"
REACTION_CORE_PROJECTION_SCHEMA_VERSION = "3.1"
REACTION_SUBSTITUENT_PROFILE_SCHEMA_VERSION = "1.0"
REACTION_CORE_EVENT_RELATION_SCHEMA_VERSION = "1.0"
REACTION_PATTERN_MATCH_SCHEMA_VERSION = "4.1"
REACTION_R_GROUP_FUNCTIONAL_CONTEXT_SCHEMA_VERSION = "1.0"
REACTION_CORE_PROJECTION_ALGORITHM_VERSION = "reaction_core_projection.v11"
REACTION_RING_CHANGE_SCHEMA_VERSION = "1.0"
REACTION_TOPOLOGY_SCHEMA_VERSION = "2.0"

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
    "correspondence_inference",
    "unresolved",
]


@dataclass(frozen=True)
class ReactionComponent:
    """Parsed reaction component with optional molecular annotations.

    This composed component is used by the public analysis and interpretation
    layers.  The structural observation receives ``ReactionStructureComponent``
    instead, so site hypotheses and motif matches cannot leak into reaction
    evidence or identity.
    """

    side: Literal["reactant", "agent", "product"]
    component_index: int
    input_smiles: str
    canonical_smiles: str
    atom_mapped: bool
    molecule_analysis: MoleculeAnalysis
    inferred_copy_of_component_index: Optional[int] = None

    def structure_only(self) -> "ReactionStructureComponent":
        """Project this component onto graph facts only."""
        return ReactionStructureComponent(
            side=self.side,
            component_index=self.component_index,
            input_smiles=self.input_smiles,
            canonical_smiles=self.canonical_smiles,
            atom_mapped=self.atom_mapped,
            molecular_structure=self.molecule_analysis.structure,
            inferred_copy_of_component_index=(
                self.inferred_copy_of_component_index
            ),
        )


@dataclass(frozen=True)
class ReactionStructureComponent:
    """One reaction component containing only parsed molecular graph facts."""

    side: Literal["reactant", "agent", "product"]
    component_index: int
    input_smiles: str
    canonical_smiles: str
    atom_mapped: bool
    molecular_structure: MolecularStructureObservation
    inferred_copy_of_component_index: Optional[int] = None


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
    schema_version: str = "1.1"


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
class RenderedReactionLabel:
    """The sole chemist-facing reaction label plus non-display provenance."""

    text: str
    status: Literal[
        "structural_equation",
        "conflicting_evidence",
        "multi_event",
        "ring_formation",
        "partial_product_correspondence",
        "unavailable",
    ]
    basis: Literal[
        "reaction_sites",
        "local_context",
        "literal_edits",
        "ring_topology",
        "partial_product_correspondence",
        "unavailable",
    ]
    evidence: str
    confidence: float
    warnings: Tuple[str, ...]
    style: str
    definition_version: str
    event_count: int = 0
    core_event_ids: Tuple[str, ...] = ()
    substituent_profile_ids: Tuple[str, ...] = ()
    pattern_ids: Tuple[str, ...] = ()
    unclassified_edit_indices: Tuple[int, ...] = ()
    schema_version: str = "1.1"

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
class ReactionRGroupPortDistance:
    """Graph distance from one unchanged motif to one R-group attachment."""

    attachment_atom_index: int
    substituent_profile_id: str
    bond_distance: int

    def __post_init__(self) -> None:
        if self.attachment_atom_index < 0 or self.bond_distance < 0:
            raise ValueError("R-group port indices and distances must be non-negative")
        if not self.substituent_profile_id.startswith("RSP1:"):
            raise ValueError("R-group port distance requires an RSP1 profile ID")


@dataclass(frozen=True)
class ReactionRGroupFunctionalGroup:
    """One unchanged functional group contained by a remote core subgraph."""

    motif_id: str
    chemist_label: str
    atom_indices: Tuple[int, ...]
    tags: Tuple[str, ...]
    nearest_reactive_site_id: Optional[str]
    distance_to_reactive_site: Optional[int]
    port_distances: Tuple[ReactionRGroupPortDistance, ...]
    unchanged_evidence: str

    def __post_init__(self) -> None:
        if not self.motif_id or not self.atom_indices:
            raise ValueError("R-group functional group requires motif provenance")
        if tuple(sorted(set(self.atom_indices))) != self.atom_indices:
            raise ValueError("R-group functional-group atoms must be unique and sorted")
        if tuple(sorted(set(self.tags))) != self.tags:
            raise ValueError("R-group functional-group tags must be unique and sorted")
        if (
            self.distance_to_reactive_site is not None
            and self.distance_to_reactive_site < 0
        ):
            raise ValueError("reactive-site distance must be non-negative")
        if not self.port_distances:
            raise ValueError("R-group functional group requires a port distance")


@dataclass(frozen=True)
class ReactionRGroupFunctionalContext:
    """Optional motif overlay for one graph-provenanced remote subgraph."""

    context_id: str
    remote_subgraph_id: str
    side: Literal["reactant"]
    component_index: int
    remote_class: ReactionCoreRemoteClass
    continuity: Literal[
        "retained", "departing", "appearing", "changed", "unresolved"
    ]
    remote_atom_indices: Tuple[int, ...]
    attachment_profile_ids: Tuple[str, ...]
    functional_groups: Tuple[ReactionRGroupFunctionalGroup, ...]
    evidence: str = "spectator_motif_contained_in_remote_subgraph"
    schema_version: str = REACTION_R_GROUP_FUNCTIONAL_CONTEXT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.context_id.startswith("RGFC1:"):
            raise ValueError("R-group functional-context IDs require RGFC1")
        if not self.remote_subgraph_id.startswith("RCR2:"):
            raise ValueError("R-group context requires a remote-subgraph ID")
        if self.component_index < 0 or not self.remote_atom_indices:
            raise ValueError("R-group context requires remote atom provenance")
        if tuple(sorted(set(self.remote_atom_indices))) != self.remote_atom_indices:
            raise ValueError("R-group context atoms must be unique and sorted")
        if tuple(sorted(set(self.attachment_profile_ids))) != (
            self.attachment_profile_ids
        ):
            raise ValueError("R-group profile IDs must be unique and sorted")
        if not self.attachment_profile_ids or not self.functional_groups:
            raise ValueError("R-group context requires ports and functional groups")


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
    """Graph-derived participating component with optional semantic role."""

    partner_id: str
    component_index: int
    role: Optional[str]
    role_confidence: float
    anchor_contexts: Tuple[str, ...]
    chemist_label: str


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
    evidence: str
    schema_version: str = "1.0"


@dataclass(frozen=True)
class ProductTransformation:
    """General product transformation supporting one or more graph edits."""

    edits: Tuple[ReactionEdit, ...]
    stereo_changes: Tuple[ReactionStereoChange, ...]
    formed_connection_labels: Tuple[str, ...]
    exact_product_verified: bool
    evidence: str


@dataclass(frozen=True)
class ReactionRingChange:
    """Minimal graph facts for one newly formed cycle."""

    change_id: str
    change_type: Literal["formed"]
    atom_references: Tuple[ReactionAtomReference, ...]
    element_sequence: Tuple[str, ...]
    bond_orders_after: Tuple[str, ...]
    source_component_indices: Tuple[int, ...]
    formed_bond_types: Tuple[str, ...]
    aromatic_after: bool
    evidence: str
    confidence: float
    schema_version: str = REACTION_RING_CHANGE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.change_id.startswith("RRG1:"):
            raise ValueError("ring-change ID must use the RRG1 namespace")
        if len(self.atom_references) < 3:
            raise ValueError("a formed ring requires at least three atoms")
        if len(self.element_sequence) != len(self.atom_references):
            raise ValueError("ring element sequence must match its atoms")
        if len(self.bond_orders_after) != len(self.atom_references):
            raise ValueError("ring bond-order sequence must close the cycle")
        if not self.source_component_indices:
            raise ValueError("ring change requires reactant component provenance")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("ring-change confidence must be in [0, 1]")

    @property
    def ring_size(self) -> int:
        return len(self.atom_references)


@dataclass(frozen=True)
class ReactionTopology:
    """Component and ring topology of an observed or reconstructed event."""

    reaction_scope: Literal[
        "intramolecular", "intermolecular", "mixed", "unimolecular", "unresolved"
    ]
    participating_component_indices: Tuple[int, ...]
    formed_bond_scopes: Tuple[Literal["intramolecular", "intermolecular"], ...]
    reactant_tether_distances: Tuple[int, ...]
    formed_ring_sizes: Tuple[int, ...]
    ring_count_delta: Optional[int]
    evidence: str
    confidence: float
    ring_changes: Tuple[ReactionRingChange, ...] = ()
    schema_version: str = REACTION_TOPOLOGY_SCHEMA_VERSION


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
    atom_correspondence: Tuple[Tuple[int, int, int, int], ...] = ()
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
    active_atom_ids: Tuple[str, ...]
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
    schema_version: str = "2.0"


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
class ReactionCoreAromaticSubstituentRelation:
    """Position of one remote aromatic substituent relative to an active atom."""

    reactive_atom_index: int
    ring_atom_index: int
    aromatic_distance: int
    positional_relation: Literal["ipso", "ortho", "meta", "para", "other"]
    substituent_attachment_atom_index: int
    substituent_element: str
    substituent_bond_order: str
    substituent_fragment_smiles: str

    def __post_init__(self) -> None:
        if self.aromatic_distance < 0:
            raise ValueError("aromatic substituent distance cannot be negative")


@dataclass(frozen=True)
class ReactionCoreSubstituentProfile:
    """Port-specific chemistry of one graph fragment omitted as an R group."""

    profile_id: str
    base_class: ReactionCoreRemoteClass
    attachment_element: str
    attachment_bond_order: str
    attachment_aromatic: bool
    attachment_hybridization: str
    carbon_substitution: Literal[
        "not_carbon",
        "not_applicable",
        "methyl",
        "primary",
        "secondary",
        "tertiary",
        "quaternary",
        "unresolved",
    ]
    cyclic: bool
    ring_sizes: Tuple[int, ...]
    benzylic: bool
    allylic: bool
    propargylic: bool
    alpha_branch_count: int
    beta_branch_count: int
    radius_1_heteroatoms: Tuple[str, ...]
    radius_2_heteroatoms: Tuple[str, ...]
    aromatic_substituent_relations: Tuple[
        ReactionCoreAromaticSubstituentRelation, ...
    ]
    local_environment_key: str
    feature_tokens: Tuple[str, ...]
    definition_version: str
    algorithm_version: str
    schema_version: str = REACTION_SUBSTITUENT_PROFILE_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.profile_id.startswith("RSP1:"):
            raise ValueError("substituent profile IDs must use the RSP1 namespace")
        if not self.local_environment_key.startswith("RSE1:"):
            raise ValueError(
                "substituent environment keys must use the RSE1 namespace"
            )
        if self.alpha_branch_count < 0 or self.beta_branch_count < 0:
            raise ValueError("substituent branching counts cannot be negative")
        if tuple(sorted(set(self.ring_sizes))) != self.ring_sizes:
            raise ValueError("substituent ring sizes must be unique and sorted")
        if tuple(sorted(set(self.feature_tokens))) != self.feature_tokens:
            raise ValueError("substituent feature tokens must be unique and sorted")


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
    radical_electrons: int
    isotope: int
    neighbor_tokens: Tuple[str, ...]
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
class ReactionCoreStateChange:
    """One explicit atom or bond property change within the observed core."""

    change_id: str
    change_type: Literal[
        "hydrogen",
        "formal_charge",
        "radical",
        "isotope",
        "aromaticity",
        "hybridization",
        "atom_stereochemistry",
        "bond_stereochemistry",
    ]
    atom_map_numbers: Tuple[int, ...]
    elements: Tuple[str, ...]
    before_value: str
    after_value: str
    evidence: str


@dataclass(frozen=True)
class ReactionCoreQuality:
    """Deterministic trust assessment for one constructed reaction core."""

    status: Literal["pass", "review", "blocked"]
    active_atom_mapping_coverage: float
    checked_edit_fraction: float
    edit_count: int
    heavy_atom_edit_count: int
    event_count: int
    passed_checks: Tuple[str, ...]
    review_reasons: Tuple[str, ...]
    blocking_reasons: Tuple[str, ...]
    definition_version: str

    def __post_init__(self) -> None:
        if not 0.0 <= self.active_atom_mapping_coverage <= 1.0:
            raise ValueError("active atom mapping coverage must be in [0, 1]")
        if not 0.0 <= self.checked_edit_fraction <= 1.0:
            raise ValueError("checked edit fraction must be in [0, 1]")
        if self.status == "pass" and (
            self.review_reasons or self.blocking_reasons
        ):
            raise ValueError("passing core quality cannot contain cautions")
        if self.status == "blocked" and not self.blocking_reasons:
            raise ValueError("blocked core quality requires a blocking reason")


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
    substituent_profile: ReactionCoreSubstituentProfile


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


@dataclass(frozen=True)
class ReactionCoreEvent:
    """One connected minimized edit event."""

    event_id: str
    transition_ids: Tuple[str, ...]
    edit_tokens: Tuple[str, ...]
    edit_indices: Tuple[int, ...]
    reactant_component_indices: Tuple[int, ...]
    product_component_indices: Tuple[int, ...]


@dataclass(frozen=True)
class ReactionCoreEventPath:
    """Shortest observed molecular path between two reaction-core events."""

    side: Literal["reactant", "product"]
    component_index: int
    start_atom_index: int
    end_atom_index: int
    atom_indices: Tuple[int, ...]
    bond_count: int

    def __post_init__(self) -> None:
        if self.bond_count < 0:
            raise ValueError("reaction event path bond count cannot be negative")
        if self.atom_indices and len(self.atom_indices) != self.bond_count + 1:
            raise ValueError("reaction event path size does not match bond count")


@dataclass(frozen=True)
class ReactionCoreEventRelation:
    """Graph-derived relationship between two minimized reaction events."""

    event_id_1: str
    event_id_2: str
    relation_type: Literal[
        "same_component", "shared_active_atom", "independent", "unresolved"
    ]
    shared_reactant_component_indices: Tuple[int, ...]
    shared_product_component_indices: Tuple[int, ...]
    shortest_paths: Tuple[ReactionCoreEventPath, ...]
    evidence: str
    schema_version: str = REACTION_CORE_EVENT_RELATION_SCHEMA_VERSION


@dataclass(frozen=True)
class ReactionCoreProjection:
    """Template-free, scaffold-abstracted view of observed reaction edits."""

    core_id: str
    exact_core_key: str
    typed_core_key: str
    shape_core_key: str
    center_transition_key: str
    mapping_equivalence_key: str
    atom_transitions: Tuple[ReactionCoreAtomTransition, ...]
    state_changes: Tuple[ReactionCoreStateChange, ...]
    events: Tuple[ReactionCoreEvent, ...]
    event_relations: Tuple[ReactionCoreEventRelation, ...]
    remote_subgraphs: Tuple[ReactionCoreRemoteSubgraph, ...]
    edit_tokens: Tuple[str, ...]
    participant_tokens: Tuple[str, ...]
    quality: ReactionCoreQuality
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
        if not self.mapping_equivalence_key.startswith("RME1:"):
            raise ValueError(
                "mapping_equivalence_key must use the RME1 namespace"
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
        event_ids = {event.event_id for event in self.events}
        if any(
            relation.event_id_1 not in event_ids
            or relation.event_id_2 not in event_ids
            for relation in self.event_relations
        ):
            raise ValueError("reaction-core event relation refers to an unknown event")


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
class ReactionPatternMatch:
    """One optional interpretation supported by an existing observation."""

    pattern_id: str
    tier: Literal["generic", "synthesis"]
    confidence: float
    specificity: int
    display_importance: int
    matched_edit_indices: Tuple[int, ...]
    matched_core_event_ids: Tuple[str, ...]
    matched_substituent_profile_ids: Tuple[str, ...]
    covered_core_event_fraction: float
    evidence: Tuple[str, ...]
    occurrence_count: int = 1
    compatible_named_families: Tuple[str, ...] = ()
    requires_condition_evidence: bool = False
    warnings: Tuple[str, ...] = ()
    schema_version: str = REACTION_PATTERN_MATCH_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not self.pattern_id:
            raise ValueError("reaction pattern requires an ID")
        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("reaction pattern confidence must be between 0 and 1")
        if self.specificity < 0 or self.display_importance < 0:
            raise ValueError("reaction pattern ranking values cannot be negative")
        if self.occurrence_count < 1:
            raise ValueError("reaction pattern occurrence count must be positive")
        if tuple(sorted(set(self.matched_edit_indices))) != self.matched_edit_indices:
            raise ValueError("reaction pattern edit indices must be unique and sorted")
        if tuple(sorted(set(self.matched_core_event_ids))) != (
            self.matched_core_event_ids
        ):
            raise ValueError("reaction pattern core-event IDs must be unique and sorted")
        if tuple(sorted(set(self.matched_substituent_profile_ids))) != (
            self.matched_substituent_profile_ids
        ):
            raise ValueError(
                "reaction pattern substituent-profile IDs must be unique and sorted"
            )
        if not 0.0 <= self.covered_core_event_fraction <= 1.0:
            raise ValueError("reaction pattern core-event coverage must be in [0, 1]")


@dataclass(frozen=True)
class ReactionObservation:
    """Interpretation-independent facts derived from one reaction graph."""

    input_reaction_smiles: str
    valid: bool
    reactants: Tuple[ReactionStructureComponent, ...] = ()
    agents: Tuple[ReactionStructureComponent, ...] = ()
    products: Tuple[ReactionStructureComponent, ...] = ()
    edits: Tuple[ReactionEdit, ...] = ()
    stereo_changes: Tuple[ReactionStereoChange, ...] = ()
    evidence_quality: str = "unresolved"
    evidence_confidence: float = 0.0
    connectivity_edit_graph: Optional[ConnectivityEditGraph] = None
    evidence_candidates: Tuple[ReactionEvidenceCandidate, ...] = ()
    edit_hypotheses: Tuple[ReactionEditHypothesis, ...] = ()
    mapped_bond_changes: Tuple[Dict[str, Any], ...] = ()
    topology: Optional[ReactionTopology] = None
    completeness: Optional[ReactionCompletenessAssessment] = None
    core: Optional[ReactionCoreProjection] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "3.0"

    def __post_init__(self) -> None:
        if not 0.0 <= self.evidence_confidence <= 1.0:
            raise ValueError("observation confidence must be between 0 and 1")
        if not self.valid and self.edits:
            raise ValueError("an invalid reaction observation cannot contain edits")


@dataclass(frozen=True)
class ReactionInterpretation:
    """Optional semantic and family interpretation of an observation."""

    pattern_matches: Tuple[ReactionPatternMatch, ...] = ()
    primary_pattern_id: Optional[str] = None
    co_occurring_pattern_ids: Tuple[str, ...] = ()
    partners: Tuple[ReactionPartner, ...] = ()
    compatible_named_families: Tuple[str, ...] = ()
    named_family: Optional[str] = None
    family_environment: Optional[ReactionFamilyEnvironment] = None
    product_connection: Optional[ProductConnection] = None
    spectator_groups: Tuple[ReactionSpectatorGroup, ...] = ()
    r_group_functional_contexts: Tuple[
        ReactionRGroupFunctionalContext, ...
    ] = ()
    evidence_quality: str = "unresolved"
    warnings: Tuple[str, ...] = ()
    schema_version: str = "7.1"

    def __post_init__(self) -> None:
        if self.named_family and self.named_family not in (
            self.compatible_named_families
        ):
            raise ValueError("named family must be one compatible family")
        pattern_ids = {pattern.pattern_id for pattern in self.pattern_matches}
        if self.primary_pattern_id and self.primary_pattern_id not in pattern_ids:
            raise ValueError("primary pattern must refer to a matched pattern")
        if any(value not in pattern_ids for value in self.co_occurring_pattern_ids):
            raise ValueError("co-occurring patterns must refer to matched patterns")


@dataclass(frozen=True)
class ReactionAnalysis:
    input_reaction_smiles: str
    valid: bool
    reactants: Tuple[ReactionComponent, ...] = ()
    agents: Tuple[ReactionComponent, ...] = ()
    products: Tuple[ReactionComponent, ...] = ()
    edit_archetype: EditArchetype = "unresolved"
    transformation_class: Optional[str] = None
    compatible_named_families: Tuple[str, ...] = ()
    named_family: Optional[str] = None
    reaction_label: Optional[RenderedReactionLabel] = None
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
    observation: Optional[ReactionObservation] = None
    interpretation: Optional[ReactionInterpretation] = None
    warnings: Tuple[str, ...] = ()
    error: Optional[str] = None
    schema_version: str = "10.0"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


__all__ = [
    "AtomStateTransition",
    "BondState",
    "BondStateKind",
    "BondTransition",
    "ConnectivityEditGraph",
    "ConnectivityObservationScope",
    "EditArchetype",
    "FragmentSourceRequirement",
    "HydrogenDelta",
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
    "REACTION_RING_CHANGE_SCHEMA_VERSION",
    "REACTION_R_GROUP_FUNCTIONAL_CONTEXT_SCHEMA_VERSION",
    "REACTION_SIGNATURE_SCHEMA_VERSION",
    "REACTION_TOPOLOGY_SCHEMA_VERSION",
    "ReactionAnalysis",
    "ReactionAtomReference",
    "ReactionCompletenessAssessment",
    "ReactionComponent",
    "ReactionStructureComponent",
    "ReactionCoreAtomState",
    "ReactionCoreAtomTransition",
    "ReactionCoreAttachmentPort",
    "ReactionCoreEvent",
    "ReactionCoreProjection",
    "ReactionCoreRemoteClass",
    "ReactionCoreRemoteSubgraph",
    "RenderedReactionLabel",
    "ReactionEdit",
    "ReactionEditHypothesis",
    "ReactionEvidenceCandidate",
    "ReactionEvent",
    "ReactionEventRelation",
    "ReactionFallbackDescriptor",
    "ReactionFamilyEnvironment",
    "ReactionInterpretation",
    "ReactionObservation",
    "ReactionPartner",
    "ReactionPatternMatch",
    "ReactionPartnerEnvironment",
    "ReactionRingChange",
    "ReactionRGroupFunctionalContext",
    "ReactionRGroupFunctionalGroup",
    "ReactionRGroupPortDistance",
    "ReactionSignature",
    "ReactionSiteReference",
    "ReactionSpectatorGroup",
    "ReactionStereoChange",
    "ReactionTopology",
]
