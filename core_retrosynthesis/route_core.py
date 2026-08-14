"""Context-preserving reaction-core projections over concrete route trees."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterator, Mapping, Optional, Sequence

from rdkit import Chem
from reactive_taxonomy import (
    build_reaction_display_projection,
    build_reaction_render_context,
    featurize_reaction,
)

from .chemistry import digest
from .route_contract import (
    MoleculeOccurrenceNode,
    ReactionRouteTree,
    RouteReactionNode,
    iter_molecule_occurrences,
)


ROUTE_CORE_SCHEMA_VERSION = "1.0"
ROUTE_CORE_ALGORITHM_VERSION = "route_core_projection.v1"
MAX_LINEAGE_CANDIDATES = 32
ROUTE_CORE_QUALITY_STATUSES = frozenset(
    {"pass", "review", "blocked", "unavailable"}
)
LINEAGE_STATUSES = frozenset(
    {
        "unique",
        "symmetry_ambiguous",
        "bounded_ambiguous",
        "component_unresolved",
        "component_ambiguous",
        "graph_mismatch",
        "unmapped",
    }
)


@dataclass(frozen=True)
class RouteCoreStep:
    """One route reaction enriched with generic structural chemistry."""

    reaction_node_id: str
    step_id: str
    source_reaction_id: Optional[str]
    retrosynthetic_depth: int
    observed_remaining_steps: int
    product_occurrence_id: str
    precursor_occurrence_ids: tuple[str, ...]
    internal_precursor_occurrence_ids: tuple[str, ...]
    terminal_precursor_occurrence_ids: tuple[str, ...]
    reaction_smiles: str
    reagents_smiles: str
    chemistry_valid: bool
    evidence_quality: str
    quality_status: str
    reaction_signature_id: Optional[str]
    reaction_core_id: Optional[str]
    exact_core_key: Optional[str]
    typed_core_key: Optional[str]
    shape_core_key: Optional[str]
    center_transition_key: Optional[str]
    transformation_class: Optional[str]
    named_family: Optional[str]
    event_count: int
    active_atom_count: int
    edit_tokens: tuple[str, ...]
    reaction_signature: Optional[dict[str, Any]]
    reaction_core: Optional[dict[str, Any]]
    minimum_reaction_smiles: Optional[str]
    render_reaction_smiles: Optional[str]
    display_definition_id: Optional[str]
    display_retained_atom_count: int
    display_removed_substituent_count: int
    display_r_group_count: int
    abstracted_reaction_smiles: str
    abstraction_status: str
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.quality_status not in ROUTE_CORE_QUALITY_STATUSES:
            raise ValueError(
                f"Unsupported route-core quality status: {self.quality_status}"
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible step projection."""

        value = asdict(self)
        value["precursor_occurrence_ids"] = list(self.precursor_occurrence_ids)
        value["internal_precursor_occurrence_ids"] = list(
            self.internal_precursor_occurrence_ids
        )
        value["terminal_precursor_occurrence_ids"] = list(
            self.terminal_precursor_occurrence_ids
        )
        value["edit_tokens"] = list(self.edit_tokens)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteCoreStep":
        """Reconstruct one serialized route-core step."""

        fields = dict(value)
        for key in (
            "precursor_occurrence_ids",
            "internal_precursor_occurrence_ids",
            "terminal_precursor_occurrence_ids",
            "edit_tokens",
            "warnings",
        ):
            fields[key] = tuple(str(item) for item in fields.get(key) or ())
        signature = fields.get("reaction_signature")
        core = fields.get("reaction_core")
        fields["reaction_signature"] = (
            dict(signature) if isinstance(signature, Mapping) else None
        )
        fields["reaction_core"] = (
            dict(core) if isinstance(core, Mapping) else None
        )
        return cls(**fields)


@dataclass(frozen=True)
class RouteAtomLineageCandidate:
    """One graph-isomorphic producer-to-consumer atom correspondence."""

    candidate_id: str
    atom_map_pairs: tuple[tuple[int, int], ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible correspondence candidate."""

        return {
            "candidate_id": self.candidate_id,
            "atom_map_pairs": [list(pair) for pair in self.atom_map_pairs],
        }

    @classmethod
    def from_dict(
        cls,
        value: Mapping[str, Any],
    ) -> "RouteAtomLineageCandidate":
        """Reconstruct one atom-lineage candidate."""

        return cls(
            candidate_id=str(value["candidate_id"]),
            atom_map_pairs=tuple(
                (int(pair[0]), int(pair[1]))
                for pair in value.get("atom_map_pairs") or ()
            ),
        )


@dataclass(frozen=True)
class RouteCoreLineageLink:
    """Cross-step identity of one observed intermediate molecule."""

    lineage_id: str
    intermediate_occurrence_id: str
    intermediate_smiles: str
    producer_reaction_node_id: str
    consumer_reaction_node_id: str
    producer_product_component_index: Optional[int]
    consumer_reactant_component_index: Optional[int]
    status: str
    candidate_count: int
    candidate_limit_reached: bool
    candidates: tuple[RouteAtomLineageCandidate, ...]
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in LINEAGE_STATUSES:
            raise ValueError(f"Unsupported route-core lineage status: {self.status}")
        if self.candidate_count < len(self.candidates):
            raise ValueError("Lineage candidate count cannot be below stored count")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible lineage link."""

        return {
            "lineage_id": self.lineage_id,
            "intermediate_occurrence_id": self.intermediate_occurrence_id,
            "intermediate_smiles": self.intermediate_smiles,
            "producer_reaction_node_id": self.producer_reaction_node_id,
            "consumer_reaction_node_id": self.consumer_reaction_node_id,
            "producer_product_component_index": (
                self.producer_product_component_index
            ),
            "consumer_reactant_component_index": (
                self.consumer_reactant_component_index
            ),
            "status": self.status,
            "candidate_count": self.candidate_count,
            "candidate_limit_reached": self.candidate_limit_reached,
            "candidates": [candidate.to_dict() for candidate in self.candidates],
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteCoreLineageLink":
        """Reconstruct one serialized cross-step lineage link."""

        return cls(
            lineage_id=str(value["lineage_id"]),
            intermediate_occurrence_id=str(value["intermediate_occurrence_id"]),
            intermediate_smiles=str(value["intermediate_smiles"]),
            producer_reaction_node_id=str(value["producer_reaction_node_id"]),
            consumer_reaction_node_id=str(value["consumer_reaction_node_id"]),
            producer_product_component_index=(
                int(value["producer_product_component_index"])
                if value.get("producer_product_component_index") is not None
                else None
            ),
            consumer_reactant_component_index=(
                int(value["consumer_reactant_component_index"])
                if value.get("consumer_reactant_component_index") is not None
                else None
            ),
            status=str(value["status"]),
            candidate_count=int(value["candidate_count"]),
            candidate_limit_reached=bool(value["candidate_limit_reached"]),
            candidates=tuple(
                RouteAtomLineageCandidate.from_dict(item)
                for item in value.get("candidates") or ()
            ),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class RouteCoreMotif:
    """One observed two- or three-event route-core subsequence."""

    motif_id: str
    abstraction_level: str
    core_keys: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.abstraction_level not in {"typed", "shape"}:
            raise ValueError("Route-core motif level must be typed or shape")
        if len(self.core_keys) not in {2, 3}:
            raise ValueError("Route-core motifs must contain two or three events")
        motif_prefix, core_prefix = (
            ("RCTMOTIF1:", "RCT2:")
            if self.abstraction_level == "typed"
            else ("RSHMOTIF1:", "RSH2:")
        )
        if not self.motif_id.startswith(motif_prefix):
            raise ValueError("Route-core motif ID contradicts its level")
        if any(not key.startswith(core_prefix) for key in self.core_keys):
            raise ValueError("Route-core motif contains an incompatible core key")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route motif."""

        return {
            "motif_id": self.motif_id,
            "abstraction_level": self.abstraction_level,
            "length": len(self.core_keys),
            "core_keys": list(self.core_keys),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteCoreMotif":
        """Reconstruct one serialized route motif."""

        motif = cls(
            motif_id=str(value["motif_id"]),
            abstraction_level=str(value["abstraction_level"]),
            core_keys=tuple(str(item) for item in value.get("core_keys") or ()),
        )
        if int(value.get("length") or len(motif.core_keys)) != len(
            motif.core_keys
        ):
            raise ValueError("Route-core motif length is inconsistent")
        return motif


@dataclass(frozen=True)
class RouteCoreProjection:
    """A non-executable, minimized chemistry projection of one route."""

    route_core_id: str
    source_tree_id: str
    source_route_id: Optional[str]
    patent_id: Optional[str]
    split: Optional[str]
    target_smiles: str
    route_kind: str
    reaction_count: int
    maximum_depth: int
    exact_route_core_key: str
    typed_route_core_key: str
    shape_route_core_key: str
    steps: tuple[RouteCoreStep, ...]
    lineage_links: tuple[RouteCoreLineageLink, ...]
    motifs: tuple[RouteCoreMotif, ...]
    typed_motif_keys: tuple[str, ...]
    shape_motif_keys: tuple[str, ...]
    fully_chemistry_resolved: bool
    fully_lineage_connected: bool
    warnings: tuple[str, ...]
    algorithm_version: str = ROUTE_CORE_ALGORITHM_VERSION
    schema_version: str = ROUTE_CORE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route-core projection."""

        return {
            "route_core_id": self.route_core_id,
            "source_tree_id": self.source_tree_id,
            "source_route_id": self.source_route_id,
            "patent_id": self.patent_id,
            "split": self.split,
            "target_smiles": self.target_smiles,
            "route_kind": self.route_kind,
            "reaction_count": self.reaction_count,
            "maximum_depth": self.maximum_depth,
            "exact_route_core_key": self.exact_route_core_key,
            "typed_route_core_key": self.typed_route_core_key,
            "shape_route_core_key": self.shape_route_core_key,
            "steps": [step.to_dict() for step in self.steps],
            "lineage_links": [link.to_dict() for link in self.lineage_links],
            "motifs": [motif.to_dict() for motif in self.motifs],
            "typed_motif_keys": list(self.typed_motif_keys),
            "shape_motif_keys": list(self.shape_motif_keys),
            "fully_chemistry_resolved": self.fully_chemistry_resolved,
            "fully_lineage_connected": self.fully_lineage_connected,
            "warnings": list(self.warnings),
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteCoreProjection":
        """Reconstruct and validate one serialized route-core projection."""

        projection = cls(
            route_core_id=str(value["route_core_id"]),
            source_tree_id=str(value["source_tree_id"]),
            source_route_id=value.get("source_route_id"),
            patent_id=value.get("patent_id"),
            split=value.get("split"),
            target_smiles=str(value["target_smiles"]),
            route_kind=str(value["route_kind"]),
            reaction_count=int(value["reaction_count"]),
            maximum_depth=int(value["maximum_depth"]),
            exact_route_core_key=str(value["exact_route_core_key"]),
            typed_route_core_key=str(value["typed_route_core_key"]),
            shape_route_core_key=str(value["shape_route_core_key"]),
            steps=tuple(
                RouteCoreStep.from_dict(item) for item in value.get("steps") or ()
            ),
            lineage_links=tuple(
                RouteCoreLineageLink.from_dict(item)
                for item in value.get("lineage_links") or ()
            ),
            motifs=tuple(
                RouteCoreMotif.from_dict(item)
                for item in value.get("motifs") or ()
            ),
            typed_motif_keys=tuple(
                str(item) for item in value.get("typed_motif_keys") or ()
            ),
            shape_motif_keys=tuple(
                str(item) for item in value.get("shape_motif_keys") or ()
            ),
            fully_chemistry_resolved=bool(value["fully_chemistry_resolved"]),
            fully_lineage_connected=bool(value["fully_lineage_connected"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            algorithm_version=str(
                value.get("algorithm_version") or ROUTE_CORE_ALGORITHM_VERSION
            ),
            schema_version=str(
                value.get("schema_version") or ROUTE_CORE_SCHEMA_VERSION
            ),
        )
        assert_valid_route_core_projection(projection)
        return projection


@dataclass(frozen=True)
class RouteCoreValidationReport:
    """Structural validation report for one route-core projection."""

    valid: bool
    issues: tuple[str, ...]


def validate_route_core_projection(
    projection: RouteCoreProjection,
) -> RouteCoreValidationReport:
    """Validate identity, counts, references, ordering and schema versions."""

    issues: list[str] = []
    if projection.schema_version != ROUTE_CORE_SCHEMA_VERSION:
        issues.append("unsupported_schema_version")
    if projection.algorithm_version != ROUTE_CORE_ALGORITHM_VERSION:
        issues.append("unsupported_algorithm_version")
    if not projection.route_core_id.startswith("ROUTECORE1:"):
        issues.append("invalid_route_core_id")
    for key, prefix, issue in (
        (projection.exact_route_core_key, "RCXROUTE1:", "invalid_exact_key"),
        (projection.typed_route_core_key, "RCTROUTE1:", "invalid_typed_key"),
        (projection.shape_route_core_key, "RSHROUTE1:", "invalid_shape_key"),
    ):
        if not key.startswith(prefix):
            issues.append(issue)
    if len(projection.steps) != projection.reaction_count:
        issues.append("reaction_count_mismatch")
    reaction_ids = tuple(step.reaction_node_id for step in projection.steps)
    if len(reaction_ids) != len(set(reaction_ids)):
        issues.append("duplicate_reaction_node_id")
    if tuple(sorted(projection.steps, key=_step_order_key)) != projection.steps:
        issues.append("step_order_mismatch")
    known_ids = set(reaction_ids)
    lineage_ids: set[str] = set()
    for link in projection.lineage_links:
        if link.lineage_id in lineage_ids:
            issues.append("duplicate_lineage_id")
        lineage_ids.add(link.lineage_id)
        if (
            link.producer_reaction_node_id not in known_ids
            or link.consumer_reaction_node_id not in known_ids
        ):
            issues.append("unknown_lineage_reaction")
    if tuple(sorted(set(projection.typed_motif_keys))) != (
        projection.typed_motif_keys
    ):
        issues.append("typed_motif_order_mismatch")
    if tuple(sorted(set(projection.shape_motif_keys))) != (
        projection.shape_motif_keys
    ):
        issues.append("shape_motif_order_mismatch")
    motif_ids = tuple(motif.motif_id for motif in projection.motifs)
    if len(motif_ids) != len(set(motif_ids)):
        issues.append("duplicate_motif_id")
    if tuple(sorted(motif_ids)) != motif_ids:
        issues.append("motif_order_mismatch")
    if {
        motif.motif_id
        for motif in projection.motifs
        if motif.abstraction_level == "typed"
    } != set(projection.typed_motif_keys):
        issues.append("typed_motif_definition_mismatch")
    if {
        motif.motif_id
        for motif in projection.motifs
        if motif.abstraction_level == "shape"
    } != set(projection.shape_motif_keys):
        issues.append("shape_motif_definition_mismatch")
    return RouteCoreValidationReport(valid=not issues, issues=tuple(issues))


def assert_valid_route_core_projection(projection: RouteCoreProjection) -> None:
    """Raise when a route-core projection violates its public contract."""

    report = validate_route_core_projection(projection)
    if not report.valid:
        raise ValueError("; ".join(report.issues))


def _step_order_key(step: RouteCoreStep) -> tuple[int, str]:
    return (-step.retrosynthetic_depth, step.reaction_node_id)


def _remaining_steps(node: MoleculeOccurrenceNode) -> int:
    if node.reaction is None:
        return 0
    return 1 + max((_remaining_steps(child) for child in node.reaction.children), default=0)


def _display_projection(analysis: Any) -> tuple[Any | None, tuple[str, ...]]:
    if analysis.observation is None or analysis.reaction_core is None:
        return None, ()
    try:
        context = build_reaction_render_context(
            observation=analysis.observation,
            reactants=analysis.reactants,
            agents=analysis.agents,
            products=analysis.products,
            signature=analysis.reaction_signature,
            partial_transformation=analysis.partial_product_transformation,
            interpretation=analysis.interpretation,
        )
        return build_reaction_display_projection(context), ()
    except (RuntimeError, TypeError, ValueError) as exc:
        return None, (f"ROUTE_CORE_DISPLAY_FAILED:{type(exc).__name__}",)


def _quality_status(analysis: Any) -> str:
    if not analysis.valid:
        return "blocked"
    core = analysis.reaction_core
    if core is None:
        return "unavailable"
    return str(core.quality.status)


def _build_step(
    molecule: MoleculeOccurrenceNode,
    reaction: RouteReactionNode,
    analysis: Any,
) -> RouteCoreStep:
    signature = analysis.reaction_signature
    core = analysis.reaction_core
    display, display_warnings = _display_projection(analysis)
    components = tuple(display.reactants + display.products) if display else ()
    warnings = tuple(
        sorted(set(tuple(analysis.warnings) + display_warnings))
    )
    return RouteCoreStep(
        reaction_node_id=reaction.reaction_node_id,
        step_id=reaction.step_id,
        source_reaction_id=reaction.evidence.source_reaction_id,
        retrosynthetic_depth=reaction.depth,
        observed_remaining_steps=_remaining_steps(molecule),
        product_occurrence_id=molecule.occurrence_id,
        precursor_occurrence_ids=tuple(child.occurrence_id for child in reaction.children),
        internal_precursor_occurrence_ids=tuple(
            child.occurrence_id for child in reaction.children if child.reaction is not None
        ),
        terminal_precursor_occurrence_ids=tuple(
            child.occurrence_id for child in reaction.children if child.reaction is None
        ),
        reaction_smiles=reaction.reaction_smiles,
        reagents_smiles=reaction.evidence.reagents_smiles,
        chemistry_valid=bool(analysis.valid),
        evidence_quality=str(analysis.evidence_quality),
        quality_status=_quality_status(analysis),
        reaction_signature_id=(signature.signature_id if signature else None),
        reaction_core_id=(core.core_id if core else None),
        exact_core_key=(core.exact_core_key if core else None),
        typed_core_key=(core.typed_core_key if core else None),
        shape_core_key=(core.shape_core_key if core else None),
        center_transition_key=(core.center_transition_key if core else None),
        transformation_class=(
            signature.transformation_class if signature else analysis.transformation_class
        ),
        named_family=(signature.named_family if signature else analysis.named_family),
        event_count=(int(core.event_count) if core else 0),
        active_atom_count=(int(core.active_atom_count) if core else 0),
        edit_tokens=tuple(core.edit_tokens) if core else (),
        reaction_signature=(asdict(signature) if signature else None),
        reaction_core=(asdict(core) if core else None),
        minimum_reaction_smiles=(display.minimum_reaction_smiles if display else None),
        render_reaction_smiles=(display.render_reaction_smiles if display else None),
        display_definition_id=(display.definition_id if display else None),
        display_retained_atom_count=sum(
            len(component.retained_atom_indices) for component in components
        ),
        display_removed_substituent_count=sum(
            component.removed_substituent_count for component in components
        ),
        display_r_group_count=sum(component.r_group_count for component in components),
        abstracted_reaction_smiles=reaction.evidence.abstracted_reaction_smiles,
        abstraction_status=reaction.evidence.abstraction_status,
        warnings=warnings,
    )


def _canonical_unmapped(smiles: str) -> Optional[str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    value = Chem.Mol(molecule)
    for atom in value.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(value, canonical=True, isomericSmiles=True)


def _matching_components(
    components: Sequence[Any],
    smiles: str,
) -> tuple[Any, ...]:
    return tuple(
        component
        for component in components
        if _canonical_unmapped(component.input_smiles) == smiles
    )


def _lineage_candidates(
    producer_smiles: str,
    consumer_smiles: str,
) -> tuple[str, int, bool, tuple[RouteAtomLineageCandidate, ...], tuple[str, ...]]:
    producer = Chem.MolFromSmiles(producer_smiles)
    consumer = Chem.MolFromSmiles(consumer_smiles)
    if producer is None or consumer is None:
        return "graph_mismatch", 0, False, (), ("LINEAGE_PARSE_FAILED",)
    if producer.GetNumAtoms() != consumer.GetNumAtoms():
        return "graph_mismatch", 0, False, (), ("LINEAGE_ATOM_COUNT_MISMATCH",)
    query = Chem.Mol(producer)
    target = Chem.Mol(consumer)
    for molecule in (query, target):
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
    raw_matches = target.GetSubstructMatches(
        query,
        uniquify=False,
        useChirality=True,
        maxMatches=MAX_LINEAGE_CANDIDATES + 1,
    )
    if not raw_matches:
        return "graph_mismatch", 0, False, (), ("LINEAGE_ISOMORPHISM_MISSING",)
    mappings: set[tuple[tuple[int, int], ...]] = set()
    unmapped = False
    for match in raw_matches:
        pairs = []
        for producer_index, consumer_index in enumerate(match):
            producer_map = producer.GetAtomWithIdx(producer_index).GetAtomMapNum()
            consumer_map = consumer.GetAtomWithIdx(consumer_index).GetAtomMapNum()
            if producer_map <= 0 or consumer_map <= 0:
                unmapped = True
                continue
            pairs.append((int(producer_map), int(consumer_map)))
        mappings.add(tuple(sorted(pairs)))
    if unmapped or not mappings:
        return "unmapped", len(mappings), False, (), ("LINEAGE_ATOM_MAP_MISSING",)
    ordered = tuple(sorted(mappings))
    limited = len(raw_matches) > MAX_LINEAGE_CANDIDATES
    stored = ordered[:MAX_LINEAGE_CANDIDATES]
    candidates = tuple(
        RouteAtomLineageCandidate(
            candidate_id=digest(
                "ROUTELINEAGECAND1",
                *(f"{left}>{right}" for left, right in mapping),
            ),
            atom_map_pairs=mapping,
        )
        for mapping in stored
    )
    status = (
        "bounded_ambiguous"
        if limited
        else "unique"
        if len(ordered) == 1
        else "symmetry_ambiguous"
    )
    warnings = (
        ("LINEAGE_CANDIDATE_LIMIT_REACHED",)
        if limited
        else ("SYMMETRY_EQUIVALENT_LINEAGE",)
        if len(ordered) > 1
        else ()
    )
    return status, len(ordered), limited, candidates, warnings


def _build_lineage_link(
    intermediate: MoleculeOccurrenceNode,
    producer: RouteReactionNode,
    consumer: RouteReactionNode,
    analyses: Mapping[str, Any],
) -> RouteCoreLineageLink:
    producer_analysis = analyses[producer.reaction_node_id]
    consumer_analysis = analyses[consumer.reaction_node_id]
    producer_components = _matching_components(
        producer_analysis.products,
        intermediate.smiles,
    )
    consumer_components = _matching_components(
        consumer_analysis.reactants,
        intermediate.smiles,
    )
    lineage_id = digest(
        "ROUTELINEAGE1",
        intermediate.occurrence_id,
        producer.reaction_node_id,
        consumer.reaction_node_id,
    )
    if not producer_components or not consumer_components:
        return RouteCoreLineageLink(
            lineage_id=lineage_id,
            intermediate_occurrence_id=intermediate.occurrence_id,
            intermediate_smiles=intermediate.smiles,
            producer_reaction_node_id=producer.reaction_node_id,
            consumer_reaction_node_id=consumer.reaction_node_id,
            producer_product_component_index=None,
            consumer_reactant_component_index=None,
            status="component_unresolved",
            candidate_count=0,
            candidate_limit_reached=False,
            candidates=(),
            warnings=("LINEAGE_COMPONENT_NOT_FOUND",),
        )
    if len(producer_components) != 1 or len(consumer_components) != 1:
        return RouteCoreLineageLink(
            lineage_id=lineage_id,
            intermediate_occurrence_id=intermediate.occurrence_id,
            intermediate_smiles=intermediate.smiles,
            producer_reaction_node_id=producer.reaction_node_id,
            consumer_reaction_node_id=consumer.reaction_node_id,
            producer_product_component_index=None,
            consumer_reactant_component_index=None,
            status="component_ambiguous",
            candidate_count=0,
            candidate_limit_reached=False,
            candidates=(),
            warnings=("LINEAGE_COMPONENT_AMBIGUOUS",),
        )
    producer_component = producer_components[0]
    consumer_component = consumer_components[0]
    status, count, limited, candidates, warnings = _lineage_candidates(
        producer_component.input_smiles,
        consumer_component.input_smiles,
    )
    return RouteCoreLineageLink(
        lineage_id=lineage_id,
        intermediate_occurrence_id=intermediate.occurrence_id,
        intermediate_smiles=intermediate.smiles,
        producer_reaction_node_id=producer.reaction_node_id,
        consumer_reaction_node_id=consumer.reaction_node_id,
        producer_product_component_index=producer_component.component_index,
        consumer_reactant_component_index=consumer_component.component_index,
        status=status,
        candidate_count=count,
        candidate_limit_reached=limited,
        candidates=candidates,
        warnings=warnings,
    )


def _route_encoding(
    node: MoleculeOccurrenceNode,
    step_by_reaction: Mapping[str, RouteCoreStep],
    field: str,
    *,
    include_leaf_identity: bool,
) -> str:
    if node.reaction is None:
        return f"leaf:{node.smiles}" if include_leaf_identity else "leaf"
    step = step_by_reaction[node.reaction.reaction_node_id]
    key = str(getattr(step, field) or f"missing:{step.quality_status}")
    children = tuple(
        sorted(
            _route_encoding(
                child,
                step_by_reaction,
                field,
                include_leaf_identity=include_leaf_identity,
            )
            for child in node.reaction.children
        )
    )
    return f"reaction:{key}({','.join(children)})"


def _path_core_keys(
    tree: ReactionRouteTree,
    step_by_reaction: Mapping[str, RouteCoreStep],
    field: str,
) -> tuple[tuple[str, ...], ...]:
    paths: list[tuple[str, ...]] = []

    def visit(node: MoleculeOccurrenceNode, keys: tuple[str, ...]) -> None:
        if node.reaction is None:
            if keys:
                paths.append(tuple(reversed(keys)))
            return
        step = step_by_reaction[node.reaction.reaction_node_id]
        key = str(getattr(step, field) or f"missing:{step.quality_status}")
        for child in node.reaction.children:
            visit(child, keys + (key,))

    visit(tree.root, ())
    return tuple(sorted(set(paths)))


def _motifs(
    paths: Sequence[Sequence[str]],
    namespace: str,
    abstraction_level: str,
) -> tuple[RouteCoreMotif, ...]:
    motifs: dict[str, RouteCoreMotif] = {}
    for path in paths:
        for length in (2, 3):
            for start in range(max(0, len(path) - length + 1)):
                motif = tuple(path[start : start + length])
                if len(motif) == length and all(
                    not key.startswith("missing:") for key in motif
                ):
                    motif_id = digest(namespace, *motif)
                    motifs[motif_id] = RouteCoreMotif(
                        motif_id=motif_id,
                        abstraction_level=abstraction_level,
                        core_keys=motif,
                    )
    return tuple(motifs[key] for key in sorted(motifs))


def build_route_core_projection(tree: ReactionRouteTree) -> RouteCoreProjection:
    """Build a chemistry-rich, non-executable projection of one route tree."""

    analyses: dict[str, Any] = {}
    steps: list[RouteCoreStep] = []
    molecule_by_reaction: dict[str, MoleculeOccurrenceNode] = {}
    for molecule in iter_molecule_occurrences(tree):
        reaction = molecule.reaction
        if reaction is None:
            continue
        analysis = featurize_reaction(reaction.reaction_smiles)
        analyses[reaction.reaction_node_id] = analysis
        molecule_by_reaction[reaction.reaction_node_id] = molecule
        steps.append(_build_step(molecule, reaction, analysis))
    ordered_steps = tuple(sorted(steps, key=_step_order_key))
    step_by_reaction = {step.reaction_node_id: step for step in ordered_steps}

    lineage_links: list[RouteCoreLineageLink] = []
    for molecule in iter_molecule_occurrences(tree):
        consumer = molecule.reaction
        if consumer is None:
            continue
        for child in consumer.children:
            if child.reaction is None:
                continue
            lineage_links.append(
                _build_lineage_link(
                    child,
                    child.reaction,
                    consumer,
                    analyses,
                )
            )
    ordered_links = tuple(sorted(lineage_links, key=lambda link: link.lineage_id))

    exact_encoding = _route_encoding(
        tree.root,
        step_by_reaction,
        "exact_core_key",
        include_leaf_identity=True,
    )
    typed_encoding = _route_encoding(
        tree.root,
        step_by_reaction,
        "typed_core_key",
        include_leaf_identity=False,
    )
    shape_encoding = _route_encoding(
        tree.root,
        step_by_reaction,
        "shape_core_key",
        include_leaf_identity=False,
    )
    exact_key = digest("RCXROUTE1", exact_encoding)
    typed_key = digest("RCTROUTE1", typed_encoding)
    shape_key = digest("RSHROUTE1", shape_encoding)
    typed_paths = _path_core_keys(tree, step_by_reaction, "typed_core_key")
    shape_paths = _path_core_keys(tree, step_by_reaction, "shape_core_key")
    typed_motifs = _motifs(typed_paths, "RCTMOTIF1", "typed")
    shape_motifs = _motifs(shape_paths, "RSHMOTIF1", "shape")
    fully_chemistry_resolved = all(
        step.quality_status in {"pass", "review"}
        and step.reaction_signature_id is not None
        and step.reaction_core_id is not None
        for step in ordered_steps
    )
    fully_lineage_connected = all(
        link.status in {"unique", "symmetry_ambiguous"}
        for link in ordered_links
    )
    warnings: list[str] = [
        "Route-core projections are display, evaluation, and retrieval evidence; "
        "they are not executable multistep templates."
    ]
    if not fully_chemistry_resolved:
        warnings.append("ROUTE_CORE_CHEMISTRY_PARTIAL")
    if not fully_lineage_connected:
        warnings.append("ROUTE_CORE_LINEAGE_PARTIAL")
    if any(link.status == "symmetry_ambiguous" for link in ordered_links):
        warnings.append("ROUTE_CORE_LINEAGE_HAS_SYMMETRY_AMBIGUITY")
    projection = RouteCoreProjection(
        route_core_id=digest("ROUTECORE1", exact_key),
        source_tree_id=tree.tree_id,
        source_route_id=tree.source_route_id,
        patent_id=tree.patent_id,
        split=tree.split,
        target_smiles=tree.target_smiles,
        route_kind=tree.route_kind,
        reaction_count=tree.reaction_count,
        maximum_depth=tree.maximum_depth,
        exact_route_core_key=exact_key,
        typed_route_core_key=typed_key,
        shape_route_core_key=shape_key,
        steps=ordered_steps,
        lineage_links=ordered_links,
        motifs=tuple(sorted(typed_motifs + shape_motifs, key=lambda item: item.motif_id)),
        typed_motif_keys=tuple(motif.motif_id for motif in typed_motifs),
        shape_motif_keys=tuple(motif.motif_id for motif in shape_motifs),
        fully_chemistry_resolved=fully_chemistry_resolved,
        fully_lineage_connected=fully_lineage_connected,
        warnings=tuple(warnings),
    )
    assert_valid_route_core_projection(projection)
    return projection


def iter_route_core_steps(
    projection: RouteCoreProjection,
) -> Iterator[RouteCoreStep]:
    """Yield route steps in forward-synthesis order."""

    yield from projection.steps


__all__ = [
    "LINEAGE_STATUSES",
    "MAX_LINEAGE_CANDIDATES",
    "ROUTE_CORE_ALGORITHM_VERSION",
    "ROUTE_CORE_QUALITY_STATUSES",
    "ROUTE_CORE_SCHEMA_VERSION",
    "RouteAtomLineageCandidate",
    "RouteCoreLineageLink",
    "RouteCoreMotif",
    "RouteCoreProjection",
    "RouteCoreStep",
    "RouteCoreValidationReport",
    "assert_valid_route_core_projection",
    "build_route_core_projection",
    "iter_route_core_steps",
    "validate_route_core_projection",
]
