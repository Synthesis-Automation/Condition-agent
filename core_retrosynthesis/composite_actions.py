"""Reusable composite retrosynthetic actions over independently validated operators.

A composite is one logical expansion with two physical reactions. Target-site
coupling is rechecked through route-core atom lineage; catalogue support is a
search prior, never evidence that a new substrate or one-pot procedure works.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from .chemistry import canonical_smiles, digest
from .coupled_strategy_search import V1_ADMITTED_RELATIONSHIPS
from .coupled_route_strategy import extract_coupled_route_strategies
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    GenericTemplateLibrary,
)
from .generic_search import disconnect_generic_target_detailed
from .route_contract import (
    MoleculeOccurrenceNode,
    PlannedRouteAction,
    ReactionRouteTree,
    RouteReactionNode,
    RouteStepEvidence,
    assert_valid_route_tree,
)
from .route_core import build_route_core_projection

COMPOSITE_ACTION_SCHEMA_VERSION = "1.0"
COMPOSITE_ACTION_ALGORITHM_VERSION = "composite_actions.v1"
COMPOSITE_CATALOG_SCHEMA_VERSION = "1.0"
# Existing workbench envelope remains readable during its v1 migration.
COUPLED_STRATEGY_EVALUATION_SCHEMA_VERSION = "1.1"
COUPLED_STRATEGY_EVALUATION_ALGORITHM_VERSION = COMPOSITE_ACTION_ALGORITHM_VERSION
SearchFunction = Callable[
    ..., tuple[tuple[GenericDisconnectionCandidate, ...], GenericSearchDiagnostics]
]


@dataclass(frozen=True)
class CompositeStrategyDefinition:
    """A recurring v1 relationship represented by two generic operators."""

    strategy_id: str
    relationship_class: str
    first_operator_id: str
    second_operator_id: str
    training_patent_ids: tuple[str, ...]
    training_occurrence_count: int
    v2_dependency_counts: tuple[tuple[str, int], ...]

    def __post_init__(self) -> None:
        if self.relationship_class not in V1_ADMITTED_RELATIONSHIPS:
            raise ValueError("only v1 structural relationships can be promoted")
        if len(set(self.training_patent_ids)) < 2:
            raise ValueError("promoted pairs require independent patent support")
        if not all((self.strategy_id, self.first_operator_id, self.second_operator_id)):
            raise ValueError("strategy and operator identities must be nonempty")
        if self.training_occurrence_count < len(set(self.training_patent_ids)):
            raise ValueError("strategy occurrence count must cover supporting patents")
        if any(count < 0 for _, count in self.v2_dependency_counts):
            raise ValueError("dependency counts must be nonnegative")

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["training_patent_ids"] = list(self.training_patent_ids)
        value["v2_dependency_counts"] = dict(self.v2_dependency_counts)
        return value


@dataclass(frozen=True)
class CompositeStrategyCatalog:
    """Versioned reusable strategy definitions, independent of evaluation cases."""

    catalog_id: str
    strategies: tuple[CompositeStrategyDefinition, ...]
    provenance: tuple[tuple[str, str], ...] = ()
    schema_version: str = COMPOSITE_CATALOG_SCHEMA_VERSION
    definition_version: str = "composite_strategy_catalog.v1"

    def __post_init__(self) -> None:
        if self.schema_version != COMPOSITE_CATALOG_SCHEMA_VERSION:
            raise ValueError("unsupported composite catalogue schema")
        if self.definition_version != "composite_strategy_catalog.v1":
            raise ValueError("unsupported composite catalogue definition")
        identities = [item.strategy_id for item in self.strategies]
        if len(set(identities)) != len(identities):
            raise ValueError("duplicate composite strategy identity")
        if self.catalog_id != _catalog_identity(self.strategies, self.provenance):
            raise ValueError("composite catalogue content identity mismatch")

    def to_dict(self) -> dict[str, Any]:
        """Serialize definitions and provenance without evaluation targets."""
        return {
            "artifact_type": "composite_strategy_catalog",
            "catalog_id": self.catalog_id,
            "schema_version": self.schema_version,
            "definition_version": self.definition_version,
            "strategies": [
                item.to_dict()
                for item in sorted(self.strategies, key=lambda x: x.strategy_id)
            ],
            "provenance": dict(self.provenance),
        }


def _catalog_identity(
    strategies: Iterable[CompositeStrategyDefinition],
    provenance: tuple[tuple[str, str], ...],
) -> str:
    payload = []
    for item in sorted(strategies, key=lambda x: x.strategy_id):
        raw = item.to_dict()
        raw["training_patent_ids"] = sorted(set(item.training_patent_ids))
        payload.append(raw)
    return digest(
        "COMPOSITECATALOG1",
        COMPOSITE_CATALOG_SCHEMA_VERSION,
        "composite_strategy_catalog.v1",
        json.dumps(
            {"strategies": payload, "provenance": dict(provenance)},
            sort_keys=True,
            separators=(",", ":"),
        ),
    )


def build_composite_strategy_catalog(
    strategies: Iterable[CompositeStrategyDefinition],
    *,
    provenance: Mapping[str, str] | None = None,
) -> CompositeStrategyCatalog:
    """Freeze reusable operator pairs into a deterministic, case-free catalogue."""
    definitions = tuple(sorted(strategies, key=lambda x: x.strategy_id))
    sources = tuple(sorted((provenance or {}).items()))
    return CompositeStrategyCatalog(
        _catalog_identity(definitions, sources), definitions, sources
    )


def load_composite_strategy_catalog(path: str | Path) -> CompositeStrategyCatalog:
    """Validate a standalone catalogue; evaluation panels require explicit export."""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if value.get("artifact_type") != "composite_strategy_catalog":
        raise ValueError(
            "expected composite strategy catalogue; export the evaluation panel first"
        )
    pairs = tuple(
        CompositeStrategyDefinition(
            strategy_id=raw["strategy_id"],
            relationship_class=raw["relationship_class"],
            first_operator_id=raw["first_operator_id"],
            second_operator_id=raw["second_operator_id"],
            training_patent_ids=tuple(raw["training_patent_ids"]),
            training_occurrence_count=raw["training_occurrence_count"],
            v2_dependency_counts=tuple(sorted(raw["v2_dependency_counts"].items())),
        )
        for raw in value["strategies"]
    )
    return CompositeStrategyCatalog(
        catalog_id=value["catalog_id"],
        strategies=pairs,
        provenance=tuple(sorted(value.get("provenance", {}).items())),
        schema_version=value["schema_version"],
        definition_version=value["definition_version"],
    )


def save_composite_strategy_catalog(
    catalog: CompositeStrategyCatalog, path: str | Path
) -> None:
    """Write a validated catalogue without local reaction or evaluation records."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(catalog.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


@dataclass(frozen=True)
class CompositeDependencyAssessment:
    """Target-specific structural coupling, including unresolved/conflicting evidence."""

    admitted: bool
    status: str
    relationship_class: str
    dependency_class: str
    lineage_status: str
    evidence: Mapping[str, Any]
    warnings: tuple[str, ...] = ()
    physical_step_signatures: tuple[Mapping[str, Any], ...] = ()
    schema_version: str = COMPOSITE_ACTION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return inspectable target-specific dependency evidence."""
        return asdict(self)


@dataclass(frozen=True)
class CompositePhysicalStep:
    """One predicted physical reaction with its own validation and precedent evidence."""

    forward_step_number: int
    reaction_smiles: str
    mapped_reaction_smiles: str
    operator_id: str
    template_id: str
    precedent_reaction_ids: tuple[str, ...]
    forward_validation_status: str
    precursor_compatibility_disposition: str
    reaction_compatibility_disposition: str
    selectivity_warnings: tuple[Mapping[str, Any], ...]
    precursor_compatibility_assessments: tuple[Mapping[str, Any], ...] = ()
    reaction_compatibility_assessments: tuple[Mapping[str, Any], ...] = ()
    evidence_kind: str = "predicted"
    condition_status: str = "not_assessed"


@dataclass(frozen=True)
class _CompositeExpansion:
    """One logical action with both validated physical steps retained."""

    strategy_id: str
    intermediate_smiles: str
    terminal_precursor_smiles: str
    first_operator_id: str
    second_operator_id: str
    first_reaction_smiles: str
    second_reaction_smiles: str
    first_forward_validation_status: str
    second_forward_validation_status: str
    score: float
    dependency: CompositeDependencyAssessment | None = None
    physical_steps: tuple[CompositePhysicalStep, ...] = ()
    route_tree: ReactionRouteTree | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CompositeRetrosyntheticAction:
    """One transferable v1 operator-pair result for an arbitrary target."""

    rank: int
    strategy_id: str
    relationship_class: str
    intermediate_smiles: str
    terminal_precursor_smiles: str
    first_operator_id: str
    second_operator_id: str
    first_reaction_smiles: str
    second_reaction_smiles: str
    first_forward_validation_status: str
    second_forward_validation_status: str
    training_patent_count: int
    training_occurrence_count: int
    v2_dependency_counts: tuple[tuple[str, int], ...]
    score: float
    dependency: CompositeDependencyAssessment | None = None
    physical_steps: tuple[CompositePhysicalStep, ...] = ()
    route_tree: ReactionRouteTree | None = None

    @property
    def physical_step_count(self) -> int:
        """Return the physical route depth represented by this logical action."""
        return 2

    @property
    def physical_step_cost(self) -> int:
        """Charge both physical reactions when comparing with single-step actions."""
        return self.physical_step_count

    @property
    def action_id(self) -> str:
        """Identify normalized chemistry independently of rank and source labels."""
        return digest(
            "COMPOSITEACTION1",
            COMPOSITE_ACTION_SCHEMA_VERSION,
            COMPOSITE_ACTION_ALGORITHM_VERSION,
            self.first_operator_id,
            self.second_operator_id,
            canonical_smiles(self.intermediate_smiles) or "",
            canonical_smiles(self.terminal_precursor_smiles) or "",
            canonical_smiles(self.second_reaction_smiles.split(">")[2]) or "",
        )

    @property
    def retrosynthetic_steps(self) -> tuple[CompositePhysicalStep, ...]:
        """Expand the physical steps in target-to-starting-material order."""
        return tuple(reversed(self.physical_steps))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible query action."""

        value = asdict(self)
        value["v2_dependency_counts"] = dict(self.v2_dependency_counts)
        value.update(
            {
                "action_id": self.action_id,
                "action_kind": "composite",
                "logical_action_count": 1,
                "physical_step_count": self.physical_step_count,
                "physical_step_cost": self.physical_step_cost,
                "condition_compatibility_status": "not_assessed",
                "one_pot_status": "not_assessed",
                "composite_schema_version": COMPOSITE_ACTION_SCHEMA_VERSION,
                "definition_id": "composite_actions.v1",
                "definition_version": "1.0",
                "route_tree": self.route_tree.to_dict() if self.route_tree else None,
                "retrosynthetic_step_numbers": [2, 1],
            }
        )
        return value


@dataclass(frozen=True)
class CompositeSearchDiagnostics:
    """Transparent bounded-search counters for one target query."""

    strategy_count: int
    capable_strategy_count: int
    capability_gap_count: int
    second_step_validation_attempt_count: int
    first_step_validation_attempt_count: int
    fallback_validation_attempt_count: int
    generated_action_count: int
    returned_action_count: int
    returned_fallback_count: int
    dependency_rejected_count: int = 0

    def to_dict(self) -> dict[str, int]:
        """Return JSON-compatible diagnostics."""

        return asdict(self)


@dataclass(frozen=True)
class CompositeSearchResult:
    """Promoted two-step actions with ordinary one-step fallbacks retained."""

    target_smiles: str
    actions: tuple[CompositeRetrosyntheticAction, ...]
    one_step_fallbacks: tuple[dict[str, Any], ...]
    diagnostics: CompositeSearchDiagnostics
    warnings: tuple[str, ...]
    dependency_reviews: tuple[dict[str, Any], ...] = ()
    schema_version: str = COUPLED_STRATEGY_EVALUATION_SCHEMA_VERSION
    algorithm_version: str = COUPLED_STRATEGY_EVALUATION_ALGORITHM_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible target-query result."""

        return {
            "artifact_type": "v1_coupled_strategy_target_query",
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
            "target_smiles": self.target_smiles,
            "valid": bool(self.actions or self.one_step_fallbacks),
            "error": (
                None
                if self.actions or self.one_step_fallbacks
                else "NO_COUPLED_STRATEGY_RESULTS"
            ),
            "actions": [item.to_dict() for item in self.actions],
            "one_step_fallbacks": list(self.one_step_fallbacks),
            "diagnostics": self.diagnostics.to_dict(),
            "warnings": list(self.warnings),
            "dependency_reviews": list(self.dependency_reviews),
        }


def _components(smiles: str) -> tuple[str, ...]:
    canonical = canonical_smiles(smiles)
    return tuple(canonical.split(".")) if canonical else ()


def _merge_terminal_precursors(
    first_precursors: str,
    second_precursors: str,
    intermediate: str,
) -> str | None:
    remaining = list(_components(second_precursors))
    try:
        remaining.remove(intermediate)
    except ValueError:
        return None
    return canonical_smiles(".".join([first_precursors, *remaining]))


def load_composite_action_policy() -> dict[str, Any]:
    """Validate versioned gates and scoring weights; no executable JSON hooks."""
    value = json.loads(
        (Path(__file__).parent / "definitions" / "composite_actions.v1.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        value.get("schema_version"),
        value.get("definition_id"),
        value.get("definition_version"),
    ) != (
        "1.0",
        "composite_actions.v1",
        "1.0",
    ):
        raise ValueError("unsupported composite action policy")
    weights = value["score_weights"]
    if (
        set(weights) != {"first_step", "second_step", "support"}
        or any(
            not isinstance(w, (int, float)) or not math.isfinite(w) or w < 0
            for w in weights.values()
        )
        or not math.isclose(sum(weights.values()), 1.0)
    ):
        raise ValueError("invalid composite scoring weights")
    if (
        set(value["allowed_relationships"]) != V1_ADMITTED_RELATIONSHIPS
        or value["physical_step_cost"] != 2
        or value["dependency_review_limit"] < 1
        or value["support_saturation_patents"] < 1
    ):
        raise ValueError("invalid composite action gates or costs")
    return value


def _physical_step(
    candidate: GenericDisconnectionCandidate, number: int
) -> CompositePhysicalStep:
    return CompositePhysicalStep(
        forward_step_number=number,
        reaction_smiles=candidate.proposed_reaction_smiles,
        mapped_reaction_smiles=candidate.condition_query_reaction_smiles,
        operator_id=candidate.operator_id,
        template_id=candidate.template_id,
        precedent_reaction_ids=candidate.precedent_reaction_ids,
        forward_validation_status=candidate.forward_validation_status,
        precursor_compatibility_disposition=candidate.precursor_compatibility_disposition,
        reaction_compatibility_disposition=candidate.reaction_compatibility_disposition,
        selectivity_warnings=tuple(asdict(w) for w in candidate.selectivity_warnings),
        precursor_compatibility_assessments=tuple(
            asdict(w) for w in candidate.precursor_compatibility_assessments
        ),
        reaction_compatibility_assessments=tuple(
            asdict(w) for w in candidate.reaction_compatibility_assessments
        ),
    )


def _reaction_matches_candidate(candidate: GenericDisconnectionCandidate) -> bool:
    for reaction in (
        candidate.proposed_reaction_smiles,
        candidate.condition_query_reaction_smiles,
    ):
        if not reaction:
            continue
        sides = reaction.split(">")
        if len(sides) != 3:
            return False
        if canonical_smiles(sides[0]) != canonical_smiles(
            candidate.precursor_smiles
        ) or canonical_smiles(sides[2]) != canonical_smiles(candidate.target_smiles):
            return False
    return True


def build_composite_route_tree(
    first: GenericDisconnectionCandidate,
    second: GenericDisconnectionCandidate,
    intermediate_smiles: str,
) -> ReactionRouteTree:
    """Expand a logical pair into an occurrence-preserving predicted route tree.

    Atom maps are retained only from validated candidate reactions. Each additional
    consumer partner stays at its own depth and is not lost in the net projection.
    Leaves remain unresolved for supply; composition never implies stock membership.
    """
    if (
        not _reaction_matches_candidate(first)
        or not _reaction_matches_candidate(second)
        or canonical_smiles(first.target_smiles) != intermediate_smiles
    ):
        raise ValueError("conflicting physical reaction structures")
    if (
        _merge_terminal_precursors(
            first.precursor_smiles, second.precursor_smiles, intermediate_smiles
        )
        is None
    ):
        raise ValueError("intermediate is not a consumer reactant")

    def leaf(smiles: str, node_id: str, depth: int) -> MoleculeOccurrenceNode:
        return MoleculeOccurrenceNode(
            node_id, smiles, depth, False, "not_assessed", "supply_not_assessed"
        )

    def reaction(
        candidate: GenericDisconnectionCandidate,
        number: int,
        depth: int,
        children: tuple[MoleculeOccurrenceNode, ...],
    ) -> RouteReactionNode:
        mapped = candidate.condition_query_reaction_smiles
        return RouteReactionNode(
            reaction_node_id=f"composite:reaction:{number}",
            step_id=f"physical:{number}",
            depth=depth,
            reaction_smiles=mapped or candidate.proposed_reaction_smiles,
            evidence=RouteStepEvidence(
                evidence_kind="predicted",
                source_dataset_id="core_retrosynthesis",
                connectivity_method="validated_operator_composition",
                warnings=("CONDITIONS_AND_SUPPLY_NOT_ASSESSED",),
            ),
            children=children,
            planned_action=PlannedRouteAction(
                candidate.operator_id,
                candidate.disconnection_site_key,
                candidate.template_id,
            ),
        )

    first_children = tuple(
        leaf(s, f"composite:first:{i}", 2)
        for i, s in enumerate(_components(first.precursor_smiles))
    )
    producer = reaction(first, 1, 2, first_children)
    consumer_children = []
    expanded = False
    for i, component in enumerate(_components(second.precursor_smiles)):
        node = leaf(component, f"composite:second:{i}", 1)
        if component == intermediate_smiles and not expanded:
            node = MoleculeOccurrenceNode(
                node.occurrence_id, component, 1, False, "expanded", None, producer
            )
            expanded = True
        consumer_children.append(node)
    consumer = reaction(second, 2, 1, tuple(consumer_children))
    target = canonical_smiles(second.target_smiles)
    if target is None:
        raise ValueError("invalid composite target")
    root = MoleculeOccurrenceNode(
        "composite:target", target, 0, False, "expanded", None, consumer
    )
    tokens = (
        first.operator_id,
        second.operator_id,
        intermediate_smiles,
        canonical_smiles(first.precursor_smiles) or "",
        canonical_smiles(second.precursor_smiles) or "",
    )
    tree = ReactionRouteTree(
        tree_id=digest(
            "COMPOSITETREE1", COMPOSITE_ACTION_SCHEMA_VERSION, target, *tokens
        ),
        route_kind="planned",
        target_smiles=target,
        root=root,
        reaction_count=2,
        maximum_depth=2,
        fingerprint_tokens=tuple(sorted(tokens)),
        connectivity_method="validated_operator_composition",
        warnings=("PREDICTED_COMPOSITE_REQUIRES_CHEMIST_REVIEW",),
    )
    assert_valid_route_tree(tree)
    return tree


def assess_composite_dependency(
    first: GenericDisconnectionCandidate,
    second: GenericDisconnectionCandidate,
    intermediate_smiles: str,
    expected_relationship: str,
) -> tuple[CompositeDependencyAssessment, ReactionRouteTree | None]:
    """Recheck the coupling on the query graph, preserving ambiguous lineage.

    Reuse the same route-core classification used to mine source strategies.
    Complete symmetry alternatives must agree; a catalogue relationship cannot
    override unresolved, independent-site or contradictory graph observations.
    """
    try:
        tree = build_composite_route_tree(first, second, intermediate_smiles)
        projection = build_route_core_projection(tree)
        occurrences = extract_coupled_route_strategies(projection)
    except ValueError as error:
        return CompositeDependencyAssessment(
            False,
            "conflicting",
            "unresolved",
            "unresolved",
            "graph_mismatch",
            {},
            (str(error),),
        ), None
    if len(occurrences) != 1:
        return CompositeDependencyAssessment(
            False,
            "unresolved",
            "unresolved",
            "unresolved",
            "component_unresolved",
            {},
            ("DEPENDENCY_OCCURRENCE_UNRESOLVED",),
        ), tree
    occurrence = occurrences[0]
    supported = occurrence.relationship_class in V1_ADMITTED_RELATIONSHIPS
    agrees = occurrence.relationship_class == expected_relationship
    admitted = supported and agrees and occurrence.evidence.ambiguity_invariant
    status = (
        "verified"
        if admitted
        else "conflicting"
        if supported and not agrees
        else "unresolved"
    )
    if occurrence.relationship_class in {
        "independent_sites",
        "shared_local_environment",
    }:
        status = "independent_sites"
    return CompositeDependencyAssessment(
        admitted,
        status,
        occurrence.relationship_class,
        occurrence.dependency_class,
        occurrence.lineage_status,
        occurrence.evidence.to_dict(),
        occurrence.warnings,
        tuple(
            step.reaction_signature
            for step in sorted(projection.steps, key=lambda s: -s.retrosynthetic_depth)
            if step.reaction_signature is not None
        ),
    ), tree


def search_composite_actions(
    target_smiles: str,
    library: GenericTemplateLibrary,
    strategies: Iterable[CompositeStrategyDefinition],
    *,
    top_k: int = 5,
    max_templates_to_apply: int = 50,
    max_candidates_to_validate: int = 12,
    include_l0: bool = True,
    use_context: bool = True,
    include_one_step_fallbacks: bool = True,
    searcher: SearchFunction = disconnect_generic_target_detailed,
) -> CompositeSearchResult:
    """Apply promoted v1 operator pairs to an arbitrary molecular target.

    Every logical result retains the two independently validated physical
    reactions. The strategy catalog supplies only operator-pair priors; graph
    execution against the query target remains the source of truth.
    """

    if min(top_k, max_templates_to_apply, max_candidates_to_validate) < 1:
        raise ValueError("coupled-strategy search limits must be positive")
    canonical_target = canonical_smiles(target_smiles)
    if canonical_target is None or "." in canonical_target:
        raise ValueError("target must be one valid molecule")
    policy = load_composite_action_policy()
    definitions = tuple(sorted(strategies, key=lambda item: item.strategy_id))
    operator_ids = {item.operator_id for item in library.operators}
    levels = ("L2", "L1", "L0") if include_l0 else ("L2", "L1")

    def query(
        target: str,
        *,
        restricted_operator_ids: Sequence[str] = (),
    ) -> tuple[tuple[GenericDisconnectionCandidate, ...], GenericSearchDiagnostics]:
        return searcher(
            target,
            library,
            operator_ids=restricted_operator_ids,
            levels=levels,
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
        )

    capability_gaps = 0
    capable = 0
    second_attempts = 0
    first_attempts = 0
    first_cache: dict[
        tuple[str, str],
        tuple[tuple[GenericDisconnectionCandidate, ...], GenericSearchDiagnostics],
    ] = {}
    generated: list[tuple[CompositeStrategyDefinition, _CompositeExpansion]] = []
    rejected_dependencies = 0
    reviews: list[dict[str, Any]] = []
    for strategy in definitions:
        if (
            strategy.first_operator_id not in operator_ids
            or strategy.second_operator_id not in operator_ids
        ):
            capability_gaps += 1
            continue
        capable += 1
        second_candidates, second_diagnostics = query(
            canonical_target,
            restricted_operator_ids=(strategy.second_operator_id,),
        )
        second_attempts += second_diagnostics.validation_attempt_count
        for second in second_candidates:
            for intermediate in _components(second.precursor_smiles):
                cache_key = (intermediate, strategy.first_operator_id)
                cached = first_cache.get(cache_key)
                if cached is None:
                    cached = query(
                        intermediate,
                        restricted_operator_ids=(strategy.first_operator_id,),
                    )
                    first_cache[cache_key] = cached
                    first_attempts += cached[1].validation_attempt_count
                for first in cached[0]:
                    terminal = _merge_terminal_precursors(
                        first.precursor_smiles,
                        second.precursor_smiles,
                        intermediate,
                    )
                    if terminal is None:
                        continue
                    if any(
                        step.forward_validation_status
                        not in policy["allowed_forward_validation_statuses"]
                        or step.precursor_compatibility_disposition
                        in policy["blocked_compatibility_dispositions"]
                        or step.reaction_compatibility_disposition
                        in policy["blocked_compatibility_dispositions"]
                        for step in (first, second)
                    ):
                        dependency = CompositeDependencyAssessment(
                            False,
                            "blocked",
                            "unresolved",
                            "unresolved",
                            "not_assessed",
                            {},
                            ("PHYSICAL_STEP_VALIDATION_OR_COMPATIBILITY_FAILED",),
                        )
                        tree = None
                    else:
                        dependency, tree = assess_composite_dependency(
                            first, second, intermediate, strategy.relationship_class
                        )
                    if not dependency.admitted:
                        rejected_dependencies += 1
                        if len(reviews) < policy["dependency_review_limit"]:
                            reviews.append(
                                {
                                    "strategy_id": strategy.strategy_id,
                                    "intermediate_smiles": intermediate,
                                    "first_reaction_smiles": first.proposed_reaction_smiles,
                                    "second_reaction_smiles": second.proposed_reaction_smiles,
                                    "dependency": dependency.to_dict(),
                                }
                            )
                        continue
                    support = min(
                        1.0,
                        math.log1p(len(set(strategy.training_patent_ids)))
                        / math.log1p(policy["support_saturation_patents"]),
                    )
                    generated.append(
                        (
                            strategy,
                            _CompositeExpansion(
                                strategy_id=strategy.strategy_id,
                                intermediate_smiles=intermediate,
                                terminal_precursor_smiles=terminal,
                                first_operator_id=first.operator_id,
                                second_operator_id=second.operator_id,
                                first_reaction_smiles=(first.proposed_reaction_smiles),
                                second_reaction_smiles=(
                                    second.proposed_reaction_smiles
                                ),
                                first_forward_validation_status=(
                                    first.forward_validation_status
                                ),
                                second_forward_validation_status=(
                                    second.forward_validation_status
                                ),
                                score=round(
                                    policy["score_weights"]["first_step"] * first.score
                                    + policy["score_weights"]["second_step"]
                                    * second.score
                                    + policy["score_weights"]["support"] * support,
                                    8,
                                ),
                                dependency=dependency,
                                physical_steps=(
                                    _physical_step(first, 1),
                                    _physical_step(second, 2),
                                ),
                                route_tree=tree,
                            ),
                        )
                    )

    unique: dict[
        tuple[str, str, str, str],
        tuple[CompositeStrategyDefinition, _CompositeExpansion],
    ] = {}
    for strategy, action in generated:
        key = (
            action.first_operator_id,
            action.second_operator_id,
            action.intermediate_smiles,
            action.terminal_precursor_smiles,
        )
        current = unique.get(key)
        if current is None or (
            -action.score,
            strategy.strategy_id,
        ) < (
            -current[1].score,
            current[0].strategy_id,
        ):
            unique[key] = (strategy, action)
    all_ranked = sorted(
        unique.values(),
        key=lambda item: (
            -item[1].score,
            item[1].intermediate_smiles,
            item[1].terminal_precursor_smiles,
            item[0].strategy_id,
        ),
    )
    first_per_strategy = []
    repeated_strategy_actions = []
    selected_strategy_ids = set()
    for item in all_ranked:
        strategy_id = item[0].strategy_id
        if strategy_id in selected_strategy_ids:
            repeated_strategy_actions.append(item)
        else:
            selected_strategy_ids.add(strategy_id)
            first_per_strategy.append(item)
    ranked = (first_per_strategy + repeated_strategy_actions)[:top_k]
    actions = tuple(
        CompositeRetrosyntheticAction(
            rank=rank,
            strategy_id=strategy.strategy_id,
            relationship_class=strategy.relationship_class,
            intermediate_smiles=action.intermediate_smiles,
            terminal_precursor_smiles=action.terminal_precursor_smiles,
            first_operator_id=action.first_operator_id,
            second_operator_id=action.second_operator_id,
            first_reaction_smiles=action.first_reaction_smiles,
            second_reaction_smiles=action.second_reaction_smiles,
            first_forward_validation_status=(action.first_forward_validation_status),
            second_forward_validation_status=(action.second_forward_validation_status),
            training_patent_count=len(strategy.training_patent_ids),
            training_occurrence_count=strategy.training_occurrence_count,
            v2_dependency_counts=strategy.v2_dependency_counts,
            score=action.score,
            dependency=action.dependency,
            physical_steps=action.physical_steps,
            route_tree=action.route_tree,
        )
        for rank, (strategy, action) in enumerate(ranked, 1)
    )

    fallback_attempts = 0
    fallback_values: tuple[dict[str, Any], ...] = ()
    if include_one_step_fallbacks:
        fallbacks, diagnostics = query(canonical_target)
        fallback_attempts = diagnostics.validation_attempt_count
        fallback_values = tuple(
            {"rank": rank, **item.to_dict()} for rank, item in enumerate(fallbacks, 1)
        )
    warnings = (
        "EXPERIMENTAL_PROMOTED_V1_OPERATOR_PAIRS",
        "TWO_PHYSICAL_STEPS_REQUIRE_CHEMIST_REVIEW",
        "STRATEGY_PAIR_DIVERSITY_APPLIED_BEFORE_REALIZATION_VARIANTS",
    )
    return CompositeSearchResult(
        target_smiles=canonical_target,
        actions=actions,
        one_step_fallbacks=fallback_values,
        diagnostics=CompositeSearchDiagnostics(
            strategy_count=len(definitions),
            capable_strategy_count=capable,
            capability_gap_count=capability_gaps,
            second_step_validation_attempt_count=second_attempts,
            first_step_validation_attempt_count=first_attempts,
            fallback_validation_attempt_count=fallback_attempts,
            generated_action_count=len(unique),
            returned_action_count=len(actions),
            returned_fallback_count=len(fallback_values),
            dependency_rejected_count=rejected_dependencies,
        ),
        warnings=warnings,
        dependency_reviews=tuple(reviews),
    )
