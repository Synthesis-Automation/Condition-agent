"""Replay observed route actions through validated single-step retrosynthesis."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import re
from typing import Any, Callable, Mapping, Optional, Sequence

from .chemistry import digest
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    OperatorLadderDiagnostics,
)
from .generic_search import disconnect_operator_ladder_detailed
from .hierarchical_ranking import CompletionPriorIndex, build_completion_prior_index
from .observed_route_action import (
    ObservedRouteActionLabel,
    build_observed_route_action_label,
    normalize_observed_reaction,
)
from .route_contract import ReactionRouteTree, iter_molecule_occurrences


ROUTE_ACTION_EVALUATION_SCHEMA_VERSION = "2.0"
ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION = "route_action_replay.v2"
_PATENT_TOKEN = re.compile(r"US0*(\d+)([A-Z]\d)", re.IGNORECASE)


@dataclass(frozen=True)
class RouteActionEvaluationConfig:
    """Deterministic search budget used for every observed route state."""

    top_k: int = 25
    max_templates_to_apply: int = 500
    max_candidates_to_validate: int = 100
    use_context: bool = True
    include_l0: bool = True
    diversify: bool = True
    use_hierarchical_ranking: bool = True
    minimum_candidates_per_level: int = 0
    lazy_validation: bool = False
    run_search: bool = True

    def __post_init__(self) -> None:
        if self.top_k < 1:
            raise ValueError("top-k must be positive")
        if self.max_templates_to_apply < 1:
            raise ValueError("maximum templates to apply must be positive")
        if self.max_candidates_to_validate < 1:
            raise ValueError("maximum candidates to validate must be positive")
        if self.max_candidates_to_validate < self.top_k:
            raise ValueError("candidate validation budget must cover top-k")
        if self.minimum_candidates_per_level < 0:
            raise ValueError("minimum candidates per level cannot be negative")

    @property
    def config_id(self) -> str:
        """Return a content-derived identity for the replay search policy."""

        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return digest("RAEC1", payload)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible configuration."""

        return {**asdict(self), "config_id": self.config_id}


@dataclass(frozen=True)
class RouteActionCandidate:
    """Compact, training-ready view of one validated candidate action."""

    candidate_rank: int
    site_rank: int
    operator_rank: int
    synthon_rank: int
    strategy_rank: int
    precursor_smiles: str
    abstraction_level: str
    template_id: str
    operator_id: str
    disconnection_site_key: str
    synthon_signature: str
    strategy_id: str
    score: float
    structural_score_band: int
    strategic_complexity_score: float
    strategic_class: str
    independent_reference_support: int
    precursor_compatibility_disposition: str
    exact_precursor_match: bool
    strategy_match: bool
    site_match: bool
    operator_match: bool
    synthon_match: bool
    supervision_label: str
    source_patent_precedent_overlap: bool
    precedent_reaction_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible candidate action."""

        value = asdict(self)
        value["precedent_reaction_ids"] = list(self.precedent_reaction_ids)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteActionCandidate":
        """Reconstruct a candidate action from JSON data."""

        fields = dict(value)
        fields["precedent_reaction_ids"] = tuple(
            str(item) for item in fields.get("precedent_reaction_ids") or ()
        )
        return cls(**fields)


@dataclass(frozen=True)
class RouteActionStepEvaluation:
    """Observed positive action and replay results for one route state."""

    evaluation_id: str
    reaction_node_id: str
    step_id: str
    source_reaction_id: Optional[str]
    product_occurrence_id: str
    retrosynthetic_depth: int
    observed_remaining_steps: int
    reaction_smiles: str
    route_target_smiles: str
    observed_action: ObservedRouteActionLabel
    search_status: str
    candidate_count: int
    exact_precursor_rank: Optional[int]
    site_rank: Optional[int]
    operator_rank: Optional[int]
    synthon_rank: Optional[int]
    strategy_rank: Optional[int]
    outcome: str
    source_patent_precedent_overlap: bool
    search_diagnostics: dict[str, Any]
    candidates: tuple[RouteActionCandidate, ...]
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.candidate_count != len(self.candidates):
            raise ValueError("candidate count contradicts retained candidates")
        if self.search_status not in {
            "searched",
            "not_search_eligible",
            "not_run",
        }:
            raise ValueError("unsupported route-action search status")
        if self.search_status == "searched" and not self.observed_action.search_eligible:
            raise ValueError("searched action lacks a learnable identity")
        if self.search_status != "searched" and self.candidates:
            raise ValueError("unsearched route actions cannot contain candidates")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible step evaluation."""

        value = asdict(self)
        value["observed_action"] = self.observed_action.to_dict()
        value["candidates"] = [candidate.to_dict() for candidate in self.candidates]
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteActionStepEvaluation":
        """Reconstruct a step evaluation from JSON data."""

        fields = dict(value)
        observed_action = fields.get("observed_action")
        if not isinstance(observed_action, Mapping):
            raise ValueError("step evaluation requires an observed-action label")
        fields["observed_action"] = ObservedRouteActionLabel.from_dict(
            observed_action
        )
        fields["candidates"] = tuple(
            RouteActionCandidate.from_dict(item)
            for item in fields.get("candidates") or ()
        )
        fields["warnings"] = tuple(
            str(item) for item in fields.get("warnings") or ()
        )
        fields["search_diagnostics"] = dict(
            fields.get("search_diagnostics") or {}
        )
        return cls(**fields)


@dataclass(frozen=True)
class RouteActionEvaluation:
    """Replay benchmark results for one occurrence-preserving route tree."""

    evaluation_id: str
    tree_id: str
    source_route_id: Optional[str]
    patent_id: Optional[str]
    split: Optional[str]
    target_smiles: str
    reaction_count: int
    maximum_depth: int
    search_config_id: str
    steps: tuple[RouteActionStepEvaluation, ...]
    schema_version: str = ROUTE_ACTION_EVALUATION_SCHEMA_VERSION
    algorithm_version: str = ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != ROUTE_ACTION_EVALUATION_SCHEMA_VERSION:
            raise ValueError("unsupported route-action evaluation schema")
        if self.algorithm_version != ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION:
            raise ValueError("unsupported route-action evaluation algorithm")
        if self.reaction_count != len(self.steps):
            raise ValueError("route-action step count contradicts route")
        if len({step.evaluation_id for step in self.steps}) != len(self.steps):
            raise ValueError("route-action evaluation IDs must be unique")

    @property
    def strategy_verified_step_count(self) -> int:
        """Return the number of steps with verified STRAT1 supervision."""

        return sum(step.observed_action.strategy_verified for step in self.steps)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route evaluation."""

        value = asdict(self)
        value["steps"] = [step.to_dict() for step in self.steps]
        return value

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteActionEvaluation":
        """Reconstruct and validate one serialized route evaluation."""

        fields = dict(value)
        fields["steps"] = tuple(
            RouteActionStepEvaluation.from_dict(item)
            for item in fields.get("steps") or ()
        )
        return cls(**fields)


RouteActionSearcher = Callable[..., tuple[Sequence[GenericDisconnectionCandidate], Any]]


def _remaining_steps(node: Any) -> int:
    if node.reaction is None:
        return 0
    return 1 + max((_remaining_steps(child) for child in node.reaction.children), default=0)


def _patent_key(value: Optional[str]) -> Optional[str]:
    match = _PATENT_TOKEN.search(str(value or ""))
    if match is None:
        return None
    return f"US{int(match.group(1))}{match.group(2).upper()}"


def _distinct_ranks(values: Sequence[str]) -> list[int]:
    ranks: dict[str, int] = {}
    result = []
    for value in values:
        if value not in ranks:
            ranks[value] = len(ranks) + 1
        result.append(ranks[value])
    return result


def _identity_rank(values: Sequence[str], expected: Optional[str]) -> Optional[int]:
    if not expected:
        return None
    distinct = []
    seen = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        distinct.append(value)
    return distinct.index(expected) + 1 if expected in seen else None


def _outcome(
    exact_rank: Optional[int],
    strategy_rank: Optional[int],
    site_rank: Optional[int],
    operator_rank: Optional[int],
    candidate_count: int,
) -> str:
    if exact_rank == 1:
        return "exact_top1"
    if exact_rank is not None:
        return "exact_lower_rank"
    if strategy_rank is not None:
        return "strategy_equivalent_alternative"
    if site_rank is not None and operator_rank is not None:
        return "correct_site_and_operator"
    if site_rank is not None:
        return "correct_site_other_operator"
    if operator_rank is not None:
        return "correct_operator_other_site"
    return "other_valid_candidates" if candidate_count else "no_candidates"


def _candidate_actions(
    candidates: Sequence[GenericDisconnectionCandidate],
    *,
    expected_precursors: Optional[str],
    expected_site: Optional[str],
    expected_operator: Optional[str],
    expected_synthon: Optional[str],
    expected_strategy: Optional[str],
    patent_id: Optional[str],
) -> tuple[RouteActionCandidate, ...]:
    site_ranks = _distinct_ranks([item.disconnection_site_key for item in candidates])
    operator_ranks = _distinct_ranks([item.operator_id for item in candidates])
    synthon_ranks = _distinct_ranks([item.synthon_signature for item in candidates])
    strategy_ranks = _distinct_ranks([item.strategy_id for item in candidates])
    patent_key = _patent_key(patent_id)
    values = []
    for index, candidate in enumerate(candidates):
        overlaps = bool(
            patent_key
            and any(
                _patent_key(reaction_id) == patent_key
                for reaction_id in candidate.precedent_reaction_ids
            )
        )
        values.append(
            RouteActionCandidate(
                candidate_rank=index + 1,
                site_rank=site_ranks[index],
                operator_rank=operator_ranks[index],
                synthon_rank=synthon_ranks[index],
                strategy_rank=strategy_ranks[index],
                precursor_smiles=candidate.precursor_smiles,
                abstraction_level=candidate.abstraction_level,
                template_id=candidate.template_id,
                operator_id=candidate.operator_id,
                disconnection_site_key=candidate.disconnection_site_key,
                synthon_signature=candidate.synthon_signature,
                strategy_id=candidate.strategy_id,
                score=float(candidate.score),
                structural_score_band=candidate.structural_score_band,
                strategic_complexity_score=float(candidate.strategic_complexity_score),
                strategic_class=candidate.strategic_class,
                independent_reference_support=candidate.independent_reference_support,
                precursor_compatibility_disposition=(
                    candidate.precursor_compatibility_disposition
                ),
                exact_precursor_match=bool(
                    expected_precursors
                    and candidate.precursor_smiles == expected_precursors
                ),
                strategy_match=bool(
                    expected_strategy and candidate.strategy_id == expected_strategy
                ),
                site_match=bool(
                    expected_site
                    and candidate.disconnection_site_key == expected_site
                ),
                operator_match=bool(
                    expected_operator and candidate.operator_id == expected_operator
                ),
                synthon_match=bool(
                    expected_synthon
                    and candidate.synthon_signature == expected_synthon
                ),
                supervision_label=(
                    "observed_exact"
                    if (
                        expected_precursors
                        and candidate.precursor_smiles == expected_precursors
                    )
                    else "strategy_equivalent"
                    if expected_strategy and candidate.strategy_id == expected_strategy
                    else "same_site_operator"
                    if (
                        expected_site
                        and expected_operator
                        and
                        candidate.disconnection_site_key == expected_site
                        and candidate.operator_id == expected_operator
                    )
                    else "same_site"
                    if (
                        expected_site
                        and candidate.disconnection_site_key == expected_site
                    )
                    else "unchosen_alternative"
                ),
                source_patent_precedent_overlap=overlaps,
                precedent_reaction_ids=candidate.precedent_reaction_ids,
            )
        )
    return tuple(values)


def evaluate_route_actions(
    tree: ReactionRouteTree,
    library: GenericTemplateLibrary,
    *,
    config: RouteActionEvaluationConfig = RouteActionEvaluationConfig(),
    completion_prior_index: Optional[CompletionPriorIndex] = None,
    searcher: RouteActionSearcher = disconnect_operator_ladder_detailed,
) -> RouteActionEvaluation:
    """Replay every observed route step through the validated operator ladder."""

    prior_index = completion_prior_index
    if config.use_hierarchical_ranking and prior_index is None:
        prior_index = build_completion_prior_index(library)
    steps = []
    for molecule in iter_molecule_occurrences(tree):
        reaction = molecule.reaction
        if reaction is None:
            continue
        evaluation_id = digest(
            "RAE1", tree.tree_id, reaction.reaction_node_id, config.config_id
        )
        observed = build_observed_route_action_label(
            reaction.reaction_smiles,
            route_product_smiles=molecule.smiles,
            reaction_id=reaction.evidence.source_reaction_id or reaction.step_id,
            reference_id=tree.patent_id or tree.source_route_id or tree.tree_id,
        )
        raw_diagnostics: Any = {}
        candidates: tuple[GenericDisconnectionCandidate, ...] = ()
        if config.run_search and observed.search_eligible and observed.target_smiles:
            searched, raw_diagnostics = searcher(
                observed.target_smiles,
                library,
                top_k=config.top_k,
                max_templates_to_apply=config.max_templates_to_apply,
                max_candidates_to_validate=config.max_candidates_to_validate,
                use_context=config.use_context,
                include_l0=config.include_l0,
                diversify=config.diversify,
                use_hierarchical_ranking=config.use_hierarchical_ranking,
                minimum_candidates_per_level=config.minimum_candidates_per_level,
                lazy_validation=config.lazy_validation,
                completion_prior_index=prior_index,
            )
            candidates = tuple(searched)
        actions = _candidate_actions(
            candidates,
            expected_precursors=(
                observed.expected_precursor_smiles
                if observed.exact_precursors_verified
                else None
            ),
            expected_site=(
                observed.disconnection_site_key
                if observed.product_site_verified
                else None
            ),
            expected_operator=(
                observed.retained_operator_id
                if observed.retained_edits_verified
                else None
            ),
            expected_synthon=(
                observed.synthon_signature
                if observed.synthon_partition_verified
                else None
            ),
            expected_strategy=(
                observed.strategy_id if observed.strategy_verified else None
            ),
            patent_id=tree.patent_id,
        )
        exact_rank = next(
            (item.candidate_rank for item in actions if item.exact_precursor_match),
            None,
        )
        site_rank = _identity_rank(
            [item.disconnection_site_key for item in actions],
            observed.disconnection_site_key if observed.product_site_verified else None,
        )
        operator_rank = _identity_rank(
            [item.operator_id for item in actions],
            observed.retained_operator_id if observed.retained_edits_verified else None,
        )
        synthon_rank = _identity_rank(
            [item.synthon_signature for item in actions],
            observed.synthon_signature if observed.synthon_partition_verified else None,
        )
        strategy_rank = _identity_rank(
            [item.strategy_id for item in actions],
            observed.strategy_id if observed.strategy_verified else None,
        )
        diagnostics = (
            raw_diagnostics.to_dict()
            if isinstance(raw_diagnostics, OperatorLadderDiagnostics)
            else dict(raw_diagnostics or {})
        )
        steps.append(
            RouteActionStepEvaluation(
                evaluation_id=evaluation_id,
                reaction_node_id=reaction.reaction_node_id,
                step_id=reaction.step_id,
                source_reaction_id=reaction.evidence.source_reaction_id,
                product_occurrence_id=molecule.occurrence_id,
                retrosynthetic_depth=molecule.depth,
                observed_remaining_steps=_remaining_steps(molecule),
                reaction_smiles=reaction.reaction_smiles,
                route_target_smiles=tree.target_smiles,
                observed_action=observed,
                search_status=(
                    "searched"
                    if config.run_search and observed.search_eligible
                    else "not_run"
                    if observed.search_eligible
                    else "not_search_eligible"
                ),
                candidate_count=len(actions),
                exact_precursor_rank=exact_rank,
                site_rank=site_rank,
                operator_rank=operator_rank,
                synthon_rank=synthon_rank,
                strategy_rank=strategy_rank,
                outcome=_outcome(
                    exact_rank,
                    strategy_rank,
                    site_rank,
                    operator_rank,
                    len(actions),
                )
                if config.run_search and observed.search_eligible
                else "not_run"
                if observed.search_eligible
                else "not_search_eligible",
                source_patent_precedent_overlap=any(
                    item.source_patent_precedent_overlap for item in actions
                ),
                search_diagnostics=diagnostics,
                candidates=actions,
            )
        )
    route_evaluation_id = digest(
        "RAER1", tree.tree_id, config.config_id, *(step.evaluation_id for step in steps)
    )
    return RouteActionEvaluation(
        evaluation_id=route_evaluation_id,
        tree_id=tree.tree_id,
        source_route_id=tree.source_route_id,
        patent_id=tree.patent_id,
        split=tree.split,
        target_smiles=tree.target_smiles,
        reaction_count=tree.reaction_count,
        maximum_depth=tree.maximum_depth,
        search_config_id=config.config_id,
        steps=tuple(steps),
    )


__all__ = [
    "ROUTE_ACTION_EVALUATION_ALGORITHM_VERSION",
    "ROUTE_ACTION_EVALUATION_SCHEMA_VERSION",
    "RouteActionCandidate",
    "RouteActionEvaluation",
    "RouteActionEvaluationConfig",
    "RouteActionStepEvaluation",
    "evaluate_route_actions",
    "normalize_observed_reaction",
]
