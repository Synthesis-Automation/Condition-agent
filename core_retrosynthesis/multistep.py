"""Small deterministic multistep search over validated one-step operators."""

from __future__ import annotations

import heapq
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from typing import Any, Callable, Optional, Protocol

from cas_tools import PrecursorRealismAssessment
from reactive_taxonomy import assess_molecule_complexity
from .chemistry import digest

from .condition_ranking import RetrosynthesisConditionEvidence
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
    OperatorLadderDiagnostics,
)
from .generic_search import disconnect_operator_ladder_detailed
from .hierarchical_ranking import build_completion_prior_index
from .multistep_ranking import load_multistep_ranking_policy
from .route_tree import (
    CanonicalRouteTree,
    build_canonical_route_tree,
    route_distance_matrix,
)
from cas_tools.molecule_index import (
    MoleculeIndexMatch,
    MoleculeIdentity,
    molecule_identity,
)


MULTISTEP_SCHEMA_VERSION = "1.4"
_TERMINAL_STOCK_ROLES = frozenset(
    {
        "reactant",
        "starting_material",
        "startingmaterial",
        "substrate",
    }
)


class LiteratureLookup(Protocol):
    """Narrow terminal-catalog contract consumed by route search."""

    def lookup(
        self,
        identity: MoleculeIdentity,
        *,
        provenance_limit: int = 5,
    ) -> Optional[MoleculeIndexMatch]: ...


@dataclass(frozen=True)
class OneStepExpansionBatch:
    """Validated actions plus staged-expansion profiling."""

    candidates: tuple[GenericDisconnectionCandidate, ...]
    diagnostics: Optional[OperatorLadderDiagnostics] = None


OneStepExpander = Callable[
    [str, int],
    tuple[GenericDisconnectionCandidate, ...] | OneStepExpansionBatch,
]
ConditionEvidenceEvaluator = Callable[[str], RetrosynthesisConditionEvidence]
PrecursorRealismScorer = Callable[
    [str],
    tuple[PrecursorRealismAssessment, ...],
]


@dataclass(frozen=True)
class StartingMaterialAssessment:
    """Terminal decision for one route leaf."""

    smiles: str
    canonical_smiles: str
    depth: int
    molecular_weight: Optional[float]
    terminal: bool
    terminal_reasons: tuple[str, ...]
    terminal_evidence: str
    catalog_role_status: str
    unresolved_reason: Optional[str] = None
    literature_match: Optional[MoleculeIndexMatch] = None
    route_node_id: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible assessment."""

        return asdict(self)


@dataclass(frozen=True)
class RetrosynthesisRouteStep:
    """One validated reaction edge in a multistep route."""

    step_id: str
    depth: int
    product_smiles: str
    precursor_smiles: tuple[str, ...]
    step_cost: float
    step_cost_components: tuple[tuple[str, float], ...]
    candidate: GenericDisconnectionCandidate
    product_node_id: str = ""
    precursor_node_ids: tuple[str, ...] = ()
    condition_evidence: Optional[RetrosynthesisConditionEvidence] = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route step."""

        return {
            "step_id": self.step_id,
            "depth": self.depth,
            "product_smiles": self.product_smiles,
            "precursor_smiles": list(self.precursor_smiles),
            "step_cost": self.step_cost,
            "step_cost_components": dict(self.step_cost_components),
            "candidate": self.candidate.to_dict(),
            "product_node_id": self.product_node_id,
            "precursor_node_ids": list(self.precursor_node_ids),
            "condition_evidence": (
                self.condition_evidence.to_dict()
                if self.condition_evidence is not None
                else None
            ),
        }


@dataclass(frozen=True)
class MultistepRouteEvidenceSummary:
    """Route-level audit of step realism, compatibility, and conditions."""

    compatibility_warning_step_count: int
    strong_compatibility_warning_step_count: int
    realism_assessed_step_count: int
    weakest_precursor_realism_score: Optional[float]
    strategic_complexity_assessed_step_count: int
    strategic_step_count: int
    tactical_step_count: int
    complexity_increasing_step_count: int
    mean_strategic_complexity_reduction_score: Optional[float]
    target_to_frontier_complexity_reduction_fraction: Optional[float]
    condition_assessed_step_count: int
    condition_supported_step_count: int
    condition_direct_step_count: int
    condition_fallback_step_count: int
    condition_insufficient_step_count: int
    condition_coverage_fraction: Optional[float]
    condition_cost_penalty: float

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route evidence summary."""

        return asdict(self)


@dataclass(frozen=True)
class MultistepRetrosynthesisRoute:
    """One solved or bounded partial synthesis route."""

    route_id: str
    target_smiles: str
    solved: bool
    route_cost: float
    reaction_count: int
    maximum_depth: int
    steps: tuple[RetrosynthesisRouteStep, ...]
    leaves: tuple[StartingMaterialAssessment, ...]
    route_tree: CanonicalRouteTree
    evidence_summary: MultistepRouteEvidenceSummary
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route."""

        return {
            "route_id": self.route_id,
            "target_smiles": self.target_smiles,
            "solved": self.solved,
            "route_cost": self.route_cost,
            "reaction_count": self.reaction_count,
            "maximum_depth": self.maximum_depth,
            "steps": [step.to_dict() for step in self.steps],
            "leaves": [leaf.to_dict() for leaf in self.leaves],
            "route_tree": self.route_tree.to_dict(),
            "evidence_summary": self.evidence_summary.to_dict(),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class MultistepSearchDiagnostics:
    """Bounded-search counters for audit and tuning."""

    expanded_states: int
    one_step_calls: int
    one_step_cache_hits: int
    generated_candidates: int
    rejected_cycles: int
    rejected_invalid_candidates: int
    duplicate_states: int
    retained_alternative_paths: int
    frontier_states: int
    solved_routes_found: int
    partial_routes_found: int
    stopped_by_expansion_limit: bool
    proposed_actions: int = 0
    validation_attempts: int = 0
    valid_actions: int = 0
    first_solution_expansion: Optional[int] = None
    beam_pruned_states: int = 0
    dead_end_states: int = 0
    expansion_level_calls: tuple[tuple[str, int], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible diagnostics."""

        return asdict(self)


@dataclass(frozen=True)
class MultistepRetrosynthesisResult:
    """Versioned solved routes, useful partial routes, and diagnostics."""

    target_smiles: str
    routes: tuple[MultistepRetrosynthesisRoute, ...]
    partial_routes: tuple[MultistepRetrosynthesisRoute, ...]
    diagnostics: MultistepSearchDiagnostics
    max_depth: int
    molecular_weight_threshold: float
    ranking_policy_definition_id: str
    precursor_realism_enabled: bool = False
    condition_availability_enabled: bool = False
    schema_version: str = MULTISTEP_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible result."""

        warnings = [
            "Literature matches mean observed in the configured corpus, "
            "not verified commercial availability.",
            "Molecular-weight terminals are a search heuristic and still "
            "require chemist review.",
            "Literature matches are terminal only when retained provenance "
            "identifies a reactant-like source role; rebuild legacy indexes "
            "with source_role provenance.",
        ]
        if self.precursor_realism_enabled:
            warnings.append(
                "Precursor realism is a versioned evidence heuristic, not a "
                "probability or a hard feasibility filter."
            )
        if self.condition_availability_enabled:
            warnings.append(
                "Condition availability reports precedent support for each "
                "retained reaction; it does not guarantee route execution."
            )
        return {
            "target_smiles": self.target_smiles,
            "routes": [route.to_dict() for route in self.routes],
            "partial_routes": [route.to_dict() for route in self.partial_routes],
            "diagnostics": self.diagnostics.to_dict(),
            "max_depth": self.max_depth,
            "molecular_weight_threshold": self.molecular_weight_threshold,
            "ranking_policy_definition_id": (self.ranking_policy_definition_id),
            "precursor_realism_enabled": self.precursor_realism_enabled,
            "condition_availability_enabled": self.condition_availability_enabled,
            "schema_version": self.schema_version,
            "route_postprocessing": {
                "definition_id": "route_distance_multiset_jaccard.v1",
                "tree_ids": [route.route_tree.tree_id for route in self.routes],
                "distance_matrix": [
                    list(row)
                    for row in route_distance_matrix(
                        route.route_tree for route in self.routes
                    )
                ],
            },
            "warnings": warnings,
        }


@dataclass(frozen=True)
class _Leaf:
    assessment: StartingMaterialAssessment
    ancestors: tuple[str, ...]


@dataclass(frozen=True)
class _RouteState:
    steps: tuple[RetrosynthesisRouteStep, ...]
    leaves: tuple[_Leaf, ...]
    cost: float


def _catalog_role_status(
    literature_match: MoleculeIndexMatch | None,
) -> str:
    if literature_match is None:
        return "no_match"
    roles = {
        str(
            record.get("source_role")
            or record.get("contextual_role")
            or record.get("role")
            or ""
        )
        .strip()
        .casefold()
        .replace("-", "_")
        .replace(" ", "_")
        for record in literature_match.source_records
    }
    roles.discard("")
    if roles & _TERMINAL_STOCK_ROLES:
        return "reactant_supported"
    if roles:
        return "nonreactant_only"
    return "untyped"


def _has_strong_stock_evidence(
    literature_match: MoleculeIndexMatch | None,
) -> bool:
    """Return whether retained provenance explicitly confirms usable stock."""

    if literature_match is None:
        return False
    return any(
        str(record.get("stock_evidence") or "").strip().casefold()
        in {"physically_available", "supplier_in_stock"}
        and str(record.get("terminal_eligible") or "").strip().casefold()
        in {"1", "true", "yes"}
        for record in literature_match.source_records
    )


def _assess_starting_material(
    smiles: str,
    *,
    depth: int,
    literature_index: LiteratureLookup,
    molecular_weight_threshold: float,
    allow_terminal: bool = True,
    allow_untyped_literature_terminal: bool = False,
    unresolved_reason: str | None = None,
) -> StartingMaterialAssessment:
    identity = molecule_identity(smiles)
    if identity is None:
        return StartingMaterialAssessment(
            smiles=smiles,
            canonical_smiles=smiles,
            depth=depth,
            molecular_weight=None,
            terminal=False,
            terminal_reasons=(),
            terminal_evidence="none",
            catalog_role_status="no_match",
            unresolved_reason=unresolved_reason or "invalid_molecule",
        )
    literature_match = literature_index.lookup(identity)
    catalog_role_status = _catalog_role_status(literature_match)
    strong_stock_evidence = _has_strong_stock_evidence(literature_match)
    reasons = []
    if identity.molecular_weight <= molecular_weight_threshold:
        reasons.append("molecular_weight_threshold")
    if strong_stock_evidence:
        reasons.append("supplier_stock_match")
    elif catalog_role_status == "reactant_supported":
        reasons.append("literature_reactant_match")
    elif catalog_role_status == "untyped" and allow_untyped_literature_terminal:
        reasons.append("literature_untyped_match")
    terminal = allow_terminal and bool(reasons)
    terminal_evidence = "none"
    if terminal:
        terminal_evidence = (
            "supplier_stock_portfolio"
            if "supplier_stock_match" in reasons
            else "role_supported_literature"
            if "literature_reactant_match" in reasons
            else "untyped_literature"
            if "literature_untyped_match" in reasons
            else "heuristic_molecular_weight"
        )
    return StartingMaterialAssessment(
        smiles=smiles,
        canonical_smiles=identity.canonical_smiles,
        depth=depth,
        molecular_weight=identity.molecular_weight,
        terminal=terminal,
        terminal_reasons=tuple(reasons) if allow_terminal else (),
        terminal_evidence=terminal_evidence,
        catalog_role_status=catalog_role_status,
        unresolved_reason=(None if terminal else unresolved_reason),
        literature_match=literature_match,
    )


def _state_signature(state: _RouteState) -> tuple[tuple[Any, ...], ...]:
    return tuple(
        sorted(
            (
                leaf.assessment.canonical_smiles,
                leaf.assessment.depth,
                leaf.assessment.terminal,
                leaf.assessment.unresolved_reason or "",
            )
            for leaf in state.leaves
        )
    )


def _expandable_leaves(
    state: _RouteState,
    max_depth: int,
) -> tuple[tuple[int, _Leaf], ...]:
    return tuple(
        (index, leaf)
        for index, leaf in enumerate(state.leaves)
        if not leaf.assessment.terminal
        and leaf.assessment.unresolved_reason is None
        and leaf.assessment.depth < max_depth
    )


def _state_priority(
    state: _RouteState,
    max_depth: int,
) -> tuple[Any, ...]:
    blocked = sum(
        not leaf.assessment.terminal
        and (
            leaf.assessment.unresolved_reason is not None
            or leaf.assessment.depth >= max_depth
        )
        for leaf in state.leaves
    )
    return (
        bool(blocked),
        round(state.cost, 8),
        len(_expandable_leaves(state, max_depth)),
        len(state.steps),
        _state_signature(state),
    )


def _select_leaf(
    expandable: tuple[tuple[int, _Leaf], ...],
) -> tuple[int, _Leaf]:
    return min(
        expandable,
        key=lambda item: (
            -(item[1].assessment.molecular_weight or 0.0),
            item[1].assessment.canonical_smiles,
            item[0],
        ),
    )


def _route_history_signature(state: _RouteState) -> tuple[str, ...]:
    return tuple(step.candidate.proposed_reaction_smiles for step in state.steps)


def _admit_state_path(
    state: _RouteState,
    retained: dict[
        tuple[tuple[Any, ...], ...],
        list[tuple[float, tuple[str, ...]]],
    ],
    maximum_paths: int,
) -> tuple[bool, bool]:
    signature = _state_signature(state)
    history = _route_history_signature(state)
    bucket = retained.setdefault(signature, [])
    for index, (cost, existing_history) in enumerate(bucket):
        if existing_history != history:
            continue
        if cost <= state.cost:
            return False, False
        bucket[index] = (state.cost, history)
        bucket.sort()
        return True, False
    alternative = bool(bucket)
    if len(bucket) < maximum_paths:
        bucket.append((state.cost, history))
        bucket.sort()
        return True, alternative
    worst = max(bucket)
    candidate = (state.cost, history)
    if candidate >= worst:
        return False, False
    bucket.remove(worst)
    bucket.append(candidate)
    bucket.sort()
    return True, alternative


def _is_retained_state_path(
    state: _RouteState,
    retained: dict[
        tuple[tuple[Any, ...], ...],
        list[tuple[float, tuple[str, ...]]],
    ],
) -> bool:
    return (
        state.cost,
        _route_history_signature(state),
    ) in retained.get(_state_signature(state), ())


def _step_cost(
    candidate: GenericDisconnectionCandidate,
    candidate_rank: int,
    child_leaves: tuple[_Leaf, ...],
) -> tuple[float, tuple[tuple[str, float], ...]]:
    policy = load_multistep_ranking_policy()
    components = {
        "base": policy.base_step_cost,
        "candidate_score_deficit": (
            policy.candidate_score_deficit_weight
            * max(0.0, 1.0 - float(candidate.score))
        ),
        "abstraction_level": policy.abstraction_penalty(candidate.abstraction_level),
        "selectivity_warnings": (
            policy.selectivity_warning_penalty * len(candidate.selectivity_warnings)
        ),
        "heuristic_terminals": (
            policy.heuristic_terminal_penalty
            * sum(
                leaf.assessment.terminal_evidence == "heuristic_molecular_weight"
                for leaf in child_leaves
            )
        ),
        "precursor_compatibility": (
            policy.precursor_compatibility_band_penalty_weight
            * candidate.precursor_compatibility_band_penalty
        ),
        "precursor_realism": (
            policy.precursor_realism_band_penalty_weight
            * candidate.precursor_realism_band_penalty
        ),
        "strategic_progress_deficit": (
            (
                policy.strategic_progress_deficit_penalty_weight
                * max(
                    0.0,
                    policy.strategic_progress_minimum_score
                    - candidate.strategic_complexity_score,
                )
                / max(policy.strategic_progress_minimum_score, 1e-12)
            )
            if candidate.strategic_complexity is not None
            else 0.0
        ),
        "complexity_increasing_step": (
            policy.complexity_increasing_step_penalty
            if candidate.strategic_class == "complexity_increasing"
            else 0.0
        ),
        "candidate_rank_tiebreak": (policy.candidate_rank_tiebreak * candidate_rank),
    }
    normalized = tuple((name, round(value, 8)) for name, value in components.items())
    return round(sum(value for _, value in normalized), 8), normalized


def _route_evidence_summary(
    steps: tuple[RetrosynthesisRouteStep, ...],
    *,
    target_smiles: str = "",
    leaves: tuple[StartingMaterialAssessment, ...] = (),
) -> MultistepRouteEvidenceSummary:
    """Aggregate per-step evidence without hiding the underlying traces."""

    realism_scores = tuple(
        float(step.candidate.precursor_realism_score)
        for step in steps
        if step.candidate.precursor_realism_score is not None
    )
    condition_evidence = tuple(
        step.condition_evidence
        for step in steps
        if step.condition_evidence is not None
    )
    direct = sum(
        evidence.status == "recommended_direct" for evidence in condition_evidence
    )
    fallback = sum(
        evidence.status == "recommended_fallback" for evidence in condition_evidence
    )
    insufficient = sum(
        evidence.status == "insufficient_evidence"
        for evidence in condition_evidence
    )
    condition_penalty = sum(
        dict(step.step_cost_components).get("condition_availability", 0.0)
        for step in steps
    )
    strategic_assessments = tuple(
        step.candidate.strategic_complexity
        for step in steps
        if step.candidate.strategic_complexity is not None
    )
    frontier_reduction = None
    if target_smiles and leaves:
        try:
            target_complexity = assess_molecule_complexity(
                target_smiles
            ).raw_complexity
            leaf_complexities = sorted(
                (
                    assess_molecule_complexity(
                        leaf.canonical_smiles
                    ).raw_complexity
                    for leaf in leaves
                ),
                reverse=True,
            )
            frontier_complexity = (
                leaf_complexities[0] + 0.15 * sum(leaf_complexities[1:])
                if leaf_complexities
                else target_complexity
            )
            frontier_reduction = round(
                (target_complexity - frontier_complexity)
                / max(target_complexity, frontier_complexity, 1e-12),
                8,
            )
        except ValueError:
            frontier_reduction = None
    return MultistepRouteEvidenceSummary(
        compatibility_warning_step_count=sum(
            step.candidate.precursor_compatibility_disposition != "pass"
            for step in steps
        ),
        strong_compatibility_warning_step_count=sum(
            step.candidate.precursor_compatibility_warning_strength == "strong"
            for step in steps
        ),
        realism_assessed_step_count=len(realism_scores),
        weakest_precursor_realism_score=(
            round(min(realism_scores), 8) if realism_scores else None
        ),
        strategic_complexity_assessed_step_count=len(strategic_assessments),
        strategic_step_count=sum(
            assessment.is_strategic for assessment in strategic_assessments
        ),
        tactical_step_count=sum(
            not assessment.is_strategic for assessment in strategic_assessments
        ),
        complexity_increasing_step_count=sum(
            assessment.strategic_class == "complexity_increasing"
            for assessment in strategic_assessments
        ),
        mean_strategic_complexity_reduction_score=(
            round(
                sum(assessment.score for assessment in strategic_assessments)
                / len(strategic_assessments),
                8,
            )
            if strategic_assessments
            else None
        ),
        target_to_frontier_complexity_reduction_fraction=frontier_reduction,
        condition_assessed_step_count=len(condition_evidence),
        condition_supported_step_count=direct + fallback,
        condition_direct_step_count=direct,
        condition_fallback_step_count=fallback,
        condition_insufficient_step_count=insufficient,
        condition_coverage_fraction=(
            round((direct + fallback) / len(condition_evidence), 8)
            if condition_evidence
            else None
        ),
        condition_cost_penalty=round(condition_penalty, 8),
    )


def _route_warnings(
    steps: tuple[RetrosynthesisRouteStep, ...],
    existing: tuple[str, ...],
) -> tuple[str, ...]:
    warnings = list(existing)
    if any(
        step.candidate.precursor_compatibility_warning_strength == "strong"
        for step in steps
    ):
        warnings.append("STRONG_INTRAMOLECULAR_COMPATIBILITY_WARNING")
    if any(
        step.candidate.precursor_realism_band_penalty > 0 for step in steps
    ):
        warnings.append("LOW_PRECURSOR_REALISM")
    assessed = tuple(
        step.candidate.strategic_complexity
        for step in steps
        if step.candidate.strategic_complexity is not None
    )
    if assessed and not any(value.is_strategic for value in assessed):
        warnings.append("NO_STRATEGIC_COMPLEXITY_REDUCTION")
    if any(
        value.strategic_class == "complexity_increasing" for value in assessed
    ):
        warnings.append("COMPLEXITY_INCREASING_RETRO_STEP")
    evidence = tuple(
        step.condition_evidence
        for step in steps
        if step.condition_evidence is not None
    )
    if any(value.status == "recommended_fallback" for value in evidence):
        warnings.append("CONDITION_FALLBACK_USED")
    if any(value.status == "insufficient_evidence" for value in evidence):
        warnings.append("CONDITION_EVIDENCE_INCOMPLETE")
    return tuple(dict.fromkeys(warnings))


def _attach_route_condition_evidence(
    route: MultistepRetrosynthesisRoute,
    evaluator: ConditionEvidenceEvaluator,
    cache: dict[str, RetrosynthesisConditionEvidence],
) -> MultistepRetrosynthesisRoute:
    """Attach cached condition evidence and its declarative route penalty."""

    policy = load_multistep_ranking_policy()
    steps = []
    route_penalty = 0.0
    for step in route.steps:
        query = (
            step.candidate.condition_query_reaction_smiles
            or step.candidate.proposed_reaction_smiles
        )
        evidence = cache.get(query)
        if evidence is None:
            evidence = evaluator(query)
            cache[query] = evidence
        penalty = round(policy.condition_status_penalty(evidence.status), 8)
        route_penalty += penalty
        components = tuple(
            (*step.step_cost_components, ("condition_availability", penalty))
        )
        steps.append(
            replace(
                step,
                step_cost=round(step.step_cost + penalty, 8),
                step_cost_components=components,
                condition_evidence=evidence,
            )
        )
    normalized_steps = tuple(steps)
    return replace(
        route,
        route_cost=round(route.route_cost + route_penalty, 8),
        steps=normalized_steps,
        evidence_summary=_route_evidence_summary(
            normalized_steps,
            target_smiles=route.target_smiles,
            leaves=route.leaves,
        ),
        warnings=_route_warnings(normalized_steps, route.warnings),
    )


def _route_from_state(
    target_smiles: str,
    state: _RouteState,
    *,
    max_depth: int,
) -> MultistepRetrosynthesisRoute:
    solved = all(leaf.assessment.terminal for leaf in state.leaves)
    step_tokens = tuple(step.candidate.proposed_reaction_smiles for step in state.steps)
    route_id = digest("ROUTE1", target_smiles, *step_tokens)
    warnings = []
    if not solved:
        warnings.append("ROUTE_INCOMPLETE")
    if any(
        leaf.assessment.unresolved_reason == "maximum_depth"
        or (not leaf.assessment.terminal and leaf.assessment.depth >= max_depth)
        for leaf in state.leaves
    ):
        warnings.append("MAXIMUM_DEPTH_REACHED")
    if any(
        leaf.assessment.terminal_evidence == "heuristic_molecular_weight"
        for leaf in state.leaves
    ):
        warnings.append("HEURISTIC_TERMINAL_MATERIALS")
    if any(
        leaf.assessment.catalog_role_status in {"nonreactant_only", "untyped"}
        for leaf in state.leaves
    ):
        warnings.append("CATALOG_MATCH_NOT_REACTANT_SUPPORTED")
    route_leaves = tuple(
        sorted(
            (leaf.assessment for leaf in state.leaves),
            key=lambda value: (
                value.canonical_smiles,
                value.depth,
                value.route_node_id,
                value.unresolved_reason or "",
            ),
        )
    )
    root_node_id = next(
        (
            step.product_node_id
            for step in state.steps
            if step.depth == 1
        ),
        next(
            (
                leaf.route_node_id
                for leaf in route_leaves
                if leaf.depth == 0
            ),
            digest("RMOL1", target_smiles, "root"),
        ),
    )
    route_tree = build_canonical_route_tree(
        target_smiles,
        root_node_id,
        state.steps,
        route_leaves,
    )
    return MultistepRetrosynthesisRoute(
        route_id=route_id,
        target_smiles=target_smiles,
        solved=solved,
        route_cost=round(state.cost, 8),
        reaction_count=len(state.steps),
        maximum_depth=max(
            (step.depth for step in state.steps),
            default=0,
        ),
        steps=state.steps,
        leaves=route_leaves,
        route_tree=route_tree,
        evidence_summary=_route_evidence_summary(
            state.steps,
            target_smiles=target_smiles,
            leaves=route_leaves,
        ),
        warnings=_route_warnings(state.steps, tuple(warnings)),
    )


def plan_multistep_routes(
    target_smiles: str,
    library: GenericTemplateLibrary,
    literature_index: LiteratureLookup,
    *,
    max_depth: int = 3,
    molecular_weight_threshold: float = 150.0,
    top_k_routes: int = 5,
    per_step_top_k: int = 5,
    beam_width: int = 20,
    max_expansions: int = 100,
    max_templates_to_apply: int = 300,
    max_candidates_to_validate: int = 50,
    use_context: bool = True,
    include_l0: bool = True,
    diversify: bool = True,
    use_hierarchical_ranking: bool = True,
    allow_untyped_literature_terminals: bool = False,
    max_paths_per_state: int | None = None,
    minimum_candidates_per_level: int | None = None,
    precursor_realism_scorer: PrecursorRealismScorer | None = None,
    condition_evidence_evaluator: ConditionEvidenceEvaluator | None = None,
    expander: OneStepExpander | None = None,
) -> MultistepRetrosynthesisResult:
    """Find short routes whose leaves pass the explicit terminal predicate."""

    if max_depth not in {2, 3}:
        raise ValueError("maximum route depth must be 2 or 3")
    if molecular_weight_threshold <= 0:
        raise ValueError("molecular-weight threshold must be positive")
    for value, name in (
        (top_k_routes, "top-k routes"),
        (per_step_top_k, "per-step top-k"),
        (beam_width, "beam width"),
        (max_expansions, "maximum expansions"),
    ):
        if value < 1:
            raise ValueError(f"{name} must be positive")
    target = molecule_identity(target_smiles)
    if target is None:
        raise ValueError("invalid target molecule")
    ranking_policy = load_multistep_ranking_policy()
    active_max_paths = (
        ranking_policy.maximum_paths_per_state
        if max_paths_per_state is None
        else max_paths_per_state
    )
    active_minimum_per_level = (
        ranking_policy.minimum_candidates_per_level
        if minimum_candidates_per_level is None
        else minimum_candidates_per_level
    )
    if active_max_paths < 1:
        raise ValueError("maximum paths per state must be positive")
    if active_minimum_per_level < 0:
        raise ValueError("minimum candidates per level cannot be negative")
    completion_prior_index = (
        build_completion_prior_index(library)
        if diversify and use_hierarchical_ranking
        else None
    )
    active_realism_scorer: PrecursorRealismScorer | None = None
    if precursor_realism_scorer is not None:
        uncached_realism_scorer = precursor_realism_scorer

        @lru_cache(maxsize=20_000)
        def cached_realism_scorer(
            precursor_smiles: str,
        ) -> tuple[PrecursorRealismAssessment, ...]:
            return uncached_realism_scorer(precursor_smiles)

        active_realism_scorer = cached_realism_scorer

    def default_expander(
        product_smiles: str,
        top_k: int,
    ) -> OneStepExpansionBatch:
        candidates, diagnostics = disconnect_operator_ladder_detailed(
            product_smiles,
            library,
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
            include_l0=include_l0,
            diversify=diversify,
            use_hierarchical_ranking=use_hierarchical_ranking,
            minimum_candidates_per_level=active_minimum_per_level,
            lazy_validation=True,
            precursor_realism_scorer=active_realism_scorer,
            completion_prior_index=completion_prior_index,
        )
        return OneStepExpansionBatch(candidates, diagnostics)

    active_expander = expander or default_expander
    root_node_id = digest("RMOL1", target.canonical_smiles, "root")
    root = _Leaf(
        assessment=replace(
            _assess_starting_material(
            target.canonical_smiles,
            depth=0,
            literature_index=literature_index,
            molecular_weight_threshold=molecular_weight_threshold,
            allow_terminal=False,
            allow_untyped_literature_terminal=(allow_untyped_literature_terminals),
            ),
            route_node_id=root_node_id,
        ),
        ancestors=(),
    )
    initial = _RouteState(steps=(), leaves=(root,), cost=0.0)
    queue: list[tuple[tuple[Any, ...], int, _RouteState]] = []
    serial = 0
    heapq.heappush(queue, (_state_priority(initial, max_depth), serial, initial))
    best_paths_by_state = {
        _state_signature(initial): [(0.0, _route_history_signature(initial))]
    }
    expansion_cache: dict[
        str,
        tuple[GenericDisconnectionCandidate, ...],
    ] = {}
    solved_states: dict[str, _RouteState] = {}
    partial_states: dict[
        tuple[tuple[tuple[Any, ...], ...], tuple[str, ...]],
        _RouteState,
    ] = {}
    expanded_states = 0
    one_step_calls = 0
    cache_hits = 0
    generated_candidates = 0
    rejected_cycles = 0
    rejected_invalid = 0
    duplicate_states = 0
    retained_alternative_paths = 0
    proposed_actions = 0
    validation_attempts = 0
    valid_actions = 0
    first_solution_expansion: int | None = None
    beam_pruned_states = 0
    dead_end_states = 0
    expansion_level_calls: dict[str, int] = {}

    while queue and expanded_states < max_expansions:
        _, _, state = heapq.heappop(queue)
        if not _is_retained_state_path(state, best_paths_by_state):
            continue
        if all(leaf.assessment.terminal for leaf in state.leaves):
            route = _route_from_state(
                target.canonical_smiles,
                state,
                max_depth=max_depth,
            )
            solved_states.setdefault(route.route_id, state)
            continue
        expandable = _expandable_leaves(state, max_depth)
        if not expandable:
            partial_states.setdefault(
                (_state_signature(state), _route_history_signature(state)),
                state,
            )
            continue
        leaf_index, leaf = _select_leaf(expandable)
        product = leaf.assessment.canonical_smiles
        candidates = expansion_cache.get(product)
        if candidates is None:
            expansion = active_expander(product, per_step_top_k)
            if isinstance(expansion, OneStepExpansionBatch):
                candidates = expansion.candidates
                expansion_diagnostics = expansion.diagnostics
                if expansion_diagnostics is not None:
                    proposed_actions += expansion_diagnostics.proposed_action_count
                    validation_attempts += (
                        expansion_diagnostics.validation_attempt_count
                    )
                    valid_actions += expansion_diagnostics.valid_action_count
                    for level in expansion_diagnostics.levels_attempted:
                        expansion_level_calls[level] = (
                            expansion_level_calls.get(level, 0) + 1
                        )
                else:
                    proposed_actions += len(candidates)
                    validation_attempts += len(candidates)
                    valid_actions += len(candidates)
            else:
                candidates = expansion
                proposed_actions += len(candidates)
                validation_attempts += len(candidates)
                valid_actions += sum(
                    candidate.forward_validation_status == "verified_signature"
                    for candidate in candidates
                )
            expansion_cache[product] = candidates
            one_step_calls += 1
        else:
            cache_hits += 1
        expanded_states += 1
        generated_candidates += len(candidates)
        if not candidates:
            dead_end_states += 1
            stopped = _Leaf(
                assessment=replace(
                    _assess_starting_material(
                        product,
                        depth=leaf.assessment.depth,
                        literature_index=literature_index,
                        molecular_weight_threshold=molecular_weight_threshold,
                        allow_untyped_literature_terminal=(
                            allow_untyped_literature_terminals
                        ),
                        unresolved_reason="no_candidates",
                    ),
                    route_node_id=leaf.assessment.route_node_id,
                ),
                ancestors=leaf.ancestors,
            )
            next_leaves = list(state.leaves)
            next_leaves[leaf_index] = stopped
            partial = _RouteState(
                steps=state.steps,
                leaves=tuple(next_leaves),
                cost=state.cost,
            )
            partial_states[
                (_state_signature(partial), _route_history_signature(partial))
            ] = partial
            continue

        for candidate_rank, candidate in enumerate(candidates, start=1):
            if candidate.forward_validation_status != "verified_signature":
                rejected_invalid += 1
                continue
            canonical_precursors = molecule_identity(candidate.precursor_smiles)
            if canonical_precursors is None:
                rejected_invalid += 1
                continue
            precursor_values = tuple(canonical_precursors.canonical_smiles.split("."))
            ancestor_set = set(leaf.ancestors) | {product}
            if any(value in ancestor_set for value in precursor_values):
                rejected_cycles += 1
                continue
            child_depth = leaf.assessment.depth + 1
            child_leaves = []
            precursor_node_ids = []
            precursor_occurrences: dict[str, int] = {}
            step_id = digest(
                "RSTEP1",
                product,
                candidate.proposed_reaction_smiles,
                str(child_depth),
                leaf.assessment.route_node_id,
            )
            for precursor in precursor_values:
                occurrence = precursor_occurrences.get(precursor, 0)
                precursor_occurrences[precursor] = occurrence + 1
                precursor_node_id = digest(
                    "RMOL1",
                    leaf.assessment.route_node_id,
                    step_id,
                    precursor,
                    str(occurrence),
                )
                precursor_node_ids.append(precursor_node_id)
                assessment = _assess_starting_material(
                    precursor,
                    depth=child_depth,
                    literature_index=literature_index,
                    molecular_weight_threshold=molecular_weight_threshold,
                    allow_untyped_literature_terminal=(
                        allow_untyped_literature_terminals
                    ),
                    unresolved_reason=(
                        "maximum_depth" if child_depth >= max_depth else None
                    ),
                )
                child_leaves.append(
                    _Leaf(
                        assessment=replace(
                            assessment,
                            route_node_id=precursor_node_id,
                        ),
                        ancestors=(*leaf.ancestors, product),
                    )
                )
            normalized_child_leaves = tuple(child_leaves)
            step_cost, step_cost_components = _step_cost(
                candidate,
                candidate_rank,
                normalized_child_leaves,
            )
            step = RetrosynthesisRouteStep(
                step_id=step_id,
                depth=child_depth,
                product_smiles=product,
                precursor_smiles=precursor_values,
                step_cost=step_cost,
                step_cost_components=step_cost_components,
                candidate=candidate,
                product_node_id=leaf.assessment.route_node_id,
                precursor_node_ids=tuple(precursor_node_ids),
            )
            next_leaves = list(state.leaves)
            next_leaves[leaf_index : leaf_index + 1] = child_leaves
            next_state = _RouteState(
                steps=(*state.steps, step),
                leaves=tuple(next_leaves),
                cost=round(
                    state.cost + step_cost,
                    8,
                ),
            )
            signature = _state_signature(next_state)
            admitted, alternative = _admit_state_path(
                next_state,
                best_paths_by_state,
                active_max_paths,
            )
            if not admitted:
                duplicate_states += 1
                continue
            if alternative:
                retained_alternative_paths += 1
            if all(value.assessment.terminal for value in next_state.leaves):
                route = _route_from_state(
                    target.canonical_smiles,
                    next_state,
                    max_depth=max_depth,
                )
                solved_states.setdefault(route.route_id, next_state)
                if first_solution_expansion is None:
                    first_solution_expansion = expanded_states
                continue
            if not _expandable_leaves(next_state, max_depth):
                partial_states[(signature, _route_history_signature(next_state))] = (
                    next_state
                )
                continue
            serial += 1
            heapq.heappush(
                queue,
                (_state_priority(next_state, max_depth), serial, next_state),
            )
        if len(queue) > beam_width:
            queue_size_before_pruning = len(queue)
            queue = [
                item
                for item in queue
                if _is_retained_state_path(item[2], best_paths_by_state)
            ]
            queue = heapq.nsmallest(beam_width, queue)
            beam_pruned_states += max(0, queue_size_before_pruning - len(queue))
            heapq.heapify(queue)

    stopped_by_limit = bool(queue and expanded_states >= max_expansions)
    if stopped_by_limit:
        retained_frontier = (
            item
            for item in queue
            if _is_retained_state_path(item[2], best_paths_by_state)
        )
        for _, _, state in heapq.nsmallest(
            top_k_routes,
            retained_frontier,
        ):
            leaves = tuple(
                _Leaf(
                    assessment=(
                        leaf.assessment
                        if leaf.assessment.terminal
                        else replace(
                            _assess_starting_material(
                                leaf.assessment.canonical_smiles,
                                depth=leaf.assessment.depth,
                                literature_index=literature_index,
                                molecular_weight_threshold=(
                                    molecular_weight_threshold
                                ),
                                allow_untyped_literature_terminal=(
                                    allow_untyped_literature_terminals
                                ),
                                unresolved_reason="search_limit",
                            ),
                            route_node_id=leaf.assessment.route_node_id,
                        )
                    ),
                    ancestors=leaf.ancestors,
                )
                for leaf in state.leaves
            )
            limited = _RouteState(state.steps, leaves, state.cost)
            partial_states.setdefault(
                (
                    _state_signature(limited),
                    _route_history_signature(limited),
                ),
                limited,
            )

    solved_route_values = tuple(
        _route_from_state(
            target.canonical_smiles,
            state,
            max_depth=max_depth,
        )
        for state in solved_states.values()
    )
    partial_route_values = tuple(
        _route_from_state(
            target.canonical_smiles,
            state,
            max_depth=max_depth,
        )
        for state in partial_states.values()
    )
    if condition_evidence_evaluator is not None:
        condition_cache: dict[str, RetrosynthesisConditionEvidence] = {}
        solved_route_values = tuple(
            _attach_route_condition_evidence(
                route,
                condition_evidence_evaluator,
                condition_cache,
            )
            for route in solved_route_values
        )
        partial_route_values = tuple(
            _attach_route_condition_evidence(
                route,
                condition_evidence_evaluator,
                condition_cache,
            )
            for route in partial_route_values
        )
    solved_routes = tuple(
        sorted(
            solved_route_values,
            key=lambda route: (
                route.route_cost,
                route.reaction_count,
                route.route_id,
            ),
        )[:top_k_routes]
    )
    partial_routes = tuple(
        sorted(
            partial_route_values,
            key=lambda route: (
                sum(not leaf.terminal for leaf in route.leaves),
                route.route_cost,
                route.reaction_count,
                route.route_id,
            ),
        )[:top_k_routes]
    )
    diagnostics = MultistepSearchDiagnostics(
        expanded_states=expanded_states,
        one_step_calls=one_step_calls,
        one_step_cache_hits=cache_hits,
        generated_candidates=generated_candidates,
        rejected_cycles=rejected_cycles,
        rejected_invalid_candidates=rejected_invalid,
        duplicate_states=duplicate_states,
        retained_alternative_paths=retained_alternative_paths,
        frontier_states=len(queue),
        solved_routes_found=len(solved_states),
        partial_routes_found=len(partial_states),
        stopped_by_expansion_limit=stopped_by_limit,
        proposed_actions=proposed_actions,
        validation_attempts=validation_attempts,
        valid_actions=valid_actions,
        first_solution_expansion=first_solution_expansion,
        beam_pruned_states=beam_pruned_states,
        dead_end_states=dead_end_states,
        expansion_level_calls=tuple(sorted(expansion_level_calls.items())),
    )
    return MultistepRetrosynthesisResult(
        target_smiles=target.canonical_smiles,
        routes=solved_routes,
        partial_routes=partial_routes,
        diagnostics=diagnostics,
        max_depth=max_depth,
        molecular_weight_threshold=molecular_weight_threshold,
        ranking_policy_definition_id=ranking_policy.definition_id,
        precursor_realism_enabled=precursor_realism_scorer is not None,
        condition_availability_enabled=condition_evidence_evaluator is not None,
    )


__all__ = [
    "MULTISTEP_SCHEMA_VERSION",
    "LiteratureLookup",
    "MultistepRouteEvidenceSummary",
    "MultistepRetrosynthesisResult",
    "MultistepRetrosynthesisRoute",
    "MultistepSearchDiagnostics",
    "OneStepExpansionBatch",
    "RetrosynthesisRouteStep",
    "StartingMaterialAssessment",
    "plan_multistep_routes",
]
