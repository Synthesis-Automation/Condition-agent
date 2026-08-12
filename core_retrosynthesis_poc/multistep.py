"""Small deterministic multistep search over validated one-step operators."""

from __future__ import annotations

import heapq
from dataclasses import asdict, dataclass
from typing import Any, Callable, Optional, Protocol

from retrosynthesis_poc.chemistry import digest

from .generic_models import (
    GenericDisconnectionCandidate,
    GenericTemplateLibrary,
)
from .generic_search import disconnect_operator_ladder
from cas_tools.molecule_index import (
    MoleculeIndexMatch,
    MoleculeIdentity,
    molecule_identity,
)


MULTISTEP_SCHEMA_VERSION = "1.0"


class LiteratureLookup(Protocol):
    """Narrow terminal-catalog contract consumed by route search."""

    def lookup(
        self,
        identity: MoleculeIdentity,
        *,
        provenance_limit: int = 5,
    ) -> Optional[MoleculeIndexMatch]: ...


OneStepExpander = Callable[
    [str, int],
    tuple[GenericDisconnectionCandidate, ...],
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
    unresolved_reason: Optional[str] = None
    literature_match: Optional[MoleculeIndexMatch] = None

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
    candidate: GenericDisconnectionCandidate

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route step."""

        return {
            "step_id": self.step_id,
            "depth": self.depth,
            "product_smiles": self.product_smiles,
            "precursor_smiles": list(self.precursor_smiles),
            "step_cost": self.step_cost,
            "candidate": self.candidate.to_dict(),
        }


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
    frontier_states: int
    solved_routes_found: int
    partial_routes_found: int
    stopped_by_expansion_limit: bool

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
    schema_version: str = MULTISTEP_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible result."""

        return {
            "target_smiles": self.target_smiles,
            "routes": [route.to_dict() for route in self.routes],
            "partial_routes": [
                route.to_dict() for route in self.partial_routes
            ],
            "diagnostics": self.diagnostics.to_dict(),
            "max_depth": self.max_depth,
            "molecular_weight_threshold": self.molecular_weight_threshold,
            "schema_version": self.schema_version,
            "warnings": [
                "Literature matches mean observed in the configured corpus, "
                "not verified commercial availability.",
                "Molecular-weight terminals are a search heuristic and still "
                "require chemist review.",
            ],
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


def _assess_starting_material(
    smiles: str,
    *,
    depth: int,
    literature_index: LiteratureLookup,
    molecular_weight_threshold: float,
    allow_terminal: bool = True,
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
            unresolved_reason=unresolved_reason or "invalid_molecule",
        )
    literature_match = literature_index.lookup(identity)
    reasons = []
    if identity.molecular_weight <= molecular_weight_threshold:
        reasons.append("molecular_weight_threshold")
    if literature_match is not None:
        reasons.append("literature_catalog_match")
    return StartingMaterialAssessment(
        smiles=smiles,
        canonical_smiles=identity.canonical_smiles,
        depth=depth,
        molecular_weight=identity.molecular_weight,
        terminal=allow_terminal and bool(reasons),
        terminal_reasons=tuple(reasons) if allow_terminal else (),
        unresolved_reason=(
            None if allow_terminal and reasons else unresolved_reason
        ),
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


def _route_from_state(
    target_smiles: str,
    state: _RouteState,
    *,
    max_depth: int,
) -> MultistepRetrosynthesisRoute:
    solved = all(leaf.assessment.terminal for leaf in state.leaves)
    step_tokens = tuple(
        step.candidate.proposed_reaction_smiles for step in state.steps
    )
    route_id = digest("ROUTE1", target_smiles, *step_tokens)
    warnings = []
    if not solved:
        warnings.append("ROUTE_INCOMPLETE")
    if any(
        leaf.assessment.unresolved_reason == "maximum_depth"
        or (
            not leaf.assessment.terminal
            and leaf.assessment.depth >= max_depth
        )
        for leaf in state.leaves
    ):
        warnings.append("MAXIMUM_DEPTH_REACHED")
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
        leaves=tuple(
            sorted(
                (leaf.assessment for leaf in state.leaves),
                key=lambda value: (
                    value.canonical_smiles,
                    value.depth,
                    value.unresolved_reason or "",
                ),
            )
        ),
        warnings=tuple(warnings),
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

    def default_expander(
        product_smiles: str,
        top_k: int,
    ) -> tuple[GenericDisconnectionCandidate, ...]:
        return disconnect_operator_ladder(
            product_smiles,
            library,
            top_k=top_k,
            max_templates_to_apply=max_templates_to_apply,
            max_candidates_to_validate=max_candidates_to_validate,
            use_context=use_context,
            include_l0=include_l0,
            diversify=diversify,
        )

    active_expander = expander or default_expander
    root = _Leaf(
        assessment=_assess_starting_material(
            target.canonical_smiles,
            depth=0,
            literature_index=literature_index,
            molecular_weight_threshold=molecular_weight_threshold,
            allow_terminal=False,
        ),
        ancestors=(),
    )
    initial = _RouteState(steps=(), leaves=(root,), cost=0.0)
    queue: list[tuple[tuple[Any, ...], int, _RouteState]] = []
    serial = 0
    heapq.heappush(queue, (_state_priority(initial, max_depth), serial, initial))
    best_cost_by_state = {_state_signature(initial): 0.0}
    expansion_cache: dict[
        str,
        tuple[GenericDisconnectionCandidate, ...],
    ] = {}
    solved_states: dict[str, _RouteState] = {}
    partial_states: dict[tuple[tuple[Any, ...], ...], _RouteState] = {}
    expanded_states = 0
    one_step_calls = 0
    cache_hits = 0
    generated_candidates = 0
    rejected_cycles = 0
    rejected_invalid = 0
    duplicate_states = 0

    while queue and expanded_states < max_expansions:
        _, _, state = heapq.heappop(queue)
        if all(leaf.assessment.terminal for leaf in state.leaves):
            route = _route_from_state(
                target.canonical_smiles,
                state,
                max_depth=max_depth,
            )
            solved_states.setdefault(route.route_id, state)
            if len(solved_states) >= top_k_routes:
                break
            continue
        expandable = _expandable_leaves(state, max_depth)
        if not expandable:
            partial_states.setdefault(_state_signature(state), state)
            continue
        leaf_index, leaf = _select_leaf(expandable)
        product = leaf.assessment.canonical_smiles
        candidates = expansion_cache.get(product)
        if candidates is None:
            candidates = active_expander(product, per_step_top_k)
            expansion_cache[product] = candidates
            one_step_calls += 1
        else:
            cache_hits += 1
        expanded_states += 1
        generated_candidates += len(candidates)
        if not candidates:
            stopped = _Leaf(
                assessment=_assess_starting_material(
                    product,
                    depth=leaf.assessment.depth,
                    literature_index=literature_index,
                    molecular_weight_threshold=molecular_weight_threshold,
                    unresolved_reason="no_candidates",
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
            partial_states[_state_signature(partial)] = partial
            continue

        for candidate_rank, candidate in enumerate(candidates, start=1):
            if candidate.forward_validation_status != "verified_signature":
                rejected_invalid += 1
                continue
            canonical_precursors = molecule_identity(candidate.precursor_smiles)
            if canonical_precursors is None:
                rejected_invalid += 1
                continue
            precursor_values = tuple(
                canonical_precursors.canonical_smiles.split(".")
            )
            ancestor_set = set(leaf.ancestors) | {product}
            if any(value in ancestor_set for value in precursor_values):
                rejected_cycles += 1
                continue
            child_depth = leaf.assessment.depth + 1
            child_leaves = []
            for precursor in precursor_values:
                assessment = _assess_starting_material(
                    precursor,
                    depth=child_depth,
                    literature_index=literature_index,
                    molecular_weight_threshold=molecular_weight_threshold,
                    unresolved_reason=(
                        "maximum_depth" if child_depth >= max_depth else None
                    ),
                )
                child_leaves.append(
                    _Leaf(
                        assessment=assessment,
                        ancestors=(*leaf.ancestors, product),
                    )
                )
            step_cost = round(
                1.0 + max(0.0, 1.0 - float(candidate.score)),
                8,
            )
            step = RetrosynthesisRouteStep(
                step_id=digest(
                    "RSTEP1",
                    product,
                    candidate.proposed_reaction_smiles,
                    str(child_depth),
                ),
                depth=child_depth,
                product_smiles=product,
                precursor_smiles=precursor_values,
                step_cost=step_cost,
                candidate=candidate,
            )
            next_leaves = list(state.leaves)
            next_leaves[leaf_index : leaf_index + 1] = child_leaves
            next_state = _RouteState(
                steps=(*state.steps, step),
                leaves=tuple(next_leaves),
                cost=round(
                    state.cost + step_cost + candidate_rank * 0.000001,
                    8,
                ),
            )
            signature = _state_signature(next_state)
            best_cost = best_cost_by_state.get(signature)
            if best_cost is not None and best_cost <= next_state.cost:
                duplicate_states += 1
                continue
            best_cost_by_state[signature] = next_state.cost
            if all(value.assessment.terminal for value in next_state.leaves):
                route = _route_from_state(
                    target.canonical_smiles,
                    next_state,
                    max_depth=max_depth,
                )
                solved_states.setdefault(route.route_id, next_state)
                continue
            if not _expandable_leaves(next_state, max_depth):
                partial_states[signature] = next_state
                continue
            serial += 1
            heapq.heappush(
                queue,
                (_state_priority(next_state, max_depth), serial, next_state),
            )
        if len(queue) > beam_width:
            queue = heapq.nsmallest(beam_width, queue)
            heapq.heapify(queue)
        if len(solved_states) >= top_k_routes:
            break

    stopped_by_limit = bool(queue and expanded_states >= max_expansions)
    if stopped_by_limit:
        for _, _, state in heapq.nsmallest(top_k_routes, queue):
            leaves = tuple(
                _Leaf(
                    assessment=(
                        leaf.assessment
                        if leaf.assessment.terminal
                        else _assess_starting_material(
                            leaf.assessment.canonical_smiles,
                            depth=leaf.assessment.depth,
                            literature_index=literature_index,
                            molecular_weight_threshold=(
                                molecular_weight_threshold
                            ),
                            unresolved_reason="search_limit",
                        )
                    ),
                    ancestors=leaf.ancestors,
                )
                for leaf in state.leaves
            )
            limited = _RouteState(state.steps, leaves, state.cost)
            partial_states.setdefault(_state_signature(limited), limited)

    solved_routes = tuple(
        sorted(
            (
                _route_from_state(
                    target.canonical_smiles,
                    state,
                    max_depth=max_depth,
                )
                for state in solved_states.values()
            ),
            key=lambda route: (
                route.route_cost,
                route.reaction_count,
                route.route_id,
            ),
        )[:top_k_routes]
    )
    partial_routes = tuple(
        sorted(
            (
                _route_from_state(
                    target.canonical_smiles,
                    state,
                    max_depth=max_depth,
                )
                for state in partial_states.values()
            ),
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
        frontier_states=len(queue),
        solved_routes_found=len(solved_states),
        partial_routes_found=len(partial_states),
        stopped_by_expansion_limit=stopped_by_limit,
    )
    return MultistepRetrosynthesisResult(
        target_smiles=target.canonical_smiles,
        routes=solved_routes,
        partial_routes=partial_routes,
        diagnostics=diagnostics,
        max_depth=max_depth,
        molecular_weight_threshold=molecular_weight_threshold,
    )


__all__ = [
    "MULTISTEP_SCHEMA_VERSION",
    "LiteratureLookup",
    "MultistepRetrosynthesisResult",
    "MultistepRetrosynthesisRoute",
    "MultistepSearchDiagnostics",
    "RetrosynthesisRouteStep",
    "StartingMaterialAssessment",
    "plan_multistep_routes",
]
