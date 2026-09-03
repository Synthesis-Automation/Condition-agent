"""Review-only strategic diversity for validated retrosynthesis actions.

This module does not generate reactions. It classifies already forward-verified
actions by how target atoms are distributed across their precursors, then
retains a bounded portfolio of complementary route-entry modes.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Protocol, Sequence, TypeVar

from rdkit import Chem

from .chemistry import digest
from .generic_models import GenericDisconnectionCandidate
from .portfolio_continuation import (
    PortfolioContinuationPolicy,
    load_portfolio_continuation_policy,
)


LATENT_STATE_ROUTE_SEARCH_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "latent_state_route_search.v1.json"
)
LATENT_STATE_ROUTE_SEARCH_VERSION = "latent_state_route_search.v1"
ROUTE_ACTION_CLASSES = (
    "convergent_assembly",
    "single_carrier_installation",
    "unary_state_change",
    "ring_reorganization",
    "multicomponent_reorganization",
    "unclassified",
)


@dataclass(frozen=True)
class LatentStateRouteSearchPolicy:
    """Validated selection limits for the strategic route portfolio."""

    definition_id: str
    schema_version: str
    candidate_pool_multiplier: int
    single_carrier_target_atom_fraction: float
    action_class_order: tuple[str, ...]
    candidate_reserve_per_class: tuple[tuple[str, int], ...]
    maximum_routes_per_first_action_class: int
    preserve_best_overall_route: bool
    ranking_influence: str

    def reserve(self, action_class: str) -> int:
        return dict(self.candidate_reserve_per_class).get(action_class, 0)


@dataclass(frozen=True)
class RouteActionClassification:
    """Graph-derived strategic shape of one validated reverse action."""

    action_class: str
    target_heavy_atom_count: int
    target_atom_coverage: float
    contributing_precursor_count: int
    precursor_component_count: int
    largest_target_atom_fraction: float
    evidence: tuple[str, ...]
    definition_id: str = LATENT_STATE_ROUTE_SEARCH_VERSION

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class LatentStateActionSelector:
    """Planner hook retaining distinct validated action classes per product."""

    def __init__(
        self,
        policy: LatentStateRouteSearchPolicy | None = None,
        continuation_policy: PortfolioContinuationPolicy | None = None,
    ) -> None:
        self.policy = policy or load_latent_state_route_search_policy()
        self.continuation_policy = (
            continuation_policy or load_portfolio_continuation_policy()
        )
        self.definition_id = self.policy.definition_id
        self.candidate_pool_multiplier = self.policy.candidate_pool_multiplier
        self.continuation_definition_id = (
            self.continuation_policy.definition_id
        )
        self.minimum_expansions_per_first_action = (
            self.continuation_policy.minimum_expansions_per_first_action
        )
        self.maximum_active_first_actions = (
            self.continuation_policy.maximum_active_first_actions
        )

    def __call__(
        self,
        product_smiles: str,
        candidates: tuple[GenericDisconnectionCandidate, ...],
        top_k: int,
    ) -> tuple[GenericDisconnectionCandidate, ...]:
        del product_smiles
        return select_latent_state_action_portfolio(
            candidates,
            top_k=top_k,
            policy=self.policy,
        )

    def state_diversity_key(self, state: Any) -> str:
        """Return the distinct target-forming action lane for a state."""

        return first_route_action_lane_id(state)

    def continuation_lane_key(self, state: Any) -> str:
        """Return the first-action identity used for fair continuation."""

        return first_route_action_lane_id(state)

    def continuation_lane_class(self, state: Any) -> str:
        """Return the graph-derived route shape for continuation diagnostics."""

        return first_route_action_class(state, policy=self.policy)

    def select_routes(
        self,
        routes: tuple[Any, ...],
        limit: int,
    ) -> tuple[Any, ...]:
        """Return a bounded final portfolio across first-action classes."""

        return select_route_class_portfolio(
            routes,
            limit=limit,
            policy=self.policy,
        )


class RouteLike(Protocol):
    steps: Sequence[Any]


_RouteValue = TypeVar("_RouteValue", bound=RouteLike)


def validate_latent_state_route_search_policy(value: Mapping[str, Any]) -> None:
    """Reject incomplete or ranking-active latent-state policies."""

    if value.get("definition_id") != LATENT_STATE_ROUTE_SEARCH_VERSION:
        raise ValueError("unexpected latent-state route-search definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported latent-state route-search schema")
    multiplier = value.get("candidate_pool_multiplier")
    if isinstance(multiplier, bool) or not isinstance(multiplier, int) or multiplier < 1:
        raise ValueError("candidate pool multiplier must be a positive integer")
    fraction = value.get("single_carrier_target_atom_fraction")
    if not isinstance(fraction, (int, float)) or isinstance(fraction, bool):
        raise ValueError("single-carrier threshold must be numeric")
    if not 0.5 < float(fraction) <= 1.0:
        raise ValueError("single-carrier threshold must be in (0.5, 1]")
    if tuple(value.get("action_class_order") or ()) != ROUTE_ACTION_CLASSES:
        raise ValueError("latent-state action-class order is unsupported")
    reserves = value.get("candidate_reserve_per_class")
    if not isinstance(reserves, Mapping) or set(reserves) != set(ROUTE_ACTION_CLASSES):
        raise ValueError("candidate reserves must cover every action class")
    if any(
        isinstance(count, bool) or not isinstance(count, int) or count < 0
        for count in reserves.values()
    ):
        raise ValueError("candidate reserves must be nonnegative integers")
    portfolio = value.get("route_portfolio")
    if not isinstance(portfolio, Mapping):
        raise ValueError("route portfolio policy is required")
    maximum = portfolio.get("maximum_routes_per_first_action_class")
    if isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1:
        raise ValueError("route class maximum must be a positive integer")
    if portfolio.get("preserve_best_overall_route") is not True:
        raise ValueError("the best overall route must be preserved")
    if portfolio.get("ranking_influence") != "review_only_opt_in":
        raise ValueError("latent-state portfolio must remain review-only")


@lru_cache(maxsize=1)
def load_latent_state_route_search_policy() -> LatentStateRouteSearchPolicy:
    """Load the canonical strategic route-portfolio policy."""

    value = json.loads(
        LATENT_STATE_ROUTE_SEARCH_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    validate_latent_state_route_search_policy(value)
    portfolio = value["route_portfolio"]
    reserves = value["candidate_reserve_per_class"]
    return LatentStateRouteSearchPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        candidate_pool_multiplier=int(value["candidate_pool_multiplier"]),
        single_carrier_target_atom_fraction=float(
            value["single_carrier_target_atom_fraction"]
        ),
        action_class_order=tuple(value["action_class_order"]),
        candidate_reserve_per_class=tuple(
            (action_class, int(reserves[action_class]))
            for action_class in ROUTE_ACTION_CLASSES
        ),
        maximum_routes_per_first_action_class=int(
            portfolio["maximum_routes_per_first_action_class"]
        ),
        preserve_best_overall_route=True,
        ranking_influence=str(portfolio["ranking_influence"]),
    )


def _heavy_map_sets(reaction_smiles: str) -> tuple[tuple[set[int], ...], set[int]]:
    if reaction_smiles.count(">>") != 1:
        return (), set()
    reactants, product = reaction_smiles.split(">>", 1)
    product_molecule = Chem.MolFromSmiles(product)
    if product_molecule is None:
        return (), set()
    product_maps = {
        int(atom.GetAtomMapNum())
        for atom in product_molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and int(atom.GetAtomMapNum()) > 0
    }
    components = []
    for component in reactants.split("."):
        molecule = Chem.MolFromSmiles(component)
        if molecule is None:
            return (), set()
        components.append(
            {
                int(atom.GetAtomMapNum())
                for atom in molecule.GetAtoms()
                if atom.GetAtomicNum() > 1 and int(atom.GetAtomMapNum()) > 0
            }
        )
    return tuple(components), product_maps


@lru_cache(maxsize=50_000)
def _classify_reaction(
    reaction_smiles: str,
    strategic_class: str,
    threshold: float,
) -> RouteActionClassification:
    precursor_maps, product_maps = _heavy_map_sets(reaction_smiles)
    if not precursor_maps or not product_maps:
        precursor_side = reaction_smiles.split(">>", 1)[0]
        component_count = len(precursor_side.split(".")) if precursor_side else 0
        action_class = (
            "unary_state_change" if component_count == 1 else "unclassified"
        )
        return RouteActionClassification(
            action_class=action_class,
            target_heavy_atom_count=len(product_maps),
            target_atom_coverage=0.0,
            contributing_precursor_count=0,
            precursor_component_count=component_count,
            largest_target_atom_fraction=0.0,
            evidence=("MAPPED_TARGET_ATOM_DISTRIBUTION_UNAVAILABLE",),
        )
    contributions = tuple(len(maps & product_maps) for maps in precursor_maps)
    contributing = tuple(value for value in contributions if value > 0)
    coverage = len(set().union(*precursor_maps) & product_maps) / len(product_maps)
    largest = max(contributing, default=0) / len(product_maps)
    if strategic_class == "ring_construction":
        action_class = "ring_reorganization"
    elif len(precursor_maps) == 1:
        action_class = "unary_state_change"
    elif largest >= threshold:
        action_class = "single_carrier_installation"
    elif len(contributing) >= 2 and coverage >= 0.999999:
        action_class = "convergent_assembly"
    else:
        action_class = "multicomponent_reorganization"
    return RouteActionClassification(
        action_class=action_class,
        target_heavy_atom_count=len(product_maps),
        target_atom_coverage=round(coverage, 8),
        contributing_precursor_count=len(contributing),
        precursor_component_count=len(precursor_maps),
        largest_target_atom_fraction=round(largest, 8),
        evidence=(
            "MAPPED_REACTION_TARGET_ATOMS",
            "TARGET_ATOM_DISTRIBUTION_ACROSS_PRECURSORS",
        ),
    )


def classify_route_action(
    candidate: GenericDisconnectionCandidate,
    *,
    policy: LatentStateRouteSearchPolicy | None = None,
) -> RouteActionClassification:
    """Classify one already validated candidate without reaction-name routing."""

    resolved = policy or load_latent_state_route_search_policy()
    reaction = (
        candidate.condition_query_reaction_smiles
        or candidate.proposed_reaction_smiles
    )
    return _classify_reaction(
        reaction,
        candidate.strategic_class,
        resolved.single_carrier_target_atom_fraction,
    )


def select_latent_state_action_portfolio(
    candidates: Iterable[GenericDisconnectionCandidate],
    *,
    top_k: int,
    policy: LatentStateRouteSearchPolicy | None = None,
) -> tuple[GenericDisconnectionCandidate, ...]:
    """Retain bounded class diversity while preserving within-pool rank order."""

    if top_k < 1:
        raise ValueError("latent-state action portfolio top-k must be positive")
    resolved = policy or load_latent_state_route_search_policy()
    values = tuple(candidates)
    if len(values) <= top_k:
        return values
    ranked_index = {id(candidate): index for index, candidate in enumerate(values)}
    by_class: dict[str, list[GenericDisconnectionCandidate]] = {
        action_class: [] for action_class in resolved.action_class_order
    }
    for candidate in values:
        classification = classify_route_action(candidate, policy=resolved)
        by_class[classification.action_class].append(candidate)
    selected: list[GenericDisconnectionCandidate] = []
    selected_ids: set[int] = set()
    for action_class in resolved.action_class_order:
        for candidate in by_class[action_class][: resolved.reserve(action_class)]:
            if len(selected) >= top_k:
                break
            selected.append(candidate)
            selected_ids.add(id(candidate))
    for candidate in values:
        if len(selected) >= top_k:
            break
        if id(candidate) in selected_ids:
            continue
        selected.append(candidate)
        selected_ids.add(id(candidate))
    return tuple(sorted(selected, key=lambda item: ranked_index[id(item)]))


def first_route_action_class(
    route: RouteLike,
    *,
    policy: LatentStateRouteSearchPolicy | None = None,
) -> str:
    """Return the graph-derived class of the route's target-forming step."""

    if not route.steps:
        return "unclassified"
    first = min(route.steps, key=lambda step: (step.depth, step.step_id))
    return classify_route_action(first.candidate, policy=policy).action_class


def first_route_action_lane_id(route: RouteLike) -> str:
    """Return a stable identity for a route's distinct first reverse action."""

    if not route.steps:
        return "root"
    first = min(route.steps, key=lambda step: (step.depth, step.step_id))
    return digest(
        "ROUTELANE1",
        first.candidate.proposed_reaction_smiles,
        first.candidate.operator_id,
        first.candidate.disconnection_site_key,
    )


def select_route_class_portfolio(
    routes: Iterable[_RouteValue],
    *,
    limit: int,
    policy: LatentStateRouteSearchPolicy | None = None,
) -> tuple[_RouteValue, ...]:
    """Interleave ranked routes across distinct target-forming action classes."""

    if limit < 1:
        raise ValueError("route portfolio limit must be positive")
    resolved = policy or load_latent_state_route_search_policy()
    values = tuple(routes)
    if len(values) <= limit:
        return values
    groups: dict[str, list[_RouteValue]] = {
        action_class: [] for action_class in resolved.action_class_order
    }
    for route in values:
        action_class = first_route_action_class(route, policy=resolved)
        groups[action_class].append(route)
    selected: list[_RouteValue] = []
    if resolved.preserve_best_overall_route and values:
        selected.append(values[0])
    selected_ids = {id(item) for item in selected}
    depth = 0
    while len(selected) < limit:
        added = False
        for action_class in resolved.action_class_order:
            group = groups[action_class]
            if depth >= min(
                len(group), resolved.maximum_routes_per_first_action_class
            ):
                continue
            route = group[depth]
            if id(route) not in selected_ids:
                selected.append(route)
                selected_ids.add(id(route))
                added = True
                if len(selected) >= limit:
                    break
        if not added:
            break
        depth += 1
    for route in values:
        if len(selected) >= limit:
            break
        if id(route) not in selected_ids:
            selected.append(route)
    return tuple(selected)


__all__ = [
    "LATENT_STATE_ROUTE_SEARCH_DEFINITION_PATH",
    "LATENT_STATE_ROUTE_SEARCH_VERSION",
    "ROUTE_ACTION_CLASSES",
    "LatentStateRouteSearchPolicy",
    "LatentStateActionSelector",
    "RouteActionClassification",
    "classify_route_action",
    "first_route_action_class",
    "first_route_action_lane_id",
    "load_latent_state_route_search_policy",
    "select_latent_state_action_portfolio",
    "select_route_class_portfolio",
    "validate_latent_state_route_search_policy",
]
