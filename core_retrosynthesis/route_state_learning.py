"""Mine literature route-state evidence and safely guide validated route search.

The catalogue in this module is learned only from molecular observations in
multistep route cores.  Reaction names are not routing keys.  The planner hooks
cannot generate a reaction or bypass validation: they only retain an observed,
already validated latent-state action and order otherwise valid search states.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import gzip
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from reactive_taxonomy import identify_reaction_patterns

from .chemistry import digest
from .coupled_route_strategy import extract_coupled_route_strategies
from .generic_library import load_generic_library
from .generic_models import GenericDisconnectionCandidate, GenericTemplateLibrary
from .html_report import reaction_svg
from .latent_state_search import LatentStateActionSelector
from .multistep import MultistepGuidanceState
from .route_core import RouteCoreProjection, RouteCoreStep
from .route_core_conversion import iter_route_core_projections


ROUTE_STATE_LEARNING_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "route_state_learning.v1.json"
)
ROUTE_STATE_LEARNING_VERSION = "route_state_learning.v1"


@dataclass(frozen=True)
class RouteStateLearningPolicy:
    """Validated declarative limits for route-state evidence."""

    definition_id: str
    schema_version: str
    state_release_pattern_ids: tuple[str, ...]
    strict_dependency_classes: tuple[str, ...]
    review_dependency_classes: tuple[str, ...]
    minimum_operator_patent_support: int
    minimum_transition_patent_support: int
    sample_size_per_split: int
    ranking_influence: str


@dataclass(frozen=True)
class RouteStateActionObservation:
    """One graph-observed route action used as auditable evidence."""

    observation_id: str
    route_core_id: str
    patent_id: str | None
    split: str
    source_reaction_id: str | None
    retrosynthetic_depth: int
    action_class: str
    pattern_ids: tuple[str, ...]
    operator_ids: tuple[str, ...]
    exact_core_key: str | None
    typed_core_key: str | None
    shape_core_key: str | None
    reaction_smiles: str
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(
        cls, value: Mapping[str, Any]
    ) -> "RouteStateActionObservation":
        return cls(
            **{
                key: item
                for key, item in value.items()
                if key not in {"pattern_ids", "operator_ids", "warnings"}
            },
            pattern_ids=tuple(value.get("pattern_ids") or ()),
            operator_ids=tuple(value.get("operator_ids") or ()),
            warnings=tuple(value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class RouteOrderingObservation:
    """One lineage-supported precedence relation between physical steps."""

    observation_id: str
    route_core_id: str
    patent_id: str | None
    split: str
    admission_class: str
    dependency_class: str
    forward_first_action_class: str
    forward_second_action_class: str
    reverse_first_operator_ids: tuple[str, ...]
    reverse_second_operator_ids: tuple[str, ...]
    first_reaction_smiles: str
    second_reaction_smiles: str
    evidence: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(
        cls, value: Mapping[str, Any]
    ) -> "RouteOrderingObservation":
        return cls(
            **{
                key: item
                for key, item in value.items()
                if key
                not in {
                    "reverse_first_operator_ids",
                    "reverse_second_operator_ids",
                    "evidence",
                }
            },
            reverse_first_operator_ids=tuple(
                value.get("reverse_first_operator_ids") or ()
            ),
            reverse_second_operator_ids=tuple(
                value.get("reverse_second_operator_ids") or ()
            ),
            evidence=tuple(value.get("evidence") or ()),
        )


@dataclass(frozen=True)
class RouteStateLearningCatalog:
    """Train-derived route guidance plus held-out coverage diagnostics."""

    definition_id: str
    schema_version: str
    source_path: str
    source_sha256: str
    operator_library_path: str
    route_counts: dict[str, int]
    action_counts: dict[str, int]
    dependency_counts: dict[str, int]
    train_state_operator_support: dict[str, int]
    train_state_operator_patent_support: dict[str, int]
    train_reverse_operator_transition_support: dict[str, int]
    train_reverse_operator_transition_patent_support: dict[str, int]
    train_action_transition_support: dict[str, int]
    heldout_metrics: dict[str, dict[str, float | int]]
    action_samples: tuple[RouteStateActionObservation, ...]
    ordering_samples: tuple[RouteOrderingObservation, ...]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RouteStateLearningCatalog":
        if value.get("definition_id") != ROUTE_STATE_LEARNING_VERSION:
            raise ValueError("unexpected route-state learning definition ID")
        if value.get("schema_version") != "1.0":
            raise ValueError("unsupported route-state learning schema")
        return cls(
            **{
                key: item
                for key, item in value.items()
                if key not in {"action_samples", "ordering_samples", "warnings"}
            },
            action_samples=tuple(
                RouteStateActionObservation.from_dict(item)
                for item in value.get("action_samples") or ()
            ),
            ordering_samples=tuple(
                RouteOrderingObservation.from_dict(item)
                for item in value.get("ordering_samples") or ()
            ),
            warnings=tuple(value.get("warnings") or ()),
        )

    @property
    def supported_state_operator_ids(self) -> frozenset[str]:
        policy = load_route_state_learning_policy()
        return frozenset(
            operator_id
            for operator_id, support in self.train_state_operator_patent_support.items()
            if support >= policy.minimum_operator_patent_support
        )


def load_route_state_learning_policy(
    source: str | Path = ROUTE_STATE_LEARNING_DEFINITION_PATH,
) -> RouteStateLearningPolicy:
    """Load and validate the versioned route-state learning policy."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if value.get("definition_id") != ROUTE_STATE_LEARNING_VERSION:
        raise ValueError("unexpected route-state learning definition ID")
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported route-state learning policy schema")
    for key in (
        "minimum_operator_patent_support",
        "minimum_transition_patent_support",
        "sample_size_per_split",
    ):
        item = value.get(key)
        if isinstance(item, bool) or not isinstance(item, int) or item < 1:
            raise ValueError(f"{key} must be a positive integer")
    if value.get("ranking_influence") != (
        "opt_in_validated_state_ordering_and_portfolio_reservation_only"
    ):
        raise ValueError("route-state evidence cannot be ranking-authoritative")
    return RouteStateLearningPolicy(
        definition_id=value["definition_id"],
        schema_version=value["schema_version"],
        state_release_pattern_ids=tuple(value["state_release_pattern_ids"]),
        strict_dependency_classes=tuple(value["strict_dependency_classes"]),
        review_dependency_classes=tuple(value["review_dependency_classes"]),
        minimum_operator_patent_support=value["minimum_operator_patent_support"],
        minimum_transition_patent_support=value[
            "minimum_transition_patent_support"
        ],
        sample_size_per_split=value["sample_size_per_split"],
        ranking_influence=value["ranking_influence"],
    )


def _source_key(projection: RouteCoreProjection, step: RouteCoreStep) -> str:
    if projection.patent_id and step.source_reaction_id:
        return f"{projection.patent_id}:{step.source_reaction_id}"
    return ""


def _operator_indexes(
    library: GenericTemplateLibrary,
) -> tuple[dict[str, tuple[str, ...]], dict[tuple[str, ...], str]]:
    values: dict[str, set[str]] = defaultdict(set)
    for template in library.templates:
        for precedent in template.precedents:
            values[precedent.reaction_id].add(template.operator_id)
    by_edits: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for operator in library.operators:
        by_edits[tuple(operator.edit_tokens)].append(operator.operator_id)
    unique_by_edits = {
        key: items[0] for key, items in by_edits.items() if len(items) == 1
    }
    return (
        {key: tuple(sorted(items)) for key, items in values.items()},
        unique_by_edits,
    )


def _step_action_class(
    step: RouteCoreStep,
    pattern_ids: Sequence[str],
    policy: RouteStateLearningPolicy,
) -> str:
    if set(pattern_ids).intersection(policy.state_release_pattern_ids):
        return "protection_state_change"
    prefixes = {token.split(":", 1)[0] for token in step.edit_tokens}
    if "formed" in prefixes and "broken" in prefixes:
        return "bond_exchange"
    if "formed" in prefixes:
        return "constructive_bond_formation"
    if "broken" in prefixes:
        return "bond_cleavage_or_state_release"
    if "order_changed" in prefixes or "hydrogen_change" in prefixes:
        return "functional_state_change"
    return "unresolved"


def _observe_action(
    projection: RouteCoreProjection,
    step: RouteCoreStep,
    operator_index: Mapping[str, tuple[str, ...]],
    unique_edit_operator_index: Mapping[tuple[str, ...], str],
    policy: RouteStateLearningPolicy,
) -> RouteStateActionObservation:
    # Protection/deprotection patterns are state-release cleavages.  Avoid the
    # comparatively expensive full pattern interpreter for constructive and
    # exchange steps that cannot satisfy that structural prerequisite.
    edit_prefixes = {token.split(":", 1)[0] for token in step.edit_tokens}
    interpretation = (
        identify_reaction_patterns(step.reaction_smiles)
        if "broken" in edit_prefixes and "formed" not in edit_prefixes
        else None
    )
    pattern_ids = tuple(
        sorted(
            {
                match.pattern_id
                for match in (
                    interpretation.pattern_matches if interpretation else ()
                )
            }
        )
    )
    operators = operator_index.get(_source_key(projection, step), ())
    mapping_from_edits = False
    if not operators:
        unique_operator = unique_edit_operator_index.get(tuple(step.edit_tokens))
        if unique_operator:
            operators = (unique_operator,)
            mapping_from_edits = True
    action_class = _step_action_class(step, pattern_ids, policy)
    warnings: list[str] = []
    if mapping_from_edits:
        warnings.append("operator_mapping_from_unique_edit_signature")
    if action_class == "protection_state_change" and not operators:
        warnings.append("no_guarded_operator_mapping")
    return RouteStateActionObservation(
        observation_id=digest(
            "RSLA1", projection.route_core_id, step.reaction_node_id
        ),
        route_core_id=projection.route_core_id,
        patent_id=projection.patent_id,
        split=str(projection.split or "unspecified").lower(),
        source_reaction_id=step.source_reaction_id,
        retrosynthetic_depth=step.retrosynthetic_depth,
        action_class=action_class,
        pattern_ids=pattern_ids,
        operator_ids=operators,
        exact_core_key=step.exact_core_key,
        typed_core_key=step.typed_core_key,
        shape_core_key=step.shape_core_key,
        reaction_smiles=step.reaction_smiles,
        warnings=tuple(warnings),
    )


def _stable_sample(
    values: Iterable[Any], size: int, *, namespace: str
) -> tuple[Any, ...]:
    return tuple(
        sorted(
            values,
            key=lambda item: digest(
                namespace, getattr(item, "observation_id", repr(item))
            ),
        )[:size]
    )


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def mine_route_state_learning_catalog(
    source: str | Path,
    operator_library: str | Path | GenericTemplateLibrary,
    *,
    max_routes: int | None = None,
    policy: RouteStateLearningPolicy | None = None,
) -> RouteStateLearningCatalog:
    """Mine train-only state/operator guidance and score untouched splits."""

    active_policy = policy or load_route_state_learning_policy()
    source_path = Path(source)
    if isinstance(operator_library, GenericTemplateLibrary):
        library = operator_library
        library_path = "<in-memory>"
    else:
        library_path = str(Path(operator_library))
        library = load_generic_library(operator_library)
    operator_index, unique_edit_operator_index = _operator_indexes(library)

    route_counts: Counter[str] = Counter()
    action_counts: Counter[str] = Counter()
    dependency_counts: Counter[str] = Counter()
    actions_by_split: dict[str, list[RouteStateActionObservation]] = defaultdict(list)
    ordering_by_split: dict[str, list[RouteOrderingObservation]] = defaultdict(list)
    state_operator_occurrences: Counter[str] = Counter()
    state_operator_patents: dict[str, set[str]] = defaultdict(set)
    reverse_transition_occurrences: Counter[str] = Counter()
    reverse_transition_patents: dict[str, set[str]] = defaultdict(set)
    action_transition_occurrences: Counter[str] = Counter()

    for route_index, projection in enumerate(iter_route_core_projections(source_path)):
        if max_routes is not None and route_index >= max_routes:
            break
        split = str(projection.split or "unspecified").lower()
        route_counts[split] += 1
        observations = {
            step.reaction_node_id: _observe_action(
                projection,
                step,
                operator_index,
                unique_edit_operator_index,
                active_policy,
            )
            for step in projection.steps
        }
        actions_by_split[split].extend(observations.values())
        for observation in observations.values():
            action_counts[f"{split}:{observation.action_class}"] += 1
            if split == "train" and observation.action_class == (
                "protection_state_change"
            ):
                for operator_id in observation.operator_ids:
                    state_operator_occurrences[operator_id] += 1
                    if observation.patent_id:
                        state_operator_patents[operator_id].add(
                            observation.patent_id
                        )

        for occurrence in extract_coupled_route_strategies(projection):
            if occurrence.admission_class not in {"strict", "review"}:
                continue
            first = observations.get(occurrence.first_reaction_node_id)
            second = observations.get(occurrence.second_reaction_node_id)
            if first is None or second is None:
                continue
            if (
                occurrence.admission_class == "strict"
                and occurrence.dependency_class
                not in active_policy.strict_dependency_classes
            ):
                continue
            if (
                occurrence.admission_class == "review"
                and occurrence.dependency_class
                not in active_policy.review_dependency_classes
            ):
                continue
            dependency_counts[
                f"{split}:{occurrence.dependency_class}"
            ] += 1
            evidence = tuple(occurrence.evidence.rationale)
            ordering = RouteOrderingObservation(
                observation_id=digest("RSLO1", occurrence.occurrence_id),
                route_core_id=projection.route_core_id,
                patent_id=projection.patent_id,
                split=split,
                admission_class=occurrence.admission_class,
                dependency_class=occurrence.dependency_class,
                forward_first_action_class=first.action_class,
                forward_second_action_class=second.action_class,
                reverse_first_operator_ids=second.operator_ids,
                reverse_second_operator_ids=first.operator_ids,
                first_reaction_smiles=first.reaction_smiles,
                second_reaction_smiles=second.reaction_smiles,
                evidence=evidence,
            )
            ordering_by_split[split].append(ordering)
            if split != "train":
                continue
            action_key = f"{second.action_class}>{first.action_class}"
            action_transition_occurrences[action_key] += 1
            for first_operator in second.operator_ids:
                for second_operator in first.operator_ids:
                    key = f"{first_operator}>{second_operator}"
                    reverse_transition_occurrences[key] += 1
                    if projection.patent_id:
                        reverse_transition_patents[key].add(projection.patent_id)

    supported_state_operators = {
        operator_id
        for operator_id, patents in state_operator_patents.items()
        if len(patents) >= active_policy.minimum_operator_patent_support
    }
    supported_transitions = {
        key
        for key, patents in reverse_transition_patents.items()
        if len(patents) >= active_policy.minimum_transition_patent_support
    }
    supported_action_transitions = set(action_transition_occurrences)
    heldout_metrics: dict[str, dict[str, float | int]] = {}
    for split in ("validation", "val", "test"):
        if not actions_by_split.get(split) and not ordering_by_split.get(split):
            continue
        state_actions = [
            item
            for item in actions_by_split[split]
            if item.action_class == "protection_state_change"
        ]
        mapped_state_actions = [item for item in state_actions if item.operator_ids]
        covered_state_actions = [
            item
            for item in mapped_state_actions
            if set(item.operator_ids).intersection(supported_state_operators)
        ]
        mapped_orderings = [
            item
            for item in ordering_by_split[split]
            if item.reverse_first_operator_ids and item.reverse_second_operator_ids
        ]
        covered_orderings = [
            item
            for item in mapped_orderings
            if any(
                f"{first}>{second}" in supported_transitions
                for first in item.reverse_first_operator_ids
                for second in item.reverse_second_operator_ids
            )
        ]
        class_covered = [
            item
            for item in ordering_by_split[split]
            if (
                f"{item.forward_second_action_class}>"
                f"{item.forward_first_action_class}"
            )
            in supported_action_transitions
        ]
        heldout_metrics[split] = {
            "state_action_count": len(state_actions),
            "mapped_state_action_count": len(mapped_state_actions),
            "covered_mapped_state_action_count": len(covered_state_actions),
            "mapped_state_action_coverage": (
                len(covered_state_actions) / len(mapped_state_actions)
                if mapped_state_actions
                else 0.0
            ),
            "ordering_count": len(ordering_by_split[split]),
            "mapped_operator_ordering_count": len(mapped_orderings),
            "covered_operator_ordering_count": len(covered_orderings),
            "mapped_operator_ordering_coverage": (
                len(covered_orderings) / len(mapped_orderings)
                if mapped_orderings
                else 0.0
            ),
            "covered_action_class_ordering_count": len(class_covered),
            "action_class_ordering_coverage": (
                len(class_covered) / len(ordering_by_split[split])
                if ordering_by_split[split]
                else 0.0
            ),
        }

    sampled_actions: list[RouteStateActionObservation] = []
    sampled_ordering: list[RouteOrderingObservation] = []
    for split in sorted(actions_by_split):
        state_values = [
            item
            for item in actions_by_split[split]
            if item.action_class == "protection_state_change"
        ]
        sampled_actions.extend(
            _stable_sample(
                state_values,
                active_policy.sample_size_per_split,
                namespace=f"RSLA-SAMPLE:{split}",
            )
        )
        sampled_ordering.extend(
            _stable_sample(
                ordering_by_split[split],
                active_policy.sample_size_per_split,
                namespace=f"RSLO-SAMPLE:{split}",
            )
        )
    warnings = (
        "Route topology supplies graph precedence, not guaranteed experimental chronology.",
        "Independent branches are intentionally not converted into ordering constraints.",
        "Protection labels are optional interpretations of structural state-release evidence.",
        "Held-out metrics measure catalogue coverage, not synthesis success.",
    )
    return RouteStateLearningCatalog(
        definition_id=active_policy.definition_id,
        schema_version="1.0",
        source_path=str(source_path),
        source_sha256=_sha256(source_path),
        operator_library_path=library_path,
        route_counts=dict(sorted(route_counts.items())),
        action_counts=dict(sorted(action_counts.items())),
        dependency_counts=dict(sorted(dependency_counts.items())),
        train_state_operator_support=dict(sorted(state_operator_occurrences.items())),
        train_state_operator_patent_support={
            key: len(value) for key, value in sorted(state_operator_patents.items())
        },
        train_reverse_operator_transition_support=dict(
            sorted(reverse_transition_occurrences.items())
        ),
        train_reverse_operator_transition_patent_support={
            key: len(value)
            for key, value in sorted(reverse_transition_patents.items())
        },
        train_action_transition_support=dict(
            sorted(action_transition_occurrences.items())
        ),
        heldout_metrics=heldout_metrics,
        action_samples=tuple(sampled_actions),
        ordering_samples=tuple(sampled_ordering),
        warnings=warnings,
    )


def save_route_state_learning_catalog(
    catalog: RouteStateLearningCatalog, destination: str | Path
) -> None:
    """Write a deterministic JSON or JSON.GZ catalogue."""

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        catalog.to_dict(), indent=2, sort_keys=True, ensure_ascii=False
    ) + "\n"
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(payload.encode("utf-8"), mtime=0))
    else:
        path.write_text(payload, encoding="utf-8")


def load_route_state_learning_catalog(
    source: str | Path,
) -> RouteStateLearningCatalog:
    """Load a route-state learning catalogue."""

    path = Path(source)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
    else:
        value = json.loads(path.read_text(encoding="utf-8"))
    return RouteStateLearningCatalog.from_dict(value)


class LiteratureRouteActionSelector:
    """Reserve train-supported latent-state actions in a validated portfolio."""

    def __init__(
        self,
        catalog: RouteStateLearningCatalog,
        base_selector: LatentStateActionSelector | None = None,
    ) -> None:
        self.catalog = catalog
        self.base_selector = base_selector or LatentStateActionSelector()
        self.definition_id = catalog.definition_id
        self.candidate_pool_multiplier = max(
            2, self.base_selector.candidate_pool_multiplier
        )
        self.continuation_definition_id = (
            self.base_selector.continuation_definition_id
        )
        self.minimum_expansions_per_first_action = (
            self.base_selector.minimum_expansions_per_first_action
        )
        self.maximum_active_first_actions = (
            self.base_selector.maximum_active_first_actions
        )
        self._supported = catalog.supported_state_operator_ids

    def __call__(
        self,
        product_smiles: str,
        candidates: tuple[GenericDisconnectionCandidate, ...],
        top_k: int,
    ) -> tuple[GenericDisconnectionCandidate, ...]:
        selected = list(self.base_selector(product_smiles, candidates, top_k))
        supported = next(
            (
                item
                for item in candidates
                if item.operator_id in self._supported
            ),
            None,
        )
        if supported is None or supported in selected or top_k < 2:
            return tuple(selected)
        if len(selected) >= top_k:
            selected[-1] = supported
        else:
            selected.append(supported)
        positions = {id(item): index for index, item in enumerate(candidates)}
        return tuple(sorted(selected, key=lambda item: positions[id(item)]))

    def state_diversity_key(self, state: Any) -> str:
        return self.base_selector.state_diversity_key(state)

    def continuation_lane_key(self, state: Any) -> str:
        return self.base_selector.continuation_lane_key(state)

    def continuation_lane_class(self, state: Any) -> str:
        return self.base_selector.continuation_lane_class(state)

    def select_routes(self, routes: tuple[Any, ...], limit: int) -> tuple[Any, ...]:
        return self.base_selector.select_routes(routes, limit)


class LiteratureRouteOrderingGuidance:
    """Order valid search states with train-supported operator sequences."""

    definition_id = ROUTE_STATE_LEARNING_VERSION

    def __init__(self, catalog: RouteStateLearningCatalog) -> None:
        policy = load_route_state_learning_policy()
        self._support = {
            key: value
            for key, value in catalog.train_reverse_operator_transition_support.items()
            if catalog.train_reverse_operator_transition_patent_support.get(key, 0)
            >= policy.minimum_transition_patent_support
        }

    def state_priority(self, state: MultistepGuidanceState) -> tuple[Any, ...]:
        if len(state.steps) < 2:
            return (0, 0.0)
        first = state.steps[-2].candidate.operator_id
        second = state.steps[-1].candidate.operator_id
        support = self._support.get(f"{first}>{second}", 0)
        # This optional priority is applied only after actions pass chemistry gates.
        return (0 if support else 1, -math.log1p(support))

    def select_leaf(
        self,
        state: MultistepGuidanceState,
        expandable_indices: tuple[int, ...],
    ) -> int:
        return min(
            expandable_indices,
            key=lambda index: (
                -(state.leaves[index].molecular_weight or 0.0),
                state.leaves[index].canonical_smiles,
                index,
            ),
        )


def _percent(value: float | int) -> str:
    return f"{100.0 * float(value):.1f}%"


def _balanced_review_sample(
    values: Sequence[Any], *, per_split: int = 6
) -> tuple[Any, ...]:
    selected: list[Any] = []
    counts: Counter[str] = Counter()
    for item in values:
        split = str(getattr(item, "split", "unspecified"))
        if counts[split] >= per_split:
            continue
        selected.append(item)
        counts[split] += 1
    return tuple(selected)


def write_route_state_learning_html(
    catalog: RouteStateLearningCatalog,
    destination: str | Path,
    *,
    title: str = "Literature route-state learning review",
) -> None:
    """Render a self-contained graphical audit of learned route evidence."""

    metrics = []
    for split, values in sorted(catalog.heldout_metrics.items()):
        metrics.append(
            "<tr>"
            f"<td>{html.escape(split)}</td>"
            f"<td>{values['state_action_count']}</td>"
            f"<td>{values['mapped_state_action_count']}</td>"
            f"<td>{_percent(values['mapped_state_action_coverage'])}</td>"
            f"<td>{values['ordering_count']}</td>"
            f"<td>{_percent(values['mapped_operator_ordering_coverage'])}</td>"
            f"<td>{_percent(values['action_class_ordering_coverage'])}</td>"
            "</tr>"
        )
    action_cards = []
    for item in _balanced_review_sample(catalog.action_samples):
        action_cards.append(
            "<article class='card'>"
            f"<div class='tag'>{html.escape(item.split)} · latent-state release</div>"
            f"<div class='graphic'>{reaction_svg(item.reaction_smiles, sub_image_size=(360, 220))}</div>"
            f"<p><b>Patterns</b> {html.escape(', '.join(item.pattern_ids) or 'none')}</p>"
            f"<p><b>Guarded operators</b> {len(item.operator_ids)}</p>"
            f"<p class='mono'>{html.escape(item.source_reaction_id or 'unknown source')}</p>"
            "</article>"
        )
    ordering_cards = []
    for item in _balanced_review_sample(catalog.ordering_samples):
        ordering_cards.append(
            "<article class='order-card'>"
            f"<div class='tag'>{html.escape(item.split)} · {html.escape(item.dependency_class)}</div>"
            "<div class='two'>"
            f"<div><h4>Forward step 1</h4>{reaction_svg(item.first_reaction_smiles, sub_image_size=(330, 190))}</div>"
            f"<div><h4>Forward step 2</h4>{reaction_svg(item.second_reaction_smiles, sub_image_size=(330, 190))}</div>"
            "</div>"
            f"<p><b>Planner order</b> {html.escape(item.forward_second_action_class)} → {html.escape(item.forward_first_action_class)}</p>"
            f"<p>{html.escape('; '.join(item.evidence[:3]))}</p>"
            "</article>"
        )
    warnings = "".join(f"<li>{html.escape(item)}</li>" for item in catalog.warnings)
    document = f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'><title>{html.escape(title)}</title>
<style>
body{{font-family:Inter,Segoe UI,sans-serif;margin:0;background:#f4f6f8;color:#18222d}}
main{{max-width:1450px;margin:auto;padding:28px}} h1{{margin-bottom:6px}}
.note,.panel,.card,.order-card{{background:white;border:1px solid #dce2e7;border-radius:12px;box-shadow:0 2px 8px #1b27331a}}
.note,.panel{{padding:18px;margin:16px 0}} .cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:14px}}
.card,.order-card{{padding:14px;overflow:hidden}} .order-card{{margin:14px 0}} .graphic svg{{width:100%;height:225px}}
.two{{display:grid;grid-template-columns:1fr 1fr;gap:14px}} .two svg{{width:100%;height:200px}}
.tag{{display:inline-block;padding:5px 9px;border-radius:999px;background:#e7f1ff;color:#174f86;font-size:12px;font-weight:700}}
table{{border-collapse:collapse;width:100%}} th,td{{padding:9px;border-bottom:1px solid #e5e9ed;text-align:left}} th{{background:#eef2f5}}
.mono{{font-family:ui-monospace,Consolas,monospace;font-size:12px;color:#596675}}
@media(max-width:760px){{.two{{grid-template-columns:1fr}}}}
</style></head><body><main>
<h1>{html.escape(title)}</h1>
<p>Structural evidence learned from training patents; validation/test are used only for coverage measurement.</p>
<section class='panel'><h2>Held-out coverage (not synthesis accuracy)</h2><table><thead><tr><th>Split</th><th>State actions</th><th>Mapped</th><th>State operator coverage</th><th>Dependencies</th><th>Operator-order coverage</th><th>Class-order coverage</th></tr></thead><tbody>{''.join(metrics)}</tbody></table></section>
<section class='note'><h2>Safety boundary</h2><ul>{warnings}</ul><p>The planner integration can only reserve or reorder actions that passed the existing graph, forward-validation, compatibility, and realism gates.</p></section>
<h2>Latent/protected state examples</h2><section class='cards'>{''.join(action_cards) or '<p>No state-action samples.</p>'}</section>
<h2>Lineage-supported ordering examples</h2><section>{''.join(ordering_cards) or '<p>No ordering samples.</p>'}</section>
</main></body></html>"""
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")


__all__ = [
    "LiteratureRouteActionSelector",
    "LiteratureRouteOrderingGuidance",
    "ROUTE_STATE_LEARNING_VERSION",
    "RouteOrderingObservation",
    "RouteStateActionObservation",
    "RouteStateLearningCatalog",
    "RouteStateLearningPolicy",
    "load_route_state_learning_catalog",
    "load_route_state_learning_policy",
    "mine_route_state_learning_catalog",
    "save_route_state_learning_catalog",
    "write_route_state_learning_html",
]
