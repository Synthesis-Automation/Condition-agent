"""Held-out evaluation of promoted v1 two-step retrosynthesis actions.

The v1 structural relationship is the admission rule.  Each promoted logical
action still executes two independently forward-validated generic operators;
v2 dependency labels are retained only as review metadata.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from .chemistry import canonical_smiles, digest
from .coupled_route_strategy import extract_coupled_route_strategies
from .coupled_strategy_search import V1_ADMITTED_RELATIONSHIPS
from .generic_compiler import generic_operator_identity_from_observation
from .generic_models import (
    GenericDisconnectionCandidate,
    GenericSearchDiagnostics,
    GenericTemplateLibrary,
)
from .generic_search import disconnect_generic_target_detailed
from .route_core_conversion import iter_route_core_projections


COUPLED_STRATEGY_EVALUATION_SCHEMA_VERSION = "1.0"
COUPLED_STRATEGY_EVALUATION_ALGORITHM_VERSION = "coupled_strategy_eval.v1"
FROZEN_V1_PANEL_SCHEMA_VERSION = "1.0"
FROZEN_V1_PANEL_ALGORITHM_VERSION = "v1_coupled_strategy_panel.v1"
_PATENT_TOKEN = re.compile(r"[^A-Z0-9]+")
SearchFunction = Callable[
    ..., tuple[tuple[GenericDisconnectionCandidate, ...], GenericSearchDiagnostics]
]


@dataclass(frozen=True)
class CoupledStrategyEvaluationConfig:
    """Fixed search and panel budgets for one paired evaluation."""

    panel_size: int = 12
    minimum_training_patents: int = 2
    top_k: int = 5
    max_templates_to_apply: int = 80
    max_candidates_to_validate: int = 20
    seed: int = 20260818

    def __post_init__(self) -> None:
        if (
            min(
                self.panel_size,
                self.minimum_training_patents,
                self.top_k,
                self.max_templates_to_apply,
                self.max_candidates_to_validate,
            )
            < 1
        ):
            raise ValueError("evaluation budgets must be positive")


@dataclass(frozen=True)
class PromotedV1OperatorPair:
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
        if len(self.training_patent_ids) < 2:
            raise ValueError("promoted pairs require independent patent support")

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["training_patent_ids"] = list(self.training_patent_ids)
        value["v2_dependency_counts"] = dict(self.v2_dependency_counts)
        return value


@dataclass(frozen=True)
class CoupledStrategyEvaluationCase:
    """One patent-held-out observed route segment used as ground truth."""

    case_id: str
    occurrence_id: str
    strategy_id: str
    patent_id: str
    split: str
    relationship_class: str
    v2_dependency_class: str
    target_smiles: str
    expected_intermediate_smiles: str
    expected_terminal_precursor_smiles: str
    observed_first_reaction_smiles: str
    observed_second_reaction_smiles: str
    exact_target_seen_in_training: bool
    target_scaffold_seen_in_training: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FrozenV1HeldoutPanel:
    """Library-independent patent-held-out operator-pair evaluation panel."""

    panel_id: str
    route_core_source: str
    route_core_sha256: str
    config: CoupledStrategyEvaluationConfig
    strategies: tuple[PromotedV1OperatorPair, ...]
    cases: tuple[CoupledStrategyEvaluationCase, ...]
    required_strategy_ids: tuple[str, ...] = ()
    schema_version: str = FROZEN_V1_PANEL_SCHEMA_VERSION
    algorithm_version: str = FROZEN_V1_PANEL_ALGORITHM_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != FROZEN_V1_PANEL_SCHEMA_VERSION:
            raise ValueError("unsupported frozen v1 panel schema")
        if self.algorithm_version != FROZEN_V1_PANEL_ALGORITHM_VERSION:
            raise ValueError("unsupported frozen v1 panel algorithm")
        if not self.panel_id.startswith("CRV1PANEL2:"):
            raise ValueError("frozen v1 panel requires a content identity")
        strategy_ids = {item.strategy_id for item in self.strategies}
        if not self.cases or len(self.cases) > self.config.panel_size:
            raise ValueError("frozen v1 panel case count is invalid")
        if any(item.strategy_id not in strategy_ids for item in self.cases):
            raise ValueError("frozen v1 panel case lacks its strategy")
        if not set(self.required_strategy_ids).issubset(strategy_ids):
            raise ValueError("required frozen-panel strategy was not selected")

    def to_dict(self) -> dict[str, Any]:
        return {
            "artifact_type": "v1_coupled_strategy_frozen_panel",
            "schema_version": self.schema_version,
            "algorithm_version": self.algorithm_version,
            "panel_id": self.panel_id,
            "route_core_source": self.route_core_source,
            "route_core_sha256": self.route_core_sha256,
            "config": asdict(self.config),
            "required_strategy_ids": list(self.required_strategy_ids),
            "strategies": [item.to_dict() for item in self.strategies],
            "cases": [item.to_dict() for item in self.cases],
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "FrozenV1HeldoutPanel":
        if value.get("artifact_type") != "v1_coupled_strategy_frozen_panel":
            raise ValueError("unexpected frozen v1 panel artifact")
        strategies = tuple(
            PromotedV1OperatorPair(
                **{
                    **raw,
                    "training_patent_ids": tuple(raw.get("training_patent_ids") or ()),
                    "v2_dependency_counts": tuple(
                        sorted((raw.get("v2_dependency_counts") or {}).items())
                    ),
                }
            )
            for raw in value.get("strategies") or ()
        )
        cases = tuple(
            CoupledStrategyEvaluationCase(**raw) for raw in value.get("cases") or ()
        )
        return cls(
            panel_id=str(value.get("panel_id") or ""),
            route_core_source=str(value.get("route_core_source") or ""),
            route_core_sha256=str(value.get("route_core_sha256") or ""),
            config=CoupledStrategyEvaluationConfig(**(value.get("config") or {})),
            strategies=strategies,
            cases=cases,
            required_strategy_ids=tuple(
                str(item) for item in value.get("required_strategy_ids") or ()
            ),
            schema_version=str(value.get("schema_version") or ""),
            algorithm_version=str(value.get("algorithm_version") or ""),
        )


@dataclass(frozen=True)
class PromotedTwoStepAction:
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

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PromotedTwoStepQueryAction:
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

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible query action."""

        value = asdict(self)
        value["v2_dependency_counts"] = dict(self.v2_dependency_counts)
        return value


@dataclass(frozen=True)
class CoupledStrategyQueryDiagnostics:
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

    def to_dict(self) -> dict[str, int]:
        """Return JSON-compatible diagnostics."""

        return asdict(self)


@dataclass(frozen=True)
class CoupledStrategyQueryResult:
    """Promoted two-step actions with ordinary one-step fallbacks retained."""

    target_smiles: str
    actions: tuple[PromotedTwoStepQueryAction, ...]
    one_step_fallbacks: tuple[dict[str, Any], ...]
    diagnostics: CoupledStrategyQueryDiagnostics
    warnings: tuple[str, ...]
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
        }


@dataclass(frozen=True)
class CoupledStrategyCaseResult:
    """Paired ordinary-depth-two and promoted-v1 result for one target."""

    case: CoupledStrategyEvaluationCase
    baseline_intermediate_rank: int | None
    baseline_operator_pair_rank: int | None
    promoted_operator_pair_rank: int | None
    baseline_top_level_candidates: tuple[dict[str, Any], ...]
    promoted_actions: tuple[PromotedTwoStepAction, ...]
    baseline_validation_attempt_count: int
    promoted_validation_attempt_count: int
    excluded_heldout_precedent_count: int
    one_step_fallback_preserved: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "case": self.case.to_dict(),
            "baseline_intermediate_rank": self.baseline_intermediate_rank,
            "baseline_operator_pair_rank": self.baseline_operator_pair_rank,
            "promoted_operator_pair_rank": self.promoted_operator_pair_rank,
            "baseline_top_level_candidates": list(self.baseline_top_level_candidates),
            "promoted_actions": [item.to_dict() for item in self.promoted_actions],
            "baseline_validation_attempt_count": (
                self.baseline_validation_attempt_count
            ),
            "promoted_validation_attempt_count": (
                self.promoted_validation_attempt_count
            ),
            "excluded_heldout_precedent_count": (self.excluded_heldout_precedent_count),
            "one_step_fallback_preserved": self.one_step_fallback_preserved,
        }


def _sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def _reaction_product(reaction_smiles: str) -> str:
    if reaction_smiles.count(">>") == 1:
        return reaction_smiles.split(">>", 1)[1]
    if reaction_smiles.count(">") == 2:
        return reaction_smiles.split(">", 2)[2]
    return ""


def _reaction_reactants(reaction_smiles: str) -> str:
    if reaction_smiles.count(">>") == 1:
        return reaction_smiles.split(">>", 1)[0]
    if reaction_smiles.count(">") == 2:
        return reaction_smiles.split(">", 2)[0]
    return ""


def _components(smiles: str) -> tuple[str, ...]:
    canonical = canonical_smiles(smiles)
    return tuple(canonical.split(".")) if canonical else ()


def _scaffold(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)


def _operator_pair(occurrence: Any, steps: Mapping[str, Any]) -> tuple[str, str] | None:
    first = steps.get(occurrence.first_reaction_node_id)
    second = steps.get(occurrence.second_reaction_node_id)
    if first is None or second is None:
        return None
    first_identity = generic_operator_identity_from_observation(
        first.reaction_smiles, first.reaction_signature or {}
    )
    second_identity = generic_operator_identity_from_observation(
        second.reaction_smiles, second.reaction_signature or {}
    )
    if first_identity is None or second_identity is None:
        return None
    return first_identity.operator_id, second_identity.operator_id


def _strategy_id(
    relationship: str, first_operator_id: str, second_operator_id: str
) -> str:
    return digest("CRV1OP1", relationship, first_operator_id, second_operator_id)


def _v1_heldout_candidate_pool(
    route_core_source: str | Path,
    config: CoupledStrategyEvaluationConfig,
) -> tuple[
    dict[str, PromotedV1OperatorPair],
    tuple[CoupledStrategyEvaluationCase, ...],
    dict[str, int],
]:
    """Build all recurrent held-out cases without consulting a library."""

    raw: list[tuple[Any, str, str]] = []
    for projection in iter_route_core_projections(route_core_source):
        steps = {item.reaction_node_id: item for item in projection.steps}
        for occurrence in extract_coupled_route_strategies(projection):
            if occurrence.relationship_class not in V1_ADMITTED_RELATIONSHIPS:
                continue
            pair = _operator_pair(occurrence, steps)
            if pair is not None:
                raw.append((occurrence, pair[0], pair[1]))

    training: dict[tuple[str, str, str], list[Any]] = defaultdict(list)
    heldout: dict[tuple[str, str, str], list[Any]] = defaultdict(list)
    for occurrence, first_operator, second_operator in raw:
        key = (occurrence.relationship_class, first_operator, second_operator)
        if occurrence.split == "train":
            training[key].append(occurrence)
        elif occurrence.split in {"validation", "test"}:
            heldout[key].append(occurrence)

    definitions: dict[str, PromotedV1OperatorPair] = {}
    cases: list[CoupledStrategyEvaluationCase] = []
    heldout_counts: dict[str, int] = {}
    for key, eval_occurrences in heldout.items():
        train_items = training.get(key, [])
        patents = tuple(
            sorted({str(item.patent_id) for item in train_items if item.patent_id})
        )
        if len(patents) < config.minimum_training_patents:
            continue
        relationship, first_operator, second_operator = key
        strategy_id = _strategy_id(relationship, first_operator, second_operator)
        dependencies = Counter(item.dependency_class for item in train_items)
        definition = PromotedV1OperatorPair(
            strategy_id=strategy_id,
            relationship_class=relationship,
            first_operator_id=first_operator,
            second_operator_id=second_operator,
            training_patent_ids=patents,
            training_occurrence_count=len(train_items),
            v2_dependency_counts=tuple(sorted(dependencies.items())),
        )
        definitions[strategy_id] = definition
        heldout_counts[strategy_id] = len(eval_occurrences)
        train_targets = {
            canonical_smiles(item.final_product_smiles)
            for item in train_items
            if canonical_smiles(item.final_product_smiles)
        }
        train_scaffolds = {_scaffold(item) for item in train_targets}
        for occurrence in eval_occurrences:
            target = canonical_smiles(occurrence.final_product_smiles)
            intermediate = canonical_smiles(occurrence.intermediate_smiles)
            terminal = canonical_smiles(
                _reaction_reactants(occurrence.overall_reaction_smiles)
            )
            patent = str(occurrence.patent_id or "")
            if not target or not intermediate or not terminal or not patent:
                continue
            if patent in patents:
                continue
            cases.append(
                CoupledStrategyEvaluationCase(
                    case_id=digest("CRV1CASE1", occurrence.occurrence_id),
                    occurrence_id=occurrence.occurrence_id,
                    strategy_id=strategy_id,
                    patent_id=patent,
                    split=str(occurrence.split),
                    relationship_class=relationship,
                    v2_dependency_class=occurrence.dependency_class,
                    target_smiles=target,
                    expected_intermediate_smiles=intermediate,
                    expected_terminal_precursor_smiles=terminal,
                    observed_first_reaction_smiles=occurrence.first_reaction_smiles,
                    observed_second_reaction_smiles=occurrence.second_reaction_smiles,
                    exact_target_seen_in_training=target in train_targets,
                    target_scaffold_seen_in_training=_scaffold(target)
                    in train_scaffolds,
                )
            )

    return definitions, tuple(cases), heldout_counts


def _select_v1_cases(
    cases: Sequence[CoupledStrategyEvaluationCase],
    config: CoupledStrategyEvaluationConfig,
    *,
    required_strategy_ids: Sequence[str] = (),
    identity_namespace: str = "CRV1PANEL2",
) -> tuple[CoupledStrategyEvaluationCase, ...]:
    ranked = sorted(
        cases,
        key=lambda item: (
            item.exact_target_seen_in_training,
            digest(identity_namespace, str(config.seed), item.case_id),
        ),
    )
    selected: list[CoupledStrategyEvaluationCase] = []
    used_patents: set[str] = set()
    used_strategies: set[str] = set()
    for strategy_id in required_strategy_ids:
        required = next(
            (
                item
                for item in ranked
                if item.strategy_id == strategy_id
                and item.patent_id not in used_patents
            ),
            None,
        )
        if required is None:
            raise ValueError(
                f"required held-out strategy has no selectable case: {strategy_id}"
            )
        selected.append(required)
        used_patents.add(required.patent_id)
        used_strategies.add(required.strategy_id)
    per_class_target = {
        relationship: config.panel_size // 2
        for relationship in sorted(V1_ADMITTED_RELATIONSHIPS)
    }
    for relationship in sorted(V1_ADMITTED_RELATIONSHIPS):
        for item in ranked:
            if (
                len([x for x in selected if x.relationship_class == relationship])
                >= per_class_target[relationship]
            ):
                break
            if (
                item.relationship_class == relationship
                and item.patent_id not in used_patents
                and item.strategy_id not in used_strategies
            ):
                selected.append(item)
                used_patents.add(item.patent_id)
                used_strategies.add(item.strategy_id)
    for item in ranked:
        if len(selected) >= config.panel_size:
            break
        if item.patent_id not in used_patents:
            selected.append(item)
            used_patents.add(item.patent_id)
            used_strategies.add(item.strategy_id)
    selected.sort(key=lambda item: item.case_id)
    return tuple(sorted(selected[: config.panel_size], key=lambda item: item.case_id))


def _library_gaps(
    definitions: Iterable[PromotedV1OperatorPair],
    library: GenericTemplateLibrary,
    heldout_counts: Mapping[str, int],
) -> tuple[dict[str, Any], ...]:
    library_operators = {item.operator_id for item in library.templates}
    gaps = []
    for definition in definitions:
        missing = tuple(
            operator
            for operator in (
                definition.first_operator_id,
                definition.second_operator_id,
            )
            if operator not in library_operators
        )
        if missing:
            gaps.append(
                {
                    "strategy_id": definition.strategy_id,
                    "relationship_class": definition.relationship_class,
                    "missing_operator_ids": list(missing),
                    "heldout_occurrence_count": int(
                        heldout_counts.get(definition.strategy_id) or 0
                    ),
                }
            )
    return tuple(sorted(gaps, key=lambda item: item["strategy_id"]))


def build_frozen_v1_heldout_panel(
    route_core_source: str | Path,
    *,
    config: CoupledStrategyEvaluationConfig | None = None,
    required_strategy_ids: Sequence[str] = (),
) -> FrozenV1HeldoutPanel:
    """Build a deterministic panel whose selection is library-independent."""

    resolved = config or CoupledStrategyEvaluationConfig()
    source = Path(route_core_source)
    definitions, cases, _ = _v1_heldout_candidate_pool(source, resolved)
    required = tuple(dict.fromkeys(str(item) for item in required_strategy_ids))
    if len(required) > resolved.panel_size:
        raise ValueError("required strategies exceed frozen panel size")
    selected = _select_v1_cases(
        cases,
        resolved,
        required_strategy_ids=required,
    )
    selected_definitions = tuple(
        definitions[key] for key in sorted({item.strategy_id for item in selected})
    )
    route_sha256 = _sha256(source)
    panel_id = digest(
        "CRV1PANEL2",
        route_sha256,
        json.dumps(asdict(resolved), sort_keys=True, separators=(",", ":")),
        json.dumps(
            [item.case_id for item in selected],
            sort_keys=True,
            separators=(",", ":"),
        ),
    )
    return FrozenV1HeldoutPanel(
        panel_id=panel_id,
        route_core_source=str(source.resolve()),
        route_core_sha256=route_sha256,
        config=resolved,
        strategies=selected_definitions,
        cases=selected,
        required_strategy_ids=required,
    )


def write_frozen_v1_heldout_panel(
    panel: FrozenV1HeldoutPanel,
    output_path: str | Path,
) -> dict[str, Any]:
    """Write one deterministic library-independent panel artifact."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(panel.to_dict(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "output_json": str(destination.resolve()),
        "panel_id": panel.panel_id,
        "panel_case_count": len(panel.cases),
        "strategy_count": len(panel.strategies),
    }


def load_frozen_v1_heldout_panel(
    source: str | Path,
) -> FrozenV1HeldoutPanel:
    """Load and validate one frozen v1 held-out panel."""

    value = json.loads(Path(source).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("frozen v1 panel must be an object")
    return FrozenV1HeldoutPanel.from_dict(value)


def build_v1_heldout_panel(
    route_core_source: str | Path,
    library: GenericTemplateLibrary,
    *,
    config: CoupledStrategyEvaluationConfig | None = None,
) -> tuple[
    tuple[PromotedV1OperatorPair, ...],
    tuple[CoupledStrategyEvaluationCase, ...],
    tuple[dict[str, Any], ...],
]:
    """Build the legacy library-covered deterministic held-out panel."""

    resolved = config or CoupledStrategyEvaluationConfig()
    definitions, cases, heldout_counts = _v1_heldout_candidate_pool(
        route_core_source,
        resolved,
    )
    gaps = _library_gaps(definitions.values(), library, heldout_counts)
    gap_ids = {str(item["strategy_id"]) for item in gaps}
    eligible_cases = tuple(item for item in cases if item.strategy_id not in gap_ids)
    selected = _select_v1_cases(
        eligible_cases,
        resolved,
        identity_namespace="CRV1PANEL1",
    )
    selected_definitions = tuple(
        definitions[key] for key in sorted({item.strategy_id for item in selected})
    )
    return selected_definitions, selected, gaps


def build_v1_operator_pair_inventory(
    route_core_source: str | Path,
    library: GenericTemplateLibrary,
    *,
    minimum_training_patents: int = 2,
    seed: int = 20260818,
) -> tuple[
    tuple[PromotedV1OperatorPair, ...],
    tuple[CoupledStrategyEvaluationCase, ...],
    tuple[dict[str, Any], ...],
]:
    """Return every recurrent patent-disjoint v1 pair covered by a library.

    Unlike the held-out panel builder, this inventory does not sample or
    balance strategies. It is intended for complete review catalogs and keeps
    every eligible example for representative-example selection.
    """

    config = CoupledStrategyEvaluationConfig(
        panel_size=1,
        minimum_training_patents=minimum_training_patents,
        seed=seed,
    )
    definitions, cases, heldout_counts = _v1_heldout_candidate_pool(
        route_core_source,
        config,
    )
    gaps = _library_gaps(definitions.values(), library, heldout_counts)
    gap_ids = {str(item["strategy_id"]) for item in gaps}
    eligible_ids = set(definitions) - gap_ids
    eligible_definitions = tuple(
        definitions[strategy_id] for strategy_id in sorted(eligible_ids)
    )
    eligible_cases = tuple(
        sorted(
            (case for case in cases if case.strategy_id in eligible_ids),
            key=lambda item: (
                item.strategy_id,
                item.patent_id,
                item.case_id,
            ),
        )
    )
    return eligible_definitions, eligible_cases, gaps


def _patent_overlap(candidate: GenericDisconnectionCandidate, patent_id: str) -> bool:
    patent = _PATENT_TOKEN.sub("", patent_id.upper())
    if not patent:
        return False
    return any(
        patent in _PATENT_TOKEN.sub("", value.upper())
        for value in candidate.precedent_reaction_ids
    )


def _without_heldout_precedents(
    candidates: Iterable[GenericDisconnectionCandidate], patent_id: str
) -> tuple[tuple[GenericDisconnectionCandidate, ...], int]:
    values = tuple(candidates)
    kept = tuple(item for item in values if not _patent_overlap(item, patent_id))
    return kept, len(values) - len(kept)


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


def _search(
    searcher: SearchFunction,
    target: str,
    library: GenericTemplateLibrary,
    config: CoupledStrategyEvaluationConfig,
    *,
    operator_ids: Sequence[str] = (),
) -> tuple[tuple[GenericDisconnectionCandidate, ...], GenericSearchDiagnostics]:
    return searcher(
        target,
        library,
        operator_ids=operator_ids,
        top_k=config.top_k,
        max_templates_to_apply=config.max_templates_to_apply,
        max_candidates_to_validate=config.max_candidates_to_validate,
        use_context=True,
    )


def search_promoted_v1_strategies(
    target_smiles: str,
    library: GenericTemplateLibrary,
    strategies: Iterable[PromotedV1OperatorPair],
    *,
    top_k: int = 5,
    max_templates_to_apply: int = 50,
    max_candidates_to_validate: int = 12,
    include_l0: bool = True,
    use_context: bool = True,
    include_one_step_fallbacks: bool = True,
    searcher: SearchFunction = disconnect_generic_target_detailed,
) -> CoupledStrategyQueryResult:
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
    generated: list[tuple[PromotedV1OperatorPair, PromotedTwoStepAction]] = []
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
                    support = min(
                        1.0,
                        math.log1p(len(strategy.training_patent_ids)) / math.log(11),
                    )
                    generated.append(
                        (
                            strategy,
                            PromotedTwoStepAction(
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
                                    0.45 * first.score
                                    + 0.45 * second.score
                                    + 0.10 * support,
                                    8,
                                ),
                            ),
                        )
                    )

    unique: dict[
        tuple[str, str, str, str],
        tuple[PromotedV1OperatorPair, PromotedTwoStepAction],
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
        PromotedTwoStepQueryAction(
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
    return CoupledStrategyQueryResult(
        target_smiles=canonical_target,
        actions=actions,
        one_step_fallbacks=fallback_values,
        diagnostics=CoupledStrategyQueryDiagnostics(
            strategy_count=len(definitions),
            capable_strategy_count=capable,
            capability_gap_count=capability_gaps,
            second_step_validation_attempt_count=second_attempts,
            first_step_validation_attempt_count=first_attempts,
            fallback_validation_attempt_count=fallback_attempts,
            generated_action_count=len(unique),
            returned_action_count=len(actions),
            returned_fallback_count=len(fallback_values),
        ),
        warnings=warnings,
    )


def evaluate_v1_case(
    case: CoupledStrategyEvaluationCase,
    strategy: PromotedV1OperatorPair,
    library: GenericTemplateLibrary,
    *,
    config: CoupledStrategyEvaluationConfig | None = None,
    searcher: SearchFunction = disconnect_generic_target_detailed,
) -> CoupledStrategyCaseResult:
    """Compare ordinary depth-two search with one promoted operator pair."""

    resolved = config or CoupledStrategyEvaluationConfig()
    if case.strategy_id != strategy.strategy_id:
        raise ValueError("case and strategy identities do not match")
    excluded = 0
    baseline_top, baseline_diagnostics = _search(
        searcher, case.target_smiles, library, resolved
    )
    baseline_top, count = _without_heldout_precedents(baseline_top, case.patent_id)
    excluded += count
    baseline_intermediate_rank = next(
        (
            rank
            for rank, candidate in enumerate(baseline_top, 1)
            if case.expected_intermediate_smiles
            in _components(candidate.precursor_smiles)
        ),
        None,
    )
    baseline_depth_two_paths: list[tuple[float, str, str, bool]] = []
    baseline_attempts = baseline_diagnostics.validation_attempt_count
    for second in baseline_top:
        for intermediate in _components(second.precursor_smiles):
            first_candidates, diagnostics = _search(
                searcher, intermediate, library, resolved
            )
            baseline_attempts += diagnostics.validation_attempt_count
            first_candidates, count = _without_heldout_precedents(
                first_candidates, case.patent_id
            )
            excluded += count
            for first in first_candidates:
                baseline_depth_two_paths.append(
                    (
                        first.score + second.score,
                        intermediate,
                        first.precursor_smiles,
                        (
                            first.operator_id == strategy.first_operator_id
                            and second.operator_id == strategy.second_operator_id
                            and intermediate == case.expected_intermediate_smiles
                        ),
                    )
                )
    baseline_depth_two_paths.sort(key=lambda item: (-item[0], item[1], item[2]))
    baseline_pair_rank = next(
        (rank for rank, item in enumerate(baseline_depth_two_paths, 1) if item[3]),
        None,
    )

    promoted_second, promoted_second_diagnostics = _search(
        searcher,
        case.target_smiles,
        library,
        resolved,
        operator_ids=(strategy.second_operator_id,),
    )
    promoted_second, count = _without_heldout_precedents(
        promoted_second, case.patent_id
    )
    excluded += count
    promoted_attempts = promoted_second_diagnostics.validation_attempt_count
    actions: list[PromotedTwoStepAction] = []
    for second in promoted_second:
        for intermediate in _components(second.precursor_smiles):
            first_candidates, diagnostics = _search(
                searcher,
                intermediate,
                library,
                resolved,
                operator_ids=(strategy.first_operator_id,),
            )
            promoted_attempts += diagnostics.validation_attempt_count
            first_candidates, count = _without_heldout_precedents(
                first_candidates, case.patent_id
            )
            excluded += count
            for first in first_candidates:
                terminal = _merge_terminal_precursors(
                    first.precursor_smiles,
                    second.precursor_smiles,
                    intermediate,
                )
                if terminal is None:
                    continue
                support = min(
                    1.0,
                    math.log1p(len(strategy.training_patent_ids)) / math.log(11),
                )
                actions.append(
                    PromotedTwoStepAction(
                        strategy_id=strategy.strategy_id,
                        intermediate_smiles=intermediate,
                        terminal_precursor_smiles=terminal,
                        first_operator_id=first.operator_id,
                        second_operator_id=second.operator_id,
                        first_reaction_smiles=first.proposed_reaction_smiles,
                        second_reaction_smiles=second.proposed_reaction_smiles,
                        first_forward_validation_status=(
                            first.forward_validation_status
                        ),
                        second_forward_validation_status=(
                            second.forward_validation_status
                        ),
                        score=round(
                            0.45 * first.score + 0.45 * second.score + 0.10 * support,
                            8,
                        ),
                    )
                )
    unique = {
        (
            item.intermediate_smiles,
            item.terminal_precursor_smiles,
            item.first_reaction_smiles,
            item.second_reaction_smiles,
        ): item
        for item in actions
    }
    ranked_actions = tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                -item.score,
                item.intermediate_smiles,
                item.terminal_precursor_smiles,
            ),
        )[: resolved.top_k]
    )
    promoted_rank = next(
        (
            rank
            for rank, item in enumerate(ranked_actions, 1)
            if item.intermediate_smiles == case.expected_intermediate_smiles
        ),
        None,
    )
    return CoupledStrategyCaseResult(
        case=case,
        baseline_intermediate_rank=baseline_intermediate_rank,
        baseline_operator_pair_rank=baseline_pair_rank,
        promoted_operator_pair_rank=promoted_rank,
        baseline_top_level_candidates=tuple(item.to_dict() for item in baseline_top),
        promoted_actions=ranked_actions,
        baseline_validation_attempt_count=baseline_attempts,
        promoted_validation_attempt_count=promoted_attempts,
        excluded_heldout_precedent_count=excluded,
    )


def run_v1_coupled_strategy_evaluation(
    route_core_source: str | Path,
    library: GenericTemplateLibrary,
    *,
    operator_library_source: str | Path | None = None,
    config: CoupledStrategyEvaluationConfig | None = None,
    frozen_panel: FrozenV1HeldoutPanel | None = None,
    searcher: SearchFunction = disconnect_generic_target_detailed,
) -> dict[str, Any]:
    """Build the held-out panel and run the paired real-example benchmark."""

    resolved = config or (
        frozen_panel.config
        if frozen_panel is not None
        else CoupledStrategyEvaluationConfig()
    )
    route_source = Path(route_core_source)
    route_sha256 = _sha256(route_source)
    if frozen_panel is None:
        definitions, cases, gaps = build_v1_heldout_panel(
            route_source, library, config=resolved
        )
    else:
        if frozen_panel.route_core_sha256 != route_sha256:
            raise ValueError("frozen panel route-core hash does not match source")
        definitions = frozen_panel.strategies
        cases = frozen_panel.cases
        selected_counts = Counter(item.strategy_id for item in cases)
        gaps = _library_gaps(definitions, library, selected_counts)
    by_id = {item.strategy_id: item for item in definitions}
    results = tuple(
        evaluate_v1_case(
            case,
            by_id[case.strategy_id],
            library,
            config=resolved,
            searcher=searcher,
        )
        for case in cases
    )
    baseline_hits = sum(
        item.baseline_operator_pair_rank is not None for item in results
    )
    promoted_hits = sum(
        item.promoted_operator_pair_rank is not None for item in results
    )
    return {
        "artifact_type": "v1_coupled_strategy_heldout_evaluation",
        "schema_version": COUPLED_STRATEGY_EVALUATION_SCHEMA_VERSION,
        "algorithm_version": COUPLED_STRATEGY_EVALUATION_ALGORITHM_VERSION,
        "policy": {
            "base": "v1_structural_relationship",
            "admitted_relationships": sorted(V1_ADMITTED_RELATIONSHIPS),
            "v2_role": "review_metadata_only",
            "physical_step_validation": "mandatory_for_both_steps",
            "fallback": "ordinary_one_step_candidates_preserved",
        },
        "config": asdict(resolved),
        "route_core_source": str(route_core_source),
        "route_core_sha256": route_sha256,
        "frozen_panel_id": frozen_panel.panel_id if frozen_panel else None,
        "panel_selection": (
            "library_independent_frozen"
            if frozen_panel is not None
            else "library_covered_legacy"
        ),
        "operator_library_source": (
            str(operator_library_source) if operator_library_source else None
        ),
        "operator_library_sha256": (
            _sha256(Path(operator_library_source)) if operator_library_source else None
        ),
        "panel_case_count": len(cases),
        "strategy_count": len(definitions),
        "relationship_counts": dict(
            sorted(Counter(item.relationship_class for item in cases).items())
        ),
        "metrics": {
            "baseline_pair_recall": (baseline_hits / len(results) if results else 0.0),
            "promoted_pair_recall": (promoted_hits / len(results) if results else 0.0),
            "baseline_pair_hit_count": baseline_hits,
            "promoted_pair_hit_count": promoted_hits,
            "baseline_validation_attempt_count": sum(
                item.baseline_validation_attempt_count for item in results
            ),
            "promoted_validation_attempt_count": sum(
                item.promoted_validation_attempt_count for item in results
            ),
            "all_fallbacks_preserved": all(
                item.one_step_fallback_preserved for item in results
            ),
        },
        "strategies": [item.to_dict() for item in definitions],
        "capability_gaps": list(gaps),
        "results": [item.to_dict() for item in results],
    }


def write_v1_coupled_strategy_evaluation(
    report: Mapping[str, Any], output_path: str | Path
) -> dict[str, Any]:
    """Write one deterministic evaluation artifact."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(dict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return {
        "output_json": str(destination.resolve()),
        "json_bytes": destination.stat().st_size,
        "panel_case_count": int(report.get("panel_case_count") or 0),
    }


__all__ = [
    "COUPLED_STRATEGY_EVALUATION_ALGORITHM_VERSION",
    "COUPLED_STRATEGY_EVALUATION_SCHEMA_VERSION",
    "FROZEN_V1_PANEL_ALGORITHM_VERSION",
    "FROZEN_V1_PANEL_SCHEMA_VERSION",
    "CoupledStrategyCaseResult",
    "CoupledStrategyEvaluationCase",
    "CoupledStrategyEvaluationConfig",
    "CoupledStrategyQueryDiagnostics",
    "CoupledStrategyQueryResult",
    "FrozenV1HeldoutPanel",
    "PromotedTwoStepAction",
    "PromotedTwoStepQueryAction",
    "PromotedV1OperatorPair",
    "build_frozen_v1_heldout_panel",
    "build_v1_operator_pair_inventory",
    "build_v1_heldout_panel",
    "evaluate_v1_case",
    "load_frozen_v1_heldout_panel",
    "run_v1_coupled_strategy_evaluation",
    "search_promoted_v1_strategies",
    "write_frozen_v1_heldout_panel",
    "write_v1_coupled_strategy_evaluation",
]
