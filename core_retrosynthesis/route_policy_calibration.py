"""Whole-route validation gate for learned route-action policy influence."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Sequence

from .chemistry import digest
from .generic_models import GenericTemplateLibrary
from .multistep import (
    LiteratureLookup,
    MultistepRetrosynthesisResult,
    plan_multistep_routes,
)
from .route_action_evaluation import RouteActionEvaluation
from .route_action_policy import (
    RouteActionPolicyModel,
    save_route_action_policy,
)


ROUTE_POLICY_CALIBRATION_SCHEMA_VERSION = "1.0"
RoutePlanner = Callable[..., MultistepRetrosynthesisResult]
ProgressCallback = Callable[[str], None]


@dataclass(frozen=True)
class RoutePolicyCalibrationConfig:
    """Fixed whole-route search and activation settings."""

    max_depth: int = 3
    molecular_weight_threshold: float = 80.0
    top_k_routes: int = 3
    per_step_top_k: int = 4
    beam_width: int = 15
    max_expansions: int = 15
    max_templates_to_apply: int = 100
    max_candidates_to_validate: int = 20
    minimum_validation_targets: int = 5

    def __post_init__(self) -> None:
        if self.max_depth not in {2, 3}:
            raise ValueError("calibration route depth must be 2 or 3")
        if self.molecular_weight_threshold <= 0.0:
            raise ValueError("calibration molecular-weight threshold must be positive")
        for value, label in (
            (self.top_k_routes, "top-k routes"),
            (self.per_step_top_k, "per-step top-k"),
            (self.beam_width, "beam width"),
            (self.max_expansions, "maximum expansions"),
            (self.max_templates_to_apply, "maximum templates"),
            (self.max_candidates_to_validate, "maximum validations"),
            (self.minimum_validation_targets, "minimum validation targets"),
        ):
            if value < 1:
                raise ValueError(f"calibration {label} must be positive")
        if self.max_candidates_to_validate < self.per_step_top_k:
            raise ValueError("calibration validation budget must cover per-step top-k")

    @property
    def config_id(self) -> str:
        """Return a content-derived identity for the route search budget."""

        payload = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return digest("RPCC1", payload)

    def planner_kwargs(self) -> dict[str, Any]:
        """Return only arguments owned by multistep route search."""

        values = asdict(self)
        values.pop("minimum_validation_targets")
        return values


@dataclass(frozen=True)
class RoutePolicyCalibrationTarget:
    """One validation route selected independently of policy outcomes."""

    route_id: str
    target_smiles: str
    observed_strategy_ids: tuple[str, ...]
    reaction_count: int
    maximum_depth: int


@dataclass(frozen=True)
class RoutePolicyTargetOutcome:
    """Planner outcome for one target at one policy residual scale."""

    route_id: str
    solved: bool
    observed_strategy_count: int
    observed_strategy_recovered_count: int
    best_terminal_fraction: float
    expanded_states: int
    first_solution_expansion: Optional[int]
    policy_reordered_expansions: int


@dataclass(frozen=True)
class RoutePolicyScaleOutcome:
    """Aggregated whole-route validation evidence for one residual scale."""

    residual_scale: float
    target_count: int
    solved_target_count: int
    observed_strategy_recovered_count: int
    mean_best_terminal_fraction: float
    first_solution_efficiency: int
    expanded_state_count: int
    policy_reordered_expansion_count: int
    targets: tuple[RoutePolicyTargetOutcome, ...]

    @property
    def selection_key(self) -> tuple[int, int, float, int, int, float]:
        """Return the validation objective, with lower influence winning ties."""

        return (
            self.solved_target_count,
            self.observed_strategy_recovered_count,
            round(self.mean_best_terminal_fraction, 8),
            self.first_solution_efficiency,
            -self.expanded_state_count,
            -self.residual_scale,
        )


@dataclass(frozen=True)
class RoutePolicyCalibrationReport:
    """Frozen route-level activation decision and its complete evidence."""

    config_id: str
    input_model_id: str
    output_model_id: str
    selected_residual_scale: float
    policy_active: bool
    activation_reason: str
    target_ids: tuple[str, ...]
    scale_outcomes: tuple[RoutePolicyScaleOutcome, ...]
    schema_version: str = ROUTE_POLICY_CALIBRATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible calibration report."""

        return asdict(self)


def build_route_policy_calibration_targets(
    evaluations: Iterable[RouteActionEvaluation],
    target_ids: Sequence[str],
    *,
    required_split: str = "validation",
) -> tuple[RoutePolicyCalibrationTarget, ...]:
    """Resolve a fixed panel to verified validation-route supervision."""

    by_route_id = {
        evaluation.source_route_id or evaluation.tree_id: evaluation
        for evaluation in evaluations
    }
    values = []
    for route_id in target_ids:
        evaluation = by_route_id.get(route_id)
        if evaluation is None:
            raise ValueError(f"calibration route is absent from replay: {route_id}")
        if evaluation.split != required_split:
            raise ValueError(
                f"calibration route {route_id} is not in {required_split} split"
            )
        if any(step.source_patent_precedent_overlap for step in evaluation.steps):
            raise ValueError(
                f"calibration route {route_id} has source-patent precedent overlap"
            )
        strategy_ids = tuple(
            sorted(
                {
                    str(step.observed_action.strategy_id)
                    for step in evaluation.steps
                    if step.observed_action.strategy_verified
                    and step.observed_action.strategy_id
                }
            )
        )
        values.append(
            RoutePolicyCalibrationTarget(
                route_id=route_id,
                target_smiles=evaluation.target_smiles,
                observed_strategy_ids=strategy_ids,
                reaction_count=evaluation.reaction_count,
                maximum_depth=evaluation.maximum_depth,
            )
        )
    if len({target.route_id for target in values}) != len(values):
        raise ValueError("calibration target IDs must be unique")
    return tuple(values)


def _returned_strategy_ids(result: MultistepRetrosynthesisResult) -> set[str]:
    routes = result.routes or result.partial_routes
    return {
        step.candidate.strategy_id
        for route in routes
        for step in route.steps
        if step.candidate.strategy_id
    }


def _best_terminal_fraction(result: MultistepRetrosynthesisResult) -> float:
    routes = result.routes or result.partial_routes
    fractions = [
        sum(leaf.terminal for leaf in route.leaves) / len(route.leaves)
        for route in routes
        if route.leaves
    ]
    return max(fractions, default=0.0)


def select_route_policy_scale(
    outcomes: Sequence[RoutePolicyScaleOutcome],
) -> RoutePolicyScaleOutcome:
    """Select whole-route validation performance, preferring zero on ties."""

    if not outcomes:
        raise ValueError("route-policy calibration requires scale outcomes")
    if not any(outcome.residual_scale == 0.0 for outcome in outcomes):
        raise ValueError("route-policy calibration requires a zero-scale baseline")
    return max(outcomes, key=lambda outcome: outcome.selection_key)


def calibrate_route_action_policy(
    model: RouteActionPolicyModel,
    targets: Sequence[RoutePolicyCalibrationTarget],
    library: GenericTemplateLibrary,
    literature_index: LiteratureLookup,
    *,
    config: RoutePolicyCalibrationConfig = RoutePolicyCalibrationConfig(),
    planner: RoutePlanner = plan_multistep_routes,
    progress: Optional[ProgressCallback] = None,
) -> tuple[RouteActionPolicyModel, RoutePolicyCalibrationReport]:
    """Freeze policy influence using only whole-route validation outcomes."""

    values = tuple(targets)
    scales = model.definition.validation_residual_scales
    outcomes = []
    for scale in scales:
        if progress is not None:
            progress(f"calibrating residual scale {scale:g}")
        scaled_model = RouteActionPolicyModel(
            definition=model.definition,
            weights=model.weights,
            residual_scale=scale,
        )
        target_outcomes = []
        for target in values:
            if progress is not None:
                progress(f"  planning {target.route_id}")
            result = planner(
                target.target_smiles,
                library,
                literature_index,
                route_action_policy=scaled_model,
                **config.planner_kwargs(),
            )
            recovered = set(target.observed_strategy_ids).intersection(
                _returned_strategy_ids(result)
            )
            target_outcomes.append(
                RoutePolicyTargetOutcome(
                    route_id=target.route_id,
                    solved=bool(result.routes),
                    observed_strategy_count=len(target.observed_strategy_ids),
                    observed_strategy_recovered_count=len(recovered),
                    best_terminal_fraction=round(_best_terminal_fraction(result), 8),
                    expanded_states=result.diagnostics.expanded_states,
                    first_solution_expansion=(
                        result.diagnostics.first_solution_expansion
                    ),
                    policy_reordered_expansions=(
                        result.diagnostics.route_policy_reordered_expansions
                    ),
                )
            )
            if progress is not None:
                progress(
                    "    solved="
                    f"{bool(result.routes)} expansions="
                    f"{result.diagnostics.expanded_states} recovered="
                    f"{len(recovered)}/{len(target.observed_strategy_ids)}"
                )
        count = len(target_outcomes)
        outcomes.append(
            RoutePolicyScaleOutcome(
                residual_scale=scale,
                target_count=count,
                solved_target_count=sum(item.solved for item in target_outcomes),
                observed_strategy_recovered_count=sum(
                    item.observed_strategy_recovered_count
                    for item in target_outcomes
                ),
                mean_best_terminal_fraction=round(
                    sum(item.best_terminal_fraction for item in target_outcomes)
                    / count
                    if count
                    else 0.0,
                    8,
                ),
                first_solution_efficiency=sum(
                    config.max_expansions + 1 - item.first_solution_expansion
                    for item in target_outcomes
                    if item.first_solution_expansion is not None
                ),
                expanded_state_count=sum(
                    item.expanded_states for item in target_outcomes
                ),
                policy_reordered_expansion_count=sum(
                    item.policy_reordered_expansions for item in target_outcomes
                ),
                targets=tuple(target_outcomes),
            )
        )
    if len(values) < config.minimum_validation_targets:
        selected_scale = 0.0
        reason = (
            "insufficient_route_validation_targets:"
            f"{len(values)}<{config.minimum_validation_targets}"
        )
    else:
        selected_scale = select_route_policy_scale(outcomes).residual_scale
        reason = (
            "whole_route_validation_selected_residual"
            if selected_scale > 0.0
            else "whole_route_validation_selected_baseline"
        )
    calibrated = RouteActionPolicyModel(
        definition=model.definition,
        weights=model.weights,
        residual_scale=selected_scale,
    )
    report = RoutePolicyCalibrationReport(
        config_id=config.config_id,
        input_model_id=model.model_id,
        output_model_id=calibrated.model_id,
        selected_residual_scale=selected_scale,
        policy_active=selected_scale > 0.0,
        activation_reason=reason,
        target_ids=tuple(target.route_id for target in values),
        scale_outcomes=tuple(outcomes),
    )
    return calibrated, report


def save_route_policy_calibration(
    model: RouteActionPolicyModel,
    report: RoutePolicyCalibrationReport,
    output_model: str | Path,
    output_report: str | Path,
    *,
    source_replay: str | Path,
    input_model: str | Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Persist a frozen model and provenance-rich calibration report."""

    model_path = Path(output_model)
    report_path = Path(output_report)
    if not overwrite:
        for path in (model_path, report_path):
            if path.exists():
                raise FileExistsError(path)
    save_route_action_policy(model, model_path)
    replay_path = Path(source_replay)
    input_path = Path(input_model)
    payload = {
        **report.to_dict(),
        "source": {
            "path": str(replay_path.resolve()),
            "sha256": hashlib.sha256(replay_path.read_bytes()).hexdigest(),
        },
        "input_model": {
            "path": str(input_path.resolve()),
            "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        },
        "output_model": {
            "path": str(model_path.resolve()),
            "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
        },
        "warnings": [
            "Calibration changes only policy influence over graph-validated actions.",
            "The fixed validation panel must be selected without policy outcomes.",
            "Patent/scaffold and close-analogue library leakage remain uncontrolled.",
        ],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return payload


__all__ = [
    "ROUTE_POLICY_CALIBRATION_SCHEMA_VERSION",
    "RoutePolicyCalibrationConfig",
    "RoutePolicyCalibrationReport",
    "RoutePolicyCalibrationTarget",
    "RoutePolicyScaleOutcome",
    "RoutePolicyTargetOutcome",
    "build_route_policy_calibration_targets",
    "calibrate_route_action_policy",
    "save_route_policy_calibration",
    "select_route_policy_scale",
]
