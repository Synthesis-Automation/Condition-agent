"""Deterministic, precedent-backed condition repair for endpoint selectivity.

The module never authors a condition or molecular structure. It evaluates only
canonical recipes already recommended for a validated route step and trains the
small choice model only on admitted reaction/reference contexts carried by those
recommendations. Unsupported, ambiguous, and conflicting evidence remains an
explicit non-actionable assessment.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from .selectivity_poc import (
    ConditionalEditChoiceModel,
    ReactionChoiceSet,
    SelectivityAssessment,
    build_reaction_choice_set,
)


CONDITION_SELECTIVITY_REPAIR_SCHEMA_VERSION = "condition_selectivity_repair.v1"
_DEFINITION_PATH = (
    Path(__file__).with_name("definitions")
    / "condition_selectivity_repair.v1.json"
)

ConditionSelectivityRepairStatus = Literal[
    "supported",
    "competing_outcome_favored",
    "ambiguous",
    "insufficient_precedent",
    "condition_conflict",
    "no_condition_evidence",
    "unsupported_topology",
]


def _stable_id(prefix: str, value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    return f"{prefix}:{hashlib.sha256(encoded).hexdigest()[:20]}"


@lru_cache(maxsize=1)
def load_condition_selectivity_repair_policy() -> dict[str, Any]:
    """Load and strictly validate the versioned repair thresholds."""

    value = json.loads(_DEFINITION_PATH.read_text(encoding="utf-8"))
    required = {
        "definition_id",
        "definition_version",
        "schema_version",
        "maximum_recommendations_assessed",
        "minimum_independent_references",
        "minimum_compatibility_score",
        "minimum_desired_probability",
        "minimum_probability_margin",
        "maximum_normalized_entropy",
        "model_feature_dimension",
        "model_epochs",
        "model_learning_rate",
        "require_direct_condition_retrieval",
    }
    if set(value) != required:
        raise ValueError("condition-selectivity repair definition fields are invalid")
    if value["definition_id"] != "condition_selectivity_repair.v1":
        raise ValueError("unexpected condition-selectivity repair definition ID")
    if value["schema_version"] != "1.0":
        raise ValueError("unsupported condition-selectivity repair schema")
    for field in (
        "maximum_recommendations_assessed",
        "minimum_independent_references",
        "model_feature_dimension",
        "model_epochs",
    ):
        if not isinstance(value[field], int) or value[field] < 1:
            raise ValueError(f"{field} must be a positive integer")
    for field in (
        "minimum_compatibility_score",
        "minimum_desired_probability",
        "maximum_normalized_entropy",
    ):
        number = float(value[field])
        if not 0.0 <= number <= 1.0:
            raise ValueError(f"{field} must be in [0, 1]")
    if float(value["minimum_probability_margin"]) < 0.0:
        raise ValueError("minimum_probability_margin cannot be negative")
    if float(value["model_learning_rate"]) <= 0.0:
        raise ValueError("model_learning_rate must be positive")
    if not isinstance(value["require_direct_condition_retrieval"], bool):
        raise ValueError("require_direct_condition_retrieval must be boolean")
    return dict(value)


@dataclass(frozen=True)
class ConditionSelectivityRepairAssessment:
    """One recipe-specific decision over the intended and competing endpoints."""

    assessment_id: str
    step_id: str
    warning_ids: tuple[str, ...]
    recipe_id: str
    recipe_core_id: str
    recommendation_rank: int | None
    status: ConditionSelectivityRepairStatus
    actionable: bool
    reason: str
    condition_evidence_status: str
    compatibility_score: float | None
    desired_candidate_id: str | None
    desired_probability: float | None
    best_competitor_probability: float | None
    probability_margin: float | None
    normalized_entropy: float | None
    exact_condition_reference_ids: tuple[str, ...]
    training_choice_set_count: int
    ranked_outcomes: tuple[dict[str, Any], ...]
    policy_definition_id: str
    policy_definition_version: str
    model_id: str | None = None
    schema_version: str = CONDITION_SELECTIVITY_REPAIR_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return an auditable JSON-compatible assessment."""

        return asdict(self)


def _warning_ids(step: Any) -> tuple[str, ...]:
    return tuple(
        f"{warning.audit_id}:{warning.code}"
        for warning in getattr(step.candidate, "selectivity_warnings", ())
    )


def _recommendation_contexts(
    recommendation: Mapping[str, Any],
) -> tuple[Mapping[str, Any], ...]:
    return tuple(
        value
        for value in (recommendation.get("precedent_reaction_contexts") or ())
        if isinstance(value, Mapping)
    )


def _training_choice_sets(
    recommendations: Sequence[Mapping[str, Any]],
) -> tuple[ReactionChoiceSet, ...]:
    rows: list[ReactionChoiceSet] = []
    seen: set[tuple[str, str, str]] = set()
    for recommendation in recommendations:
        recipe = recommendation.get("resolved_recipe")
        if not isinstance(recipe, Mapping):
            continue
        recipe_id = str(recommendation.get("recipe_id") or "")
        for context in _recommendation_contexts(recommendation):
            reaction_smiles = str(context.get("reaction_smiles") or "")
            reference_id = str(context.get("reference_id") or "")
            identity = (recipe_id, reaction_smiles, reference_id)
            if not reaction_smiles or identity in seen:
                continue
            seen.add(identity)
            try:
                choice_set = build_reaction_choice_set(
                    reaction_smiles,
                    recipe,
                    reference_id=reference_id,
                    label_strength=1.0,
                )
            except ValueError:
                continue
            if len(choice_set.candidates) >= 2:
                rows.append(choice_set)
    return tuple(sorted(rows, key=lambda value: value.choice_set_id))


def _direct_reference_ids(
    query: ReactionChoiceSet,
    training: Sequence[ReactionChoiceSet],
) -> tuple[str, ...]:
    selected_signature = query.selected_candidate.site_signature
    return tuple(
        sorted(
            {
                row.reference_id
                for row in training
                if row.reference_id
                and row.condition_tokens == query.condition_tokens
                and row.selected_candidate.site_signature == selected_signature
            }
        )
    )


def _ranked_outcomes(
    assessment: SelectivityAssessment,
) -> tuple[dict[str, Any], ...]:
    return tuple(asdict(value) for value in assessment.ranked_outcomes)


def _empty_assessment(
    step: Any,
    *,
    status: ConditionSelectivityRepairStatus,
    reason: str,
    condition_status: str,
) -> ConditionSelectivityRepairAssessment:
    policy = load_condition_selectivity_repair_policy()
    identity = {
        "step_id": step.step_id,
        "warning_ids": _warning_ids(step),
        "status": status,
        "reason": reason,
        "policy": policy["definition_id"],
    }
    return ConditionSelectivityRepairAssessment(
        assessment_id=_stable_id("CSRA1", identity),
        step_id=step.step_id,
        warning_ids=_warning_ids(step),
        recipe_id="",
        recipe_core_id="",
        recommendation_rank=None,
        status=status,
        actionable=False,
        reason=reason,
        condition_evidence_status=condition_status,
        compatibility_score=None,
        desired_candidate_id=None,
        desired_probability=None,
        best_competitor_probability=None,
        probability_margin=None,
        normalized_entropy=None,
        exact_condition_reference_ids=(),
        training_choice_set_count=0,
        ranked_outcomes=(),
        policy_definition_id=str(policy["definition_id"]),
        policy_definition_version=str(policy["definition_version"]),
    )


def assess_condition_selectivity_repairs(
    step: Any,
) -> tuple[ConditionSelectivityRepairAssessment, ...]:
    """Assess existing recommended recipes for one selectivity-warning step."""

    warnings = _warning_ids(step)
    if not warnings:
        return ()
    existing = getattr(step, "condition_selectivity_assessment", None)
    if existing is not None and existing.actionable:
        return ()
    evidence = getattr(step, "condition_evidence", None)
    condition_status = str(getattr(evidence, "status", "unavailable"))
    recommendations = tuple(
        value
        for value in (getattr(evidence, "recommendations", ()) or ())
        if isinstance(value, Mapping)
    )
    if not recommendations:
        return (
            _empty_assessment(
                step,
                status="no_condition_evidence",
                reason="No deterministic condition recommendation is attached.",
                condition_status=condition_status,
            ),
        )
    policy = load_condition_selectivity_repair_policy()
    limited = recommendations[: int(policy["maximum_recommendations_assessed"])]
    training = _training_choice_sets(limited)
    if not training:
        return (
            _empty_assessment(
                step,
                status="unsupported_topology",
                reason=(
                    "Attached precedents do not yield compatible endpoint choice sets."
                ),
                condition_status=condition_status,
            ),
        )
    model = ConditionalEditChoiceModel(
        feature_dimension=int(policy["model_feature_dimension"])
    )
    model.fit(
        training,
        epochs=int(policy["model_epochs"]),
        learning_rate=float(policy["model_learning_rate"]),
    )
    reaction_smiles = str(
        getattr(step.candidate, "condition_query_reaction_smiles", "")
        or getattr(step.candidate, "proposed_reaction_smiles", "")
    )
    assessments: list[ConditionSelectivityRepairAssessment] = []
    for recommendation in limited:
        recipe = recommendation.get("resolved_recipe")
        if not isinstance(recipe, Mapping):
            continue
        recipe_id = str(recommendation.get("recipe_id") or "")
        recipe_core_id = str(recommendation.get("recipe_core_id") or "")
        try:
            query = build_reaction_choice_set(reaction_smiles, recipe)
        except ValueError:
            continue
        model_assessment = model.assess(query)
        reference_ids = _direct_reference_ids(query, training)
        compatibility = (
            float(recommendation["compatibility_score"])
            if recommendation.get("compatibility_score") is not None
            else None
        )
        ranked_first_is_desired = bool(
            model_assessment.ranked_outcomes
            and model_assessment.ranked_outcomes[0].is_selected
        )
        direct_required = bool(policy["require_direct_condition_retrieval"])
        if not recipe_id or not recipe_core_id:
            status: ConditionSelectivityRepairStatus = "insufficient_precedent"
            reason = "The recommendation lacks canonical recipe identity."
        elif compatibility is None or compatibility < float(
            policy["minimum_compatibility_score"]
        ):
            status = "condition_conflict"
            reason = "The recommended recipe does not pass the repair compatibility threshold."
        elif direct_required and condition_status != "recommended_direct":
            status = "insufficient_precedent"
            reason = "Only fallback condition retrieval supports this recipe."
        elif not ranked_first_is_desired:
            status = "competing_outcome_favored"
            reason = "The condition-aware model favors a competing endpoint."
        elif len(reference_ids) < int(policy["minimum_independent_references"]):
            status = "insufficient_precedent"
            reason = "Too few independent exact-condition endpoint precedents are available."
        elif (
            model_assessment.desired_probability
            < float(policy["minimum_desired_probability"])
            or model_assessment.probability_margin
            < float(policy["minimum_probability_margin"])
            or model_assessment.normalized_entropy
            > float(policy["maximum_normalized_entropy"])
        ):
            status = "ambiguous"
            reason = "The intended endpoint does not clear the calibrated confidence gates."
        else:
            status = "supported"
            reason = "Independent exact-condition precedents support the intended endpoint."
        identity = {
            "step_id": step.step_id,
            "warning_ids": warnings,
            "recipe_id": recipe_id,
            "recipe_core_id": recipe_core_id,
            "choice_set_id": model_assessment.choice_set_id,
            "status": status,
            "reference_ids": reference_ids,
            "policy": policy["definition_id"],
        }
        assessments.append(
            ConditionSelectivityRepairAssessment(
                assessment_id=_stable_id("CSRA1", identity),
                step_id=step.step_id,
                warning_ids=warnings,
                recipe_id=recipe_id,
                recipe_core_id=recipe_core_id,
                recommendation_rank=(
                    int(recommendation["rank"])
                    if recommendation.get("rank") is not None
                    else None
                ),
                status=status,
                actionable=status == "supported",
                reason=reason,
                condition_evidence_status=condition_status,
                compatibility_score=compatibility,
                desired_candidate_id=model_assessment.desired_candidate_id,
                desired_probability=model_assessment.desired_probability,
                best_competitor_probability=(
                    model_assessment.best_competitor_probability
                ),
                probability_margin=model_assessment.probability_margin,
                normalized_entropy=model_assessment.normalized_entropy,
                exact_condition_reference_ids=reference_ids,
                training_choice_set_count=len(training),
                ranked_outcomes=_ranked_outcomes(model_assessment),
                policy_definition_id=str(policy["definition_id"]),
                policy_definition_version=str(policy["definition_version"]),
                model_id=model_assessment.model_id,
            )
        )
    if assessments:
        return tuple(
            sorted(
                assessments,
                key=lambda value: (
                    not value.actionable,
                    -(value.probability_margin or -1.0),
                    value.recommendation_rank or 10**9,
                    value.recipe_id,
                ),
            )
        )
    return (
        _empty_assessment(
            step,
            status="unsupported_topology",
            reason="The route reaction cannot be assessed under the attached recipes.",
            condition_status=condition_status,
        ),
    )


__all__ = [
    "CONDITION_SELECTIVITY_REPAIR_SCHEMA_VERSION",
    "ConditionSelectivityRepairAssessment",
    "ConditionSelectivityRepairStatus",
    "assess_condition_selectivity_repairs",
    "load_condition_selectivity_repair_policy",
]
