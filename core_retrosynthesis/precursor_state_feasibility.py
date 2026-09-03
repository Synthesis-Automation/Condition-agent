"""Review-only evidence for the chemical support of latent precursor states."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .chemistry import digest
from .latent_transition import LatentStateTransition
from .route_contract import MoleculeOccurrenceNode
from .step_precedents import StepPrecedentLookupResult


PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION = "1.0"
PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION = (
    "precursor_state_feasibility.v1"
)
PRECURSOR_STATE_FEASIBILITY_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "precursor_state_feasibility.v1.json"
)
PRECURSOR_STATE_EVIDENCE_LEVELS = ("E4", "E3", "E2", "E1", "E0")
PRECURSOR_STATE_SUPPORT_STATUSES = (
    "supported",
    "supported_with_cautions",
    "insufficient_evidence",
    "contradicted",
)
PRECURSOR_STATE_PROMOTION_RECOMMENDATIONS = (
    "eligible_for_route_review",
    "hold_for_evidence",
    "not_promotable",
)
PRECURSOR_STATE_REACTANT_SUPPORT = (
    "exact_reactant_state",
    "close_reactant_state",
    "template_only",
    "operator_only",
    "contradicted",
)
PRECURSOR_STATE_DIRECTIONALITY_STATUSES = (
    "reactant_precedent_supported",
    "template_only",
    "not_precedent_supported",
    "contradicted",
)


@dataclass(frozen=True)
class PrecursorStateFeasibilityPolicy:
    """Validated evidence thresholds for latent precursor-state assessment."""

    definition_id: str
    schema_version: str
    exact_product_minimum: float
    exact_precursor_minimum: float
    close_product_minimum: float
    close_precursor_minimum: float
    strong_terminal_evidence: tuple[str, ...]
    weak_terminal_evidence: tuple[str, ...]
    supported_evidence_levels: tuple[str, ...]
    hold_evidence_levels: tuple[str, ...]
    contradicted_evidence_levels: tuple[str, ...]
    ranking_influence: str


@dataclass(frozen=True)
class PrecursorStateFeasibility:
    """Evidence that one proposed set of precursor states is chemically usable."""

    assessment_id: str
    step_id: str
    transition_id: str | None
    input_state_ids: tuple[str, ...]
    output_state_id: str | None
    tactical_input_count: int
    evidence_level: str
    support_status: str
    promotion_recommendation: str
    reactant_state_support: str
    directionality_status: str
    best_product_similarity: float | None
    best_precursor_similarity: float | None
    independent_precedent_count: int
    derived_input_count: int
    strong_terminal_input_count: int
    weak_terminal_input_count: int
    unresolved_input_count: int
    condition_status: str
    compatibility_dispositions: tuple[str, ...]
    reasons: tuple[str, ...]
    warnings: tuple[str, ...]
    policy_definition_id: str
    algorithm_version: str = PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION
    schema_version: str = PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION:
            raise ValueError("unsupported precursor-state schema")
        if self.algorithm_version != PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION:
            raise ValueError("unsupported precursor-state algorithm")
        if self.evidence_level not in PRECURSOR_STATE_EVIDENCE_LEVELS:
            raise ValueError("unsupported precursor-state evidence level")
        if self.support_status not in PRECURSOR_STATE_SUPPORT_STATUSES:
            raise ValueError("unsupported precursor-state support status")
        if (
            self.promotion_recommendation
            not in PRECURSOR_STATE_PROMOTION_RECOMMENDATIONS
        ):
            raise ValueError(
                "unsupported precursor-state promotion recommendation"
            )
        if self.reactant_state_support not in PRECURSOR_STATE_REACTANT_SUPPORT:
            raise ValueError("unsupported reactant-state support")
        if self.directionality_status not in PRECURSOR_STATE_DIRECTIONALITY_STATUSES:
            raise ValueError("unsupported precursor-state directionality status")
        counts = (
            self.independent_precedent_count,
            self.derived_input_count,
            self.strong_terminal_input_count,
            self.weak_terminal_input_count,
            self.unresolved_input_count,
            self.tactical_input_count,
        )
        if any(value < 0 for value in counts):
            raise ValueError("precursor-state evidence counts must be nonnegative")
        similarities = (
            self.best_product_similarity,
            self.best_precursor_similarity,
        )
        if any(
            value is not None and not 0.0 <= value <= 1.0
            for value in similarities
        ):
            raise ValueError("precursor-state similarities must be in [0, 1]")

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible precursor-state evidence."""

        value = asdict(self)
        value["input_state_ids"] = list(self.input_state_ids)
        value["compatibility_dispositions"] = list(
            self.compatibility_dispositions
        )
        value["reasons"] = list(self.reasons)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(
        cls,
        value: Mapping[str, Any],
    ) -> "PrecursorStateFeasibility":
        """Reconstruct serialized precursor-state evidence."""

        return cls(
            assessment_id=str(value["assessment_id"]),
            step_id=str(value["step_id"]),
            transition_id=(
                str(value["transition_id"])
                if value.get("transition_id") is not None
                else None
            ),
            input_state_ids=tuple(
                str(item) for item in value.get("input_state_ids") or ()
            ),
            output_state_id=(
                str(value["output_state_id"])
                if value.get("output_state_id") is not None
                else None
            ),
            tactical_input_count=int(value["tactical_input_count"]),
            evidence_level=str(value["evidence_level"]),
            support_status=str(value["support_status"]),
            promotion_recommendation=str(value["promotion_recommendation"]),
            reactant_state_support=str(value["reactant_state_support"]),
            directionality_status=str(value["directionality_status"]),
            best_product_similarity=(
                float(value["best_product_similarity"])
                if value.get("best_product_similarity") is not None
                else None
            ),
            best_precursor_similarity=(
                float(value["best_precursor_similarity"])
                if value.get("best_precursor_similarity") is not None
                else None
            ),
            independent_precedent_count=int(value["independent_precedent_count"]),
            derived_input_count=int(value["derived_input_count"]),
            strong_terminal_input_count=int(
                value["strong_terminal_input_count"]
            ),
            weak_terminal_input_count=int(value["weak_terminal_input_count"]),
            unresolved_input_count=int(value["unresolved_input_count"]),
            condition_status=str(value["condition_status"]),
            compatibility_dispositions=tuple(
                str(item)
                for item in value.get("compatibility_dispositions") or ()
            ),
            reasons=tuple(str(item) for item in value.get("reasons") or ()),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            policy_definition_id=str(value["policy_definition_id"]),
            algorithm_version=str(
                value.get("algorithm_version")
                or PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION
            ),
            schema_version=str(
                value.get("schema_version")
                or PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION
            ),
        )


@dataclass(frozen=True)
class PrecursorStateRouteFeasibility:
    """Weakest-link aggregation of precursor-state evidence over a route."""

    assessment_id: str
    status: str
    promotion_recommendation: str
    supported_step_count: int
    held_step_count: int
    contradicted_step_count: int
    strong_terminal_input_count: int
    weak_terminal_input_count: int
    unresolved_input_count: int
    weakest_step_id: str | None
    weakest_evidence_level: str
    warnings: tuple[str, ...]
    policy_definition_id: str
    algorithm_version: str = PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION
    schema_version: str = PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION:
            raise ValueError("unsupported route precursor-state schema")
        if self.algorithm_version != PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION:
            raise ValueError("unsupported route precursor-state algorithm")
        if self.status not in PRECURSOR_STATE_SUPPORT_STATUSES:
            raise ValueError("unsupported route precursor-state status")
        if (
            self.promotion_recommendation
            not in PRECURSOR_STATE_PROMOTION_RECOMMENDATIONS
        ):
            raise ValueError(
                "unsupported route precursor-state promotion recommendation"
            )
        if self.weakest_evidence_level not in PRECURSOR_STATE_EVIDENCE_LEVELS:
            raise ValueError("unsupported weakest precursor-state evidence")
        counts = (
            self.supported_step_count,
            self.held_step_count,
            self.contradicted_step_count,
            self.strong_terminal_input_count,
            self.weak_terminal_input_count,
            self.unresolved_input_count,
        )
        if any(value < 0 for value in counts):
            raise ValueError(
                "route precursor-state evidence counts must be nonnegative"
            )

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible route precursor-state evidence."""

        value = asdict(self)
        value["warnings"] = list(self.warnings)
        return value

    @classmethod
    def from_dict(
        cls,
        value: Mapping[str, Any],
    ) -> "PrecursorStateRouteFeasibility":
        """Reconstruct serialized route precursor-state evidence."""

        return cls(
            assessment_id=str(value["assessment_id"]),
            status=str(value["status"]),
            promotion_recommendation=str(value["promotion_recommendation"]),
            supported_step_count=int(value["supported_step_count"]),
            held_step_count=int(value["held_step_count"]),
            contradicted_step_count=int(value["contradicted_step_count"]),
            strong_terminal_input_count=int(
                value["strong_terminal_input_count"]
            ),
            weak_terminal_input_count=int(value["weak_terminal_input_count"]),
            unresolved_input_count=int(value["unresolved_input_count"]),
            weakest_step_id=(
                str(value["weakest_step_id"])
                if value.get("weakest_step_id") is not None
                else None
            ),
            weakest_evidence_level=str(value["weakest_evidence_level"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            policy_definition_id=str(value["policy_definition_id"]),
            algorithm_version=str(
                value.get("algorithm_version")
                or PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION
            ),
            schema_version=str(
                value.get("schema_version")
                or PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION
            ),
        )


def validate_precursor_state_feasibility_policy(
    value: Mapping[str, Any],
) -> list[str]:
    """Return stable validation errors for a feasibility policy."""

    errors: list[str] = []
    required = {
        "definition_id",
        "schema_version",
        "precedent_similarity",
        "strong_terminal_evidence",
        "weak_terminal_evidence",
        "supported_evidence_levels",
        "hold_evidence_levels",
        "contradicted_evidence_levels",
        "ranking_influence",
    }
    if set(value) != required:
        errors.append("invalid_fields")
    if value.get("definition_id") != "precursor_state_feasibility.v1":
        errors.append("invalid_definition_id")
    if value.get("schema_version") != "1.0":
        errors.append("invalid_schema_version")
    thresholds = value.get("precedent_similarity")
    threshold_fields = {
        "exact_product_minimum",
        "exact_precursor_minimum",
        "close_product_minimum",
        "close_precursor_minimum",
    }
    if not isinstance(thresholds, Mapping) or set(thresholds) != threshold_fields:
        errors.append("invalid_precedent_similarity")
    else:
        numbers = tuple(thresholds[field] for field in threshold_fields)
        if any(
            isinstance(number, bool)
            or not isinstance(number, (int, float))
            or not 0.0 <= float(number) <= 1.0
            for number in numbers
        ):
            errors.append("invalid_precedent_similarity")
        elif (
            float(thresholds["exact_product_minimum"])
            < float(thresholds["close_product_minimum"])
            or float(thresholds["exact_precursor_minimum"])
            < float(thresholds["close_precursor_minimum"])
        ):
            errors.append("invalid_precedent_similarity")
    strong = tuple(value.get("strong_terminal_evidence") or ())
    weak = tuple(value.get("weak_terminal_evidence") or ())
    if not strong or not weak or set(strong).intersection(weak):
        errors.append("invalid_terminal_evidence")
    evidence_groups = (
        tuple(value.get("supported_evidence_levels") or ()),
        tuple(value.get("hold_evidence_levels") or ()),
        tuple(value.get("contradicted_evidence_levels") or ()),
    )
    flattened = tuple(item for group in evidence_groups for item in group)
    if (
        set(flattened) != set(PRECURSOR_STATE_EVIDENCE_LEVELS)
        or len(flattened) != len(set(flattened))
    ):
        errors.append("invalid_evidence_levels")
    if value.get("ranking_influence") != "none_review_only":
        errors.append("invalid_ranking_influence")
    return sorted(set(errors))


@lru_cache(maxsize=1)
def load_precursor_state_feasibility_policy() -> PrecursorStateFeasibilityPolicy:
    """Load the canonical precursor-state feasibility policy."""

    value = json.loads(
        PRECURSOR_STATE_FEASIBILITY_POLICY_PATH.read_text(encoding="utf-8")
    )
    errors = validate_precursor_state_feasibility_policy(value)
    if errors:
        raise ValueError(
            "invalid precursor-state feasibility policy: " + ", ".join(errors)
        )
    thresholds = value["precedent_similarity"]
    return PrecursorStateFeasibilityPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        exact_product_minimum=float(thresholds["exact_product_minimum"]),
        exact_precursor_minimum=float(thresholds["exact_precursor_minimum"]),
        close_product_minimum=float(thresholds["close_product_minimum"]),
        close_precursor_minimum=float(thresholds["close_precursor_minimum"]),
        strong_terminal_evidence=tuple(value["strong_terminal_evidence"]),
        weak_terminal_evidence=tuple(value["weak_terminal_evidence"]),
        supported_evidence_levels=tuple(value["supported_evidence_levels"]),
        hold_evidence_levels=tuple(value["hold_evidence_levels"]),
        contradicted_evidence_levels=tuple(
            value["contradicted_evidence_levels"]
        ),
        ranking_influence=str(value["ranking_influence"]),
    )


def _precedent_support(
    lookup: StepPrecedentLookupResult | None,
    *,
    structurally_valid: bool,
    hard_incompatible: bool,
    policy: PrecursorStateFeasibilityPolicy,
) -> tuple[str, str, float | None, float | None, int]:
    if not structurally_valid or hard_incompatible:
        return "E0", "contradicted", None, None, 0
    if lookup is None or not lookup.matches:
        return "E1", "operator_only", None, None, 0
    best = lookup.matches[0]
    independent = len({item.reference_id for item in lookup.matches})
    if (
        best.product_similarity >= policy.exact_product_minimum
        and best.precursor_similarity >= policy.exact_precursor_minimum
    ):
        level = "E4"
        support = "exact_reactant_state"
    elif (
        best.product_similarity >= policy.close_product_minimum
        and best.precursor_similarity >= policy.close_precursor_minimum
    ):
        level = "E3"
        support = "close_reactant_state"
    else:
        level = "E2"
        support = "template_only"
    return (
        level,
        support,
        best.product_similarity,
        best.precursor_similarity,
        independent,
    )


def assess_precursor_state_feasibility(
    *,
    step_id: str,
    transition: LatentStateTransition | None,
    precursor_nodes: Iterable[MoleculeOccurrenceNode],
    precedent_lookup: StepPrecedentLookupResult | None,
    structurally_valid: bool,
    hard_incompatible: bool,
    compatibility_dispositions: Iterable[str],
    condition_status: str,
    policy: PrecursorStateFeasibilityPolicy | None = None,
) -> PrecursorStateFeasibility:
    """Assess reactant-state support without changing route ranking."""

    resolved = policy or load_precursor_state_feasibility_policy()
    level, reactant_support, product_similarity, precursor_similarity, independent = (
        _precedent_support(
            precedent_lookup,
            structurally_valid=structurally_valid,
            hard_incompatible=hard_incompatible,
            policy=resolved,
        )
    )
    derived = 0
    strong = 0
    weak = 0
    unresolved = 0
    for node in precursor_nodes:
        if node.reaction is not None:
            derived += 1
        elif node.terminal_evidence in resolved.strong_terminal_evidence:
            strong += 1
        elif node.terminal_evidence in resolved.weak_terminal_evidence:
            weak += 1
        else:
            unresolved += 1
    reasons = {
        "FORWARD_PRODUCT_RECONSTRUCTION_VALIDATED"
        if structurally_valid
        else "FORWARD_PRODUCT_RECONSTRUCTION_INVALID"
    }
    warnings: set[str] = set()
    if level == "E4":
        reasons.add("EXACT_REACTANT_STATE_PRECEDENT")
    elif level == "E3":
        reasons.add("CLOSE_REACTANT_STATE_PRECEDENT")
    elif level == "E2":
        reasons.add("SAME_TEMPLATE_PRECEDENT_ONLY")
        warnings.add("TEMPLATE_ONLY_PRECURSOR_STATE_SUPPORT")
    elif level == "E1":
        reasons.add("VALIDATED_OPERATOR_WITHOUT_REACTANT_STATE_PRECEDENT")
        warnings.add("NO_REACTANT_STATE_PRECEDENT")
    else:
        reasons.add("STRUCTURAL_OR_COMPATIBILITY_CONTRADICTION")
        warnings.add("PRECURSOR_STATE_CONTRADICTED")
    dispositions = set(compatibility_dispositions)
    if dispositions.intersection({"warn", "demote"}):
        warnings.add("PRECURSOR_OR_REACTION_COMPATIBILITY_CAUTION")
    if transition is None:
        warnings.add("LATENT_STATE_TRANSITION_UNAVAILABLE")
    if weak:
        warnings.add("WEAK_TERMINAL_INPUT_EVIDENCE")
    if unresolved:
        warnings.add("UNRESOLVED_INPUT_STATE")
    if condition_status == "insufficient_evidence":
        warnings.add("CONDITION_SUPPORT_ABSENT")
    if level in resolved.contradicted_evidence_levels:
        support_status = "contradicted"
        promotion = "not_promotable"
        directionality = "contradicted"
    elif level in resolved.hold_evidence_levels or transition is None:
        support_status = "insufficient_evidence"
        promotion = "hold_for_evidence"
        directionality = (
            "template_only" if level == "E2" else "not_precedent_supported"
        )
    else:
        support_status = (
            "supported_with_cautions" if warnings else "supported"
        )
        promotion = "eligible_for_route_review"
        directionality = "reactant_precedent_supported"
    assessment_id = digest(
        "PSTATE5E",
        step_id,
        transition.transition_id if transition is not None else "none",
        level,
        support_status,
        promotion,
        resolved.definition_id,
    )
    return PrecursorStateFeasibility(
        assessment_id=assessment_id,
        step_id=step_id,
        transition_id=(transition.transition_id if transition is not None else None),
        input_state_ids=(
            transition.input_state_ids if transition is not None else ()
        ),
        output_state_id=(
            transition.output_state_id if transition is not None else None
        ),
        tactical_input_count=(
            len(transition.tactical_input_occurrence_ids)
            if transition is not None
            else 0
        ),
        evidence_level=level,
        support_status=support_status,
        promotion_recommendation=promotion,
        reactant_state_support=reactant_support,
        directionality_status=directionality,
        best_product_similarity=product_similarity,
        best_precursor_similarity=precursor_similarity,
        independent_precedent_count=independent,
        derived_input_count=derived,
        strong_terminal_input_count=strong,
        weak_terminal_input_count=weak,
        unresolved_input_count=unresolved,
        condition_status=condition_status,
        compatibility_dispositions=tuple(sorted(dispositions)),
        reasons=tuple(sorted(reasons)),
        warnings=tuple(sorted(warnings)),
        policy_definition_id=resolved.definition_id,
    )


def aggregate_precursor_state_route_feasibility(
    assessments: Iterable[PrecursorStateFeasibility],
    *,
    policy: PrecursorStateFeasibilityPolicy | None = None,
) -> PrecursorStateRouteFeasibility:
    """Aggregate step feasibility using a conservative weakest-link rule."""

    resolved = policy or load_precursor_state_feasibility_policy()
    values = tuple(assessments)
    weakest = (
        max(
            values,
            key=lambda item: (
                PRECURSOR_STATE_EVIDENCE_LEVELS.index(item.evidence_level),
                item.step_id,
            ),
        )
        if values
        else None
    )
    supported = sum(
        item.evidence_level in resolved.supported_evidence_levels
        for item in values
    )
    held = sum(
        item.promotion_recommendation == "hold_for_evidence"
        for item in values
    )
    contradicted = sum(
        item.promotion_recommendation == "not_promotable"
        for item in values
    )
    strong = sum(item.strong_terminal_input_count for item in values)
    weak = sum(item.weak_terminal_input_count for item in values)
    unresolved = sum(item.unresolved_input_count for item in values)
    warnings: set[str] = set()
    if contradicted:
        status = "contradicted"
        promotion = "not_promotable"
        warnings.add("CONTRADICTED_PRECURSOR_STATE_PRESENT")
    elif not values or held or supported < len(values):
        status = "insufficient_evidence"
        promotion = "hold_for_evidence"
        warnings.add("PRECURSOR_STATE_SUPPORT_INCOMPLETE")
    elif weak or unresolved:
        status = "supported_with_cautions"
        promotion = "hold_for_evidence"
    else:
        status = (
            "supported_with_cautions"
            if any(item.warnings for item in values)
            else "supported"
        )
        promotion = "eligible_for_route_review"
    if weak:
        warnings.add("ROUTE_HAS_WEAK_TERMINAL_INPUT_EVIDENCE")
    if unresolved:
        warnings.add("ROUTE_HAS_UNRESOLVED_INPUT_STATES")
    if any(item.condition_status == "insufficient_evidence" for item in values):
        warnings.add("ROUTE_CONDITION_SUPPORT_INCOMPLETE")
    weakest_level = weakest.evidence_level if weakest is not None else "E0"
    assessment_id = digest(
        "PSTATEROUTE5E",
        *(item.assessment_id for item in values),
        status,
        promotion,
        resolved.definition_id,
    )
    return PrecursorStateRouteFeasibility(
        assessment_id=assessment_id,
        status=status,
        promotion_recommendation=promotion,
        supported_step_count=supported,
        held_step_count=held,
        contradicted_step_count=contradicted,
        strong_terminal_input_count=strong,
        weak_terminal_input_count=weak,
        unresolved_input_count=unresolved,
        weakest_step_id=weakest.step_id if weakest is not None else None,
        weakest_evidence_level=weakest_level,
        warnings=tuple(sorted(warnings)),
        policy_definition_id=resolved.definition_id,
    )


__all__ = [
    "PRECURSOR_STATE_EVIDENCE_LEVELS",
    "PRECURSOR_STATE_FEASIBILITY_ALGORITHM_VERSION",
    "PRECURSOR_STATE_FEASIBILITY_POLICY_PATH",
    "PRECURSOR_STATE_FEASIBILITY_SCHEMA_VERSION",
    "PRECURSOR_STATE_PROMOTION_RECOMMENDATIONS",
    "PRECURSOR_STATE_REACTANT_SUPPORT",
    "PRECURSOR_STATE_DIRECTIONALITY_STATUSES",
    "PRECURSOR_STATE_SUPPORT_STATUSES",
    "PrecursorStateFeasibility",
    "PrecursorStateFeasibilityPolicy",
    "PrecursorStateRouteFeasibility",
    "aggregate_precursor_state_route_feasibility",
    "assess_precursor_state_feasibility",
    "load_precursor_state_feasibility_policy",
    "validate_precursor_state_feasibility_policy",
]
