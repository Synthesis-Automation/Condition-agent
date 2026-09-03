"""Review-only precedent, condition, and feasibility assessment for partitions.

The assessment consumes the durable Phase 3 route-tree contract. It does not
change route generation, realization status, route cost, or production ranking.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Mapping

from .chemistry import canonical_smiles, digest, split_reaction_smiles
from .condition_ranking import RetrosynthesisConditionEvidence
from .generic_models import GenericTemplateLibrary
from .multistep import ConditionEvidenceEvaluator
from .partition_realization import (
    InterfaceRealization,
    PartitionRealizationResult,
    StrategicRouteRealization,
)
from .precursor_compatibility import assess_precursor_compatibility
from .reaction_compatibility import assess_candidate_reaction_compatibility
from .route_contract import (
    MoleculeOccurrenceNode,
    RouteReactionNode,
    validate_route_tree,
)
from .step_precedents import (
    StepPrecedentLookupResult,
    lookup_reaction_precedents,
)


PARTITION_ASSESSMENT_SCHEMA_VERSION = "1.0"
PARTITION_ASSESSMENT_ALGORITHM_VERSION = "partition_evidence_review.v1"
PARTITION_ASSESSMENT_POLICY_PATH = (
    Path(__file__).with_name("definitions")
    / "synthetic_partition_assessment_policy.v1.json"
)
_EVIDENCE_ORDER = ("E4", "E3", "E2", "E1", "E0")
_ASSESSMENT_STATUSES = (
    "supported",
    "supported_with_cautions",
    "insufficient_evidence",
    "hard_incompatible",
)


@dataclass(frozen=True)
class PartitionAssessmentPolicy:
    """Validated thresholds for review-only route evidence assessment."""

    definition_id: str
    schema_version: str
    exact_product_minimum: float
    exact_precursor_minimum: float
    close_product_minimum: float
    close_precursor_minimum: float
    hard_incompatible_dispositions: tuple[str, ...]
    caution_dispositions: tuple[str, ...]
    latent_atom_classifications: tuple[str, ...]
    review_status_order: tuple[str, ...]
    ranking_influence: str


@dataclass(frozen=True)
class PartitionStepAssessment:
    """Evidence and deterministic cautions for one persisted route reaction."""

    assessment_id: str
    step_id: str
    depth: int
    reaction_smiles: str
    operator_id: str
    template_id: str
    structural_status: str
    structural_issues: tuple[str, ...]
    precedent_evidence_level: str
    precedent_lookup: StepPrecedentLookupResult | None
    independent_precedent_count: int
    condition_evidence: RetrosynthesisConditionEvidence
    precursor_compatibility: Mapping[str, Any]
    reaction_compatibility: Mapping[str, Any]
    hard_incompatible: bool
    cautions: tuple[str, ...]
    warnings: tuple[str, ...]
    schema_version: str = PARTITION_ASSESSMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.structural_status not in {"validated", "invalid"}:
            raise ValueError("unsupported step structural status")
        if self.precedent_evidence_level not in _EVIDENCE_ORDER:
            raise ValueError("unsupported step precedent evidence level")

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible step evidence."""

        return {
            "assessment_id": self.assessment_id,
            "step_id": self.step_id,
            "depth": self.depth,
            "reaction_smiles": self.reaction_smiles,
            "operator_id": self.operator_id,
            "template_id": self.template_id,
            "structural_status": self.structural_status,
            "structural_issues": list(self.structural_issues),
            "precedent_evidence_level": self.precedent_evidence_level,
            "precedent_lookup": (
                self.precedent_lookup.to_dict()
                if self.precedent_lookup is not None
                else None
            ),
            "independent_precedent_count": self.independent_precedent_count,
            "condition_evidence": self.condition_evidence.to_dict(),
            "precursor_compatibility": dict(self.precursor_compatibility),
            "reaction_compatibility": dict(self.reaction_compatibility),
            "hard_incompatible": self.hard_incompatible,
            "cautions": list(self.cautions),
            "warnings": list(self.warnings),
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionStepAssessment":
        """Reconstruct serialized step evidence."""

        precedent = value.get("precedent_lookup")
        return cls(
            assessment_id=str(value["assessment_id"]),
            step_id=str(value["step_id"]),
            depth=int(value["depth"]),
            reaction_smiles=str(value["reaction_smiles"]),
            operator_id=str(value.get("operator_id") or ""),
            template_id=str(value.get("template_id") or ""),
            structural_status=str(value["structural_status"]),
            structural_issues=tuple(
                str(item) for item in value.get("structural_issues") or ()
            ),
            precedent_evidence_level=str(value["precedent_evidence_level"]),
            precedent_lookup=(
                StepPrecedentLookupResult.from_dict(dict(precedent))
                if isinstance(precedent, Mapping)
                else None
            ),
            independent_precedent_count=int(value["independent_precedent_count"]),
            condition_evidence=RetrosynthesisConditionEvidence.from_dict(
                dict(value["condition_evidence"])
            ),
            precursor_compatibility=dict(value.get("precursor_compatibility") or {}),
            reaction_compatibility=dict(value.get("reaction_compatibility") or {}),
            hard_incompatible=bool(value["hard_incompatible"]),
            cautions=tuple(str(item) for item in value.get("cautions") or ()),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            schema_version=str(
                value.get("schema_version") or PARTITION_ASSESSMENT_SCHEMA_VERSION
            ),
        )


@dataclass(frozen=True)
class PartitionInterfaceAssessment:
    """Aggregated step evidence for one requested strategic interface."""

    interface_id: str
    realization_status: str
    route_step_ids: tuple[str, ...]
    assessed_step_ids: tuple[str, ...]
    precedent_match_ids: tuple[str, ...]
    weakest_precedent_evidence: str
    condition_statuses: tuple[str, ...]
    hard_incompatible: bool
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible interface evidence."""

        value = asdict(self)
        for field in (
            "route_step_ids",
            "assessed_step_ids",
            "precedent_match_ids",
            "condition_statuses",
            "warnings",
        ):
            value[field] = list(value[field])
        return value

    @classmethod
    def from_dict(
        cls,
        value: Mapping[str, Any],
    ) -> "PartitionInterfaceAssessment":
        """Reconstruct serialized interface evidence."""

        return cls(
            interface_id=str(value["interface_id"]),
            realization_status=str(value["realization_status"]),
            route_step_ids=tuple(
                str(item) for item in value.get("route_step_ids") or ()
            ),
            assessed_step_ids=tuple(
                str(item) for item in value.get("assessed_step_ids") or ()
            ),
            precedent_match_ids=tuple(
                str(item) for item in value.get("precedent_match_ids") or ()
            ),
            weakest_precedent_evidence=str(value["weakest_precedent_evidence"]),
            condition_statuses=tuple(
                str(item) for item in value.get("condition_statuses") or ()
            ),
            hard_incompatible=bool(value["hard_incompatible"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class PartitionRouteAssessment:
    """Review-only whole-route assessment without feasibility overclaiming."""

    assessment_id: str
    realization_id: str
    source_realization_status: str
    status: str
    step_assessments: tuple[PartitionStepAssessment, ...]
    interface_assessments: tuple[PartitionInterfaceAssessment, ...]
    weakest_step_id: str | None
    weakest_interface_id: str | None
    structurally_valid_step_count: int
    precedent_supported_step_count: int
    condition_supported_step_count: int
    hard_incompatible_step_count: int
    caution_step_count: int
    latent_atom_class_counts: tuple[tuple[str, int], ...]
    protection_burden_count: int
    unresolved_latent_atom_count: int
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.status not in _ASSESSMENT_STATUSES:
            raise ValueError("unsupported partition route assessment status")

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-compatible whole-route evidence."""

        return {
            "assessment_id": self.assessment_id,
            "realization_id": self.realization_id,
            "source_realization_status": self.source_realization_status,
            "status": self.status,
            "step_assessments": [item.to_dict() for item in self.step_assessments],
            "interface_assessments": [
                item.to_dict() for item in self.interface_assessments
            ],
            "weakest_step_id": self.weakest_step_id,
            "weakest_interface_id": self.weakest_interface_id,
            "structurally_valid_step_count": self.structurally_valid_step_count,
            "precedent_supported_step_count": self.precedent_supported_step_count,
            "condition_supported_step_count": self.condition_supported_step_count,
            "hard_incompatible_step_count": self.hard_incompatible_step_count,
            "caution_step_count": self.caution_step_count,
            "latent_atom_class_counts": [
                list(item) for item in self.latent_atom_class_counts
            ],
            "protection_burden_count": self.protection_burden_count,
            "unresolved_latent_atom_count": self.unresolved_latent_atom_count,
            "warnings": list(self.warnings),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionRouteAssessment":
        """Reconstruct a serialized route assessment."""

        return cls(
            assessment_id=str(value["assessment_id"]),
            realization_id=str(value["realization_id"]),
            source_realization_status=str(value["source_realization_status"]),
            status=str(value["status"]),
            step_assessments=tuple(
                PartitionStepAssessment.from_dict(item)
                for item in value.get("step_assessments") or ()
            ),
            interface_assessments=tuple(
                PartitionInterfaceAssessment.from_dict(item)
                for item in value.get("interface_assessments") or ()
            ),
            weakest_step_id=(
                str(value["weakest_step_id"])
                if value.get("weakest_step_id") is not None
                else None
            ),
            weakest_interface_id=(
                str(value["weakest_interface_id"])
                if value.get("weakest_interface_id") is not None
                else None
            ),
            structurally_valid_step_count=int(value["structurally_valid_step_count"]),
            precedent_supported_step_count=int(value["precedent_supported_step_count"]),
            condition_supported_step_count=int(value["condition_supported_step_count"]),
            hard_incompatible_step_count=int(value["hard_incompatible_step_count"]),
            caution_step_count=int(value["caution_step_count"]),
            latent_atom_class_counts=tuple(
                (str(item[0]), int(item[1]))
                for item in value.get("latent_atom_class_counts") or ()
            ),
            protection_burden_count=int(value["protection_burden_count"]),
            unresolved_latent_atom_count=int(value["unresolved_latent_atom_count"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
        )


@dataclass(frozen=True)
class PartitionAssessmentResult:
    """Phase 4 evidence attached to an unchanged Phase 3 result."""

    source: PartitionRealizationResult
    status: str
    route_assessments: tuple[PartitionRouteAssessment, ...]
    policy_definition_id: str
    ranking_influence: str
    warnings: tuple[str, ...]
    algorithm_version: str = PARTITION_ASSESSMENT_ALGORITHM_VERSION
    schema_version: str = PARTITION_ASSESSMENT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.status not in _ASSESSMENT_STATUSES:
            raise ValueError("unsupported partition assessment result status")
        if self.ranking_influence != "none_review_only":
            raise ValueError("partition assessment cannot influence ranking")

    def to_dict(self) -> dict[str, Any]:
        """Return a self-contained JSON-compatible Phase 4 artifact."""

        return {
            "source": self.source.to_dict(),
            "status": self.status,
            "route_assessments": [item.to_dict() for item in self.route_assessments],
            "policy_definition_id": self.policy_definition_id,
            "ranking_influence": self.ranking_influence,
            "warnings": list(self.warnings),
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PartitionAssessmentResult":
        """Reconstruct a serialized Phase 4 assessment."""

        return cls(
            source=PartitionRealizationResult.from_dict(value["source"]),
            status=str(value["status"]),
            route_assessments=tuple(
                PartitionRouteAssessment.from_dict(item)
                for item in value.get("route_assessments") or ()
            ),
            policy_definition_id=str(value["policy_definition_id"]),
            ranking_influence=str(value["ranking_influence"]),
            warnings=tuple(str(item) for item in value.get("warnings") or ()),
            algorithm_version=str(value["algorithm_version"]),
            schema_version=str(value["schema_version"]),
        )


def validate_partition_assessment_policy(value: Mapping[str, Any]) -> list[str]:
    """Return stable validation errors for the Phase 4 policy."""

    errors: list[str] = []
    required = {
        "definition_id",
        "schema_version",
        "precedent_similarity",
        "precedent_evidence_levels",
        "compatibility",
        "latent_atom_classifications",
        "review_status_order",
        "ranking_influence",
    }
    if set(value) != required:
        errors.append("invalid_fields")
    if value.get("definition_id") != "synthetic_partition_assessment_policy.v1":
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
        numbers = tuple(thresholds[field] for field in sorted(threshold_fields))
        if any(
            isinstance(number, bool)
            or not isinstance(number, (int, float))
            or not 0.0 <= float(number) <= 1.0
            for number in numbers
        ):
            errors.append("invalid_precedent_similarity")
        elif float(thresholds["exact_product_minimum"]) < float(
            thresholds["close_product_minimum"]
        ) or float(thresholds["exact_precursor_minimum"]) < float(
            thresholds["close_precursor_minimum"]
        ):
            errors.append("invalid_precedent_similarity")
    levels = value.get("precedent_evidence_levels")
    if not isinstance(levels, Mapping) or dict(levels) != {
        "exact": "E4",
        "close": "E3",
        "same_template": "E2",
        "validated_operator_only": "E1",
        "structurally_invalid": "E0",
    }:
        errors.append("invalid_precedent_evidence_levels")
    compatibility = value.get("compatibility")
    if not isinstance(compatibility, Mapping):
        errors.append("invalid_compatibility")
    else:
        hard = tuple(compatibility.get("hard_incompatible_dispositions") or ())
        caution = tuple(compatibility.get("caution_dispositions") or ())
        allowed = {"pass", "warn", "demote", "reject"}
        if not hard or not set((*hard, *caution)).issubset(allowed):
            errors.append("invalid_compatibility")
    classifications = tuple(value.get("latent_atom_classifications") or ())
    if not classifications or len(classifications) != len(set(classifications)):
        errors.append("invalid_latent_atom_classifications")
    if tuple(value.get("review_status_order") or ()) != _ASSESSMENT_STATUSES:
        errors.append("invalid_review_status_order")
    if value.get("ranking_influence") != "none_review_only":
        errors.append("invalid_ranking_influence")
    return sorted(set(errors))


@lru_cache(maxsize=1)
def load_partition_assessment_policy() -> PartitionAssessmentPolicy:
    """Load and validate the canonical Phase 4 assessment policy."""

    value = json.loads(PARTITION_ASSESSMENT_POLICY_PATH.read_text("utf-8"))
    errors = validate_partition_assessment_policy(value)
    if errors:
        raise ValueError("invalid partition assessment policy: " + ", ".join(errors))
    thresholds = value["precedent_similarity"]
    compatibility = value["compatibility"]
    return PartitionAssessmentPolicy(
        definition_id=str(value["definition_id"]),
        schema_version=str(value["schema_version"]),
        exact_product_minimum=float(thresholds["exact_product_minimum"]),
        exact_precursor_minimum=float(thresholds["exact_precursor_minimum"]),
        close_product_minimum=float(thresholds["close_product_minimum"]),
        close_precursor_minimum=float(thresholds["close_precursor_minimum"]),
        hard_incompatible_dispositions=tuple(
            str(item) for item in compatibility["hard_incompatible_dispositions"]
        ),
        caution_dispositions=tuple(
            str(item) for item in compatibility["caution_dispositions"]
        ),
        latent_atom_classifications=tuple(
            str(item) for item in value["latent_atom_classifications"]
        ),
        review_status_order=tuple(str(item) for item in value["review_status_order"]),
        ranking_influence=str(value["ranking_influence"]),
    )


def _unavailable_condition_evidence(
    reaction_smiles: str,
    warning: str,
) -> RetrosynthesisConditionEvidence:
    return RetrosynthesisConditionEvidence(
        status="insufficient_evidence",
        query_reaction_smiles=reaction_smiles,
        recommender_valid=False,
        recommendation_mode="not_assessed",
        retrieval_level=None,
        uses_type_agnostic_fallback=False,
        candidate_count=0,
        independent_candidate_count=0,
        compatible_candidate_count=0,
        independent_compatible_candidate_count=0,
        excluded_candidate_count=0,
        best_recipe_score=None,
        best_recipe_compatibility_score=None,
        best_recipe_reference_support=0,
        recommendations=(),
        warnings=(warning,),
        error=None,
    )


def _reaction_occurrences(
    root: MoleculeOccurrenceNode,
) -> tuple[tuple[MoleculeOccurrenceNode, RouteReactionNode], ...]:
    values: list[tuple[MoleculeOccurrenceNode, RouteReactionNode]] = []

    def visit(node: MoleculeOccurrenceNode) -> None:
        if node.reaction is None:
            return
        values.append((node, node.reaction))
        for child in node.reaction.children:
            visit(child)

    visit(root)
    return tuple(values)


def _structural_issues(
    owner: MoleculeOccurrenceNode,
    reaction: RouteReactionNode,
    library: GenericTemplateLibrary,
) -> tuple[str, ...]:
    issues: list[str] = []
    split = split_reaction_smiles(reaction.reaction_smiles)
    if split is None:
        return ("invalid_reaction_smiles",)
    precursors, product = split
    expected_precursors = canonical_smiles(
        ".".join(child.smiles for child in reaction.children)
    )
    if canonical_smiles(precursors) != expected_precursors:
        issues.append("reaction_precursors_disagree_with_route_tree")
    if canonical_smiles(product) != canonical_smiles(owner.smiles):
        issues.append("reaction_product_disagrees_with_route_tree")
    if reaction.planned_action is None:
        issues.append("planned_action_missing")
        return tuple(issues)
    template = next(
        (
            item
            for item in library.templates
            if item.template_id == reaction.planned_action.template_id
        ),
        None,
    )
    if template is None:
        issues.append("template_missing_from_library")
    elif template.operator_id != reaction.planned_action.operator_id:
        issues.append("operator_template_mismatch")
    return tuple(issues)


def _precedent_level(
    lookup: StepPrecedentLookupResult | None,
    *,
    structurally_valid: bool,
    policy: PartitionAssessmentPolicy,
) -> str:
    if not structurally_valid:
        return "E0"
    if lookup is None or not lookup.matches:
        return "E1"
    best = lookup.matches[0]
    if (
        best.product_similarity >= policy.exact_product_minimum
        and best.precursor_similarity >= policy.exact_precursor_minimum
    ):
        return "E4"
    if (
        best.product_similarity >= policy.close_product_minimum
        and best.precursor_similarity >= policy.close_precursor_minimum
    ):
        return "E3"
    return "E2"


def _assess_step(
    owner: MoleculeOccurrenceNode,
    reaction: RouteReactionNode,
    library: GenericTemplateLibrary,
    policy: PartitionAssessmentPolicy,
    condition_evaluator: ConditionEvidenceEvaluator | None,
    condition_cache: dict[str, RetrosynthesisConditionEvidence],
    precedent_limit: int,
) -> PartitionStepAssessment:
    structural_issues = _structural_issues(owner, reaction, library)
    structurally_valid = not structural_issues
    action = reaction.planned_action
    operator_id = action.operator_id if action is not None else ""
    template_id = str(action.template_id or "") if action is not None else ""
    split = split_reaction_smiles(reaction.reaction_smiles)
    precursor_smiles = split[0] if split is not None else ""
    precedent_lookup = None
    warnings: set[str] = set()
    if structurally_valid:
        try:
            precedent_lookup = lookup_reaction_precedents(
                step_id=reaction.step_id,
                template_id=template_id,
                operator_id=operator_id,
                product_smiles=owner.smiles,
                precursor_smiles=precursor_smiles,
                library=library,
                limit=precedent_limit,
            )
        except ValueError:
            structural_issues = (*structural_issues, "precedent_lookup_contract_failed")
            structurally_valid = False
    precedent_level = _precedent_level(
        precedent_lookup,
        structurally_valid=structurally_valid,
        policy=policy,
    )
    if precedent_lookup is None or not precedent_lookup.matches:
        warnings.add("NO_ADMITTED_STEP_PRECEDENT")
    if condition_evaluator is None:
        condition = _unavailable_condition_evidence(
            reaction.reaction_smiles,
            "CONDITION_EVIDENCE_NOT_REQUESTED",
        )
    elif not structurally_valid:
        condition = _unavailable_condition_evidence(
            reaction.reaction_smiles,
            "CONDITION_QUERY_BLOCKED_BY_STRUCTURAL_VALIDATION",
        )
    else:
        condition = condition_cache.get(reaction.reaction_smiles)
        if condition is None:
            try:
                condition = condition_evaluator(reaction.reaction_smiles)
            except Exception:
                condition = _unavailable_condition_evidence(
                    reaction.reaction_smiles,
                    "CONDITION_RECOMMENDATION_FAILED",
                )
            condition_cache[reaction.reaction_smiles] = condition
    if structurally_valid:
        precursor_compatibility = assess_precursor_compatibility(
            precursor_smiles
        ).to_dict()
        reaction_compatibility = assess_candidate_reaction_compatibility(
            reaction.reaction_smiles
        ).to_dict()
    else:
        precursor_compatibility = {"disposition": "not_assessed"}
        reaction_compatibility = {"disposition": "not_assessed"}
    dispositions = {
        str(precursor_compatibility.get("disposition") or ""),
        str(reaction_compatibility.get("disposition") or ""),
    }
    hard = not structurally_valid or bool(
        dispositions.intersection(policy.hard_incompatible_dispositions)
    )
    cautions = tuple(
        sorted(
            {
                str(item.get("message") or item.get("rule_id") or "compatibility")
                for result in (precursor_compatibility, reaction_compatibility)
                for item in result.get("assessments") or ()
                if isinstance(item, Mapping)
            }
        )
    )
    if dispositions.intersection(policy.caution_dispositions):
        warnings.add("COMPATIBILITY_REVIEW_REQUIRED")
    if condition.status == "insufficient_evidence":
        warnings.add("CONDITION_EVIDENCE_INCOMPLETE")
    if condition.status == "recommended_fallback":
        warnings.add("CONDITION_FALLBACK_USED")
    if hard:
        warnings.add("HARD_INCOMPATIBILITY")
    independent_count = (
        len({item.reference_id for item in precedent_lookup.matches})
        if precedent_lookup is not None
        else 0
    )
    assessment_id = digest(
        "SPSTEP4",
        reaction.step_id,
        precedent_level,
        condition.status,
        "hard" if hard else "reviewable",
        policy.definition_id,
    )
    return PartitionStepAssessment(
        assessment_id=assessment_id,
        step_id=reaction.step_id,
        depth=reaction.depth,
        reaction_smiles=reaction.reaction_smiles,
        operator_id=operator_id,
        template_id=template_id,
        structural_status="validated" if structurally_valid else "invalid",
        structural_issues=tuple(structural_issues),
        precedent_evidence_level=precedent_level,
        precedent_lookup=precedent_lookup,
        independent_precedent_count=independent_count,
        condition_evidence=condition,
        precursor_compatibility=precursor_compatibility,
        reaction_compatibility=reaction_compatibility,
        hard_incompatible=hard,
        cautions=cautions,
        warnings=tuple(sorted(warnings)),
    )


def _interface_assessment(
    interface: InterfaceRealization,
    by_step: Mapping[str, PartitionStepAssessment],
) -> PartitionInterfaceAssessment:
    steps = tuple(
        by_step[step_id] for step_id in interface.route_step_ids if step_id in by_step
    )
    levels = tuple(item.precedent_evidence_level for item in steps)
    weakest = max(levels, key=_EVIDENCE_ORDER.index) if levels else "E0"
    matches = tuple(
        sorted(
            {
                match.match_id
                for step in steps
                if step.precedent_lookup is not None
                for match in step.precedent_lookup.matches
            }
        )
    )
    warnings = set(interface.warnings)
    if interface.status != "realized":
        warnings.add("INTERFACE_NOT_REALIZED")
    if not steps:
        warnings.add("INTERFACE_HAS_NO_ASSESSED_STEP")
    if not matches:
        warnings.add("INTERFACE_HAS_NO_ADMITTED_PRECEDENT")
    return PartitionInterfaceAssessment(
        interface_id=interface.interface_id,
        realization_status=interface.status,
        route_step_ids=interface.route_step_ids,
        assessed_step_ids=tuple(item.step_id for item in steps),
        precedent_match_ids=matches,
        weakest_precedent_evidence=weakest,
        condition_statuses=tuple(
            sorted({item.condition_evidence.status for item in steps})
        ),
        hard_incompatible=any(item.hard_incompatible for item in steps),
        warnings=tuple(sorted(warnings)),
    )


def _step_weakness_key(step: PartitionStepAssessment) -> tuple[Any, ...]:
    return (
        0 if step.hard_incompatible else 1,
        -_EVIDENCE_ORDER.index(step.precedent_evidence_level),
        0 if step.condition_evidence.status == "insufficient_evidence" else 1,
        -len(step.cautions),
        step.step_id,
    )


def _route_assessment(
    realization: StrategicRouteRealization,
    library: GenericTemplateLibrary,
    policy: PartitionAssessmentPolicy,
    condition_evaluator: ConditionEvidenceEvaluator | None,
    condition_cache: dict[str, RetrosynthesisConditionEvidence],
    precedent_limit: int,
) -> PartitionRouteAssessment:
    tree_report = validate_route_tree(realization.route_tree)
    occurrences = _reaction_occurrences(realization.route_tree.root)
    steps = tuple(
        _assess_step(
            owner,
            reaction,
            library,
            policy,
            condition_evaluator,
            condition_cache,
            precedent_limit,
        )
        for owner, reaction in occurrences
    )
    by_step = {item.step_id: item for item in steps}
    interfaces = tuple(
        _interface_assessment(item, by_step)
        for item in realization.interface_realizations
    )
    class_counts = {key: 0 for key in policy.latent_atom_classifications}
    for state in realization.frontier_states:
        for atom in state.non_target_atoms:
            class_counts.setdefault(atom.classification, 0)
            class_counts[atom.classification] += 1
    hard_count = sum(item.hard_incompatible for item in steps)
    condition_supported = sum(
        item.condition_evidence.status in {"recommended_direct", "recommended_fallback"}
        for item in steps
    )
    precedent_supported = sum(
        item.precedent_evidence_level in {"E4", "E3", "E2"} for item in steps
    )
    caution_count = sum(bool(item.cautions or item.warnings) for item in steps)
    warnings = set(realization.warnings)
    if not tree_report.valid:
        warnings.add("ROUTE_TREE_INVALID")
    if hard_count or any(item.hard_incompatible for item in interfaces):
        status = "hard_incompatible"
        warnings.add("HARD_INCOMPATIBLE_ROUTE_NOT_PROMOTABLE")
    elif not steps or precedent_supported < len(steps):
        status = "insufficient_evidence"
        warnings.add("STEP_PRECEDENT_COVERAGE_INCOMPLETE")
    elif (
        condition_supported < len(steps)
        or caution_count
        or realization.status != "fully_realized"
    ):
        status = "supported_with_cautions"
    else:
        status = "supported"
    weakest_step = min(steps, key=_step_weakness_key) if steps else None
    weakest_interface = (
        max(
            interfaces,
            key=lambda item: (
                item.hard_incompatible,
                _EVIDENCE_ORDER.index(item.weakest_precedent_evidence),
                item.interface_id,
            ),
        )
        if interfaces
        else None
    )
    protection_burden = class_counts.get("PROTECTING_GROUP", 0) + class_counts.get(
        "AUXILIARY", 0
    )
    unresolved_atoms = class_counts.get("UNKNOWN", 0)
    if protection_burden:
        warnings.add("PROTECTION_OR_AUXILIARY_BURDEN_PRESENT")
    if unresolved_atoms:
        warnings.add("LATENT_ATOM_CLASSIFICATION_INCOMPLETE")
    assessment_id = digest(
        "SPROUTE4",
        realization.realization_id,
        status,
        *(item.assessment_id for item in steps),
        policy.definition_id,
    )
    return PartitionRouteAssessment(
        assessment_id=assessment_id,
        realization_id=realization.realization_id,
        source_realization_status=realization.status,
        status=status,
        step_assessments=steps,
        interface_assessments=interfaces,
        weakest_step_id=weakest_step.step_id if weakest_step is not None else None,
        weakest_interface_id=(
            weakest_interface.interface_id if weakest_interface is not None else None
        ),
        structurally_valid_step_count=sum(
            item.structural_status == "validated" for item in steps
        ),
        precedent_supported_step_count=precedent_supported,
        condition_supported_step_count=condition_supported,
        hard_incompatible_step_count=hard_count,
        caution_step_count=caution_count,
        latent_atom_class_counts=tuple(
            (key, class_counts[key]) for key in sorted(class_counts)
        ),
        protection_burden_count=protection_burden,
        unresolved_latent_atom_count=unresolved_atoms,
        warnings=tuple(sorted(warnings)),
    )


def assess_partition_realizations(
    result: PartitionRealizationResult,
    library: GenericTemplateLibrary,
    *,
    condition_evaluator: ConditionEvidenceEvaluator | None = None,
    precedent_limit: int = 5,
) -> PartitionAssessmentResult:
    """Attach Phase 4 evidence without changing Phase 3 ordering or status."""

    if precedent_limit < 1 or precedent_limit > 20:
        raise ValueError("precedent limit must be between one and twenty")
    policy = load_partition_assessment_policy()
    condition_cache: dict[str, RetrosynthesisConditionEvidence] = {}
    assessments = tuple(
        _route_assessment(
            realization,
            library,
            policy,
            condition_evaluator,
            condition_cache,
            precedent_limit,
        )
        for realization in result.realizations
    )
    if not assessments:
        status = "insufficient_evidence"
    else:
        status = max(
            (item.status for item in assessments),
            key=policy.review_status_order.index,
        )
    warnings = {
        "REVIEW_ONLY_DOES_NOT_CHANGE_ROUTE_RANKING",
        "PRECEDENT_AND_CONDITION_SUPPORT_DO_NOT_GUARANTEE_LABORATORY_SUCCESS",
    }
    if condition_evaluator is None:
        warnings.add("CONDITION_EVIDENCE_NOT_REQUESTED")
    if not assessments:
        warnings.add("NO_ROUTE_REALIZATION_TO_ASSESS")
    return PartitionAssessmentResult(
        source=result,
        status=status,
        route_assessments=assessments,
        policy_definition_id=policy.definition_id,
        ranking_influence=policy.ranking_influence,
        warnings=tuple(sorted(warnings)),
    )


__all__ = [
    "PARTITION_ASSESSMENT_ALGORITHM_VERSION",
    "PARTITION_ASSESSMENT_POLICY_PATH",
    "PARTITION_ASSESSMENT_SCHEMA_VERSION",
    "PartitionAssessmentPolicy",
    "PartitionAssessmentResult",
    "PartitionInterfaceAssessment",
    "PartitionRouteAssessment",
    "PartitionStepAssessment",
    "assess_partition_realizations",
    "load_partition_assessment_policy",
    "validate_partition_assessment_policy",
]
