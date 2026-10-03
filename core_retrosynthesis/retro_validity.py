"""Precedent-backed validity evidence for one concrete retrosynthetic step.

Structural gates, analogue support, recipe compatibility and forward challenges
remain independent evidence. No result is a calibrated success probability and
this tool does not change canonical route admission or search ranking.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from condition_recommender.generic_indexing import GenericReactionIndex
from condition_recommender.reaction_precedent_support import (
    ReactionPrecedentSupport,
    assess_reaction_precedent_support,
    qualify_reaction_precedent,
    reaction_evidence_core,
    summarize_reaction_support,
)
from condition_recommender.recipe_assessment import ReactionRecipeAssessment, assess_reaction_recipe
from condition_recommender.shared_core_index import SharedCoreIndex
from forward_synthesis import RouteStepForwardAssessment

from .chemistry import canonical_smiles, digest
from .external_proposal_assessment import (
    ExternalRetrosynthesisAssessment,
    ExternalRetrosynthesisProposal,
    assess_external_retrosynthesis_proposal,
    load_external_proposal_admission_policy,
)
from .generic_models import GenericTemplateLibrary


POLICY_PATH = Path(__file__).with_name("definitions") / "retro_validity.v1.json"


@dataclass(frozen=True)
class RetroValidityPolicy:
    """Versioned evidence grades and advisory agent decisions."""

    definition_id: str
    schema_version: str
    evidence_levels: tuple[tuple[str, str, int], ...]
    structural_gates: tuple[str, ...]
    suggested_actions: tuple[tuple[str, str], ...]
    candidate_limit: int
    match_limit: int


def validate_retro_validity_policy(value: Mapping[str, Any]) -> None:
    """Reject malformed definitions or promotion to uncalibrated admission."""
    expected_levels = (
        "whole_reaction", "observed_local", "retained_local", "retained_typed", "operator_only",
    )
    if (
        value.get("definition_id") != "retro_validity.v1"
        or value.get("schema_version") != "1.0"
        or value.get("status") != "development_pending_independent_review"
        or value.get("ranking_influence") != "none_advisory_only"
        or value.get("score_semantics") != "ordinal_evidence_not_success_probability"
    ):
        raise ValueError("invalid retro validity definition or evidence boundary")
    levels = value.get("evidence_levels")
    if not isinstance(levels, list) or len(levels) != 5:
        raise ValueError("retro validity evidence levels must be complete")
    grades = []
    for item, level, rank in zip(levels, expected_levels, (4, 3, 2, 1, 0)):
        if (not isinstance(item, Mapping) or item.get("level") != level
                or type(item.get("rank")) is not int or item["rank"] != rank
                or not isinstance(item.get("grade"), str) or not item["grade"]):
            raise ValueError("invalid retro validity evidence grade")
        grades.append(item["grade"])
    if len(set(grades)) != 5:
        raise ValueError("retro validity grades must be distinct")
    if value.get("structural_gates") != [
        "input_structure", "reaction_side_consistency", "atom_correspondence",
        "reaction_completeness", "reaction_signature",
    ]:
        raise ValueError("retro validity structural gates must not be weakened")
    actions = value.get("suggested_actions")
    if not isinstance(actions, Mapping) or set(actions) != {
        "supported", "supported_with_cautions", "insufficient_evidence", "contradicted", "invalid_input",
    } or any(not isinstance(item, str) or not item for item in actions.values()):
        raise ValueError("invalid retro validity suggested actions")
    limits = value.get("limits")
    if not isinstance(limits, Mapping):
        raise ValueError("retro validity requires bounded work limits")
    for name, maximum in (("candidate_limit", 512), ("match_limit", 50)):
        if type(limits.get(name)) is not int or not 1 <= limits[name] <= maximum:
            raise ValueError(f"invalid retro validity {name}")


def load_retro_validity_policy() -> RetroValidityPolicy:
    """Load validated definitions without sharing mutable configuration."""
    value = json.loads(POLICY_PATH.read_text("utf-8"))
    validate_retro_validity_policy(value)
    return RetroValidityPolicy(
        definition_id=value["definition_id"], schema_version=value["schema_version"],
        evidence_levels=tuple((item["level"], item["grade"], item["rank"])
                              for item in value["evidence_levels"]),
        structural_gates=tuple(value["structural_gates"]),
        suggested_actions=tuple(sorted(value["suggested_actions"].items())),
        **value["limits"],
    )


@dataclass(frozen=True)
class RetroForwardEvidence:
    """Separate forward audit bound to the same structures and supplied recipe."""

    starting_materials: str
    intended_product: str
    validity: str
    targeted_replay_status: str
    intended_match: str
    best_competitor_product: str | None
    warnings: tuple[str, ...]
    audited_recipe: Mapping[str, Any] | None = None
    checks: tuple[Mapping[str, Any], ...] = ()
    intended_product_rank: int | None = None
    score_margin: float | None = None
    ranking_definition_id: str | None = None
    audit_schema_version: str | None = None
    schema_version: str = "retro_forward_evidence.v1"

    def __post_init__(self) -> None:
        if self.validity not in {
            "structurally_supported", "structurally_supported_with_competition",
            "inconclusive", "contradicted", "out_of_scope",
        }:
            raise ValueError("invalid forward validity evidence")

    @classmethod
    def from_assessment(
        cls, assessment: RouteStepForwardAssessment, *,
        audited_recipe: Mapping[str, Any] | None = None,
    ) -> RetroForwardEvidence:
        """Project a canonical audit without counting its ranking as probability."""
        return cls(
            starting_materials=assessment.starting_materials,
            intended_product=assessment.intended_product, validity=assessment.validity,
            targeted_replay_status=assessment.targeted_replay_status,
            intended_match=assessment.intended_match,
            best_competitor_product=assessment.best_competitor_product,
            warnings=assessment.warnings, audited_recipe=audited_recipe,
            checks=tuple(check.to_dict() for check in assessment.checks),
            intended_product_rank=assessment.intended_product_rank,
            score_margin=assessment.score_margin,
            ranking_definition_id=assessment.blind_prediction.ranking_definition_id,
            audit_schema_version=assessment.schema_version,
        )


@dataclass(frozen=True)
class RetroValidityAssessment:
    """An advisory decision with independently inspectable evidence axes."""

    assessment_id: str
    status: str
    structural_status: str
    precedent_grade: str
    evidence_rank: int | None
    suggested_action: str
    precedent_query_reaction_smiles: str | None
    precedent_query_evidence: str
    corpus_query_reaction_smiles: str | None
    structural_assessment: ExternalRetrosynthesisAssessment
    operator_precedent_support: ReactionPrecedentSupport
    corpus_precedent_support: ReactionPrecedentSupport | None
    recipe_assessment: ReactionRecipeAssessment | None
    forward_evidence: RetroForwardEvidence | None
    forward_execution_status: str
    forward_status: str
    cautions: tuple[str, ...]
    unresolved_checks: tuple[str, ...]
    warnings: tuple[str, ...]
    definition_id: str
    definition_version: str
    experimental_feasibility: str = "not_established"
    ranking_influence: str = "none_advisory_only"
    score_semantics: str = "ordinal_evidence_not_success_probability"
    schema_version: str = "retro_validity_assessment.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize all evidence, including contradictions and bounded coverage."""
        return asdict(self)


def assess_retro_validity(
    proposal: ExternalRetrosynthesisProposal, operator_library: GenericTemplateLibrary,
    *, condition_index: GenericReactionIndex | None = None,
    shared_core_index: SharedCoreIndex | None = None,
    forward_evidence: RetroForwardEvidence | None = None,
    forward_execution_status: str = "not_run",
    candidate_limit: int | None = None, match_limit: int | None = None,
) -> RetroValidityAssessment:
    """Assess a concrete realization; missing evidence never means impossibility.

    Corpus lookup is optional and source-bound. The forward challenge is supplied
    separately so agents can reserve expensive prediction for consequential gaps.
    This function never expands routes or modifies their admission/ranking.
    """
    policy = load_retro_validity_policy()
    candidates = policy.candidate_limit if candidate_limit is None else candidate_limit
    matches_limit = policy.match_limit if match_limit is None else match_limit
    for name, limit, maximum in (("candidate_limit", candidates, 512), ("match_limit", matches_limit, 50)):
        if type(limit) is not int or not 1 <= limit <= maximum:
            raise ValueError(f"{name} must be an integer between 1 and {maximum}")
    if (condition_index is None) != (shared_core_index is None):
        raise ValueError("condition_index and shared_core_index must be supplied together")
    if forward_execution_status not in {"not_run", "completed", "error", "timed_out", "cancelled"}:
        raise ValueError("invalid forward execution status")
    if forward_evidence is not None:
        if forward_execution_status not in {"not_run", "completed"}:
            raise ValueError("Incomplete forward execution cannot supply valid evidence")
        forward_execution_status = "completed"
    elif forward_execution_status == "completed":
        raise ValueError("Completed forward execution requires assessment evidence")
    assessment = assess_external_retrosynthesis_proposal(proposal, operator_library)
    if forward_evidence is not None:
        if (
            canonical_smiles(forward_evidence.starting_materials) != assessment.canonical_precursor_smiles
            or canonical_smiles(forward_evidence.intended_product) != assessment.canonical_target_smiles
            or not assessment.canonical_target_smiles or not assessment.canonical_precursor_smiles
            or forward_evidence.audited_recipe != proposal.proposed_conditions
        ):
            raise ValueError("Forward evidence must match this realization and supplied recipe exactly")
    gates = {gate.gate_id: gate for gate in assessment.gates}
    if gates["input_structure"].status == "fail":
        structural = "invalid_input"
    elif any(gates[key].status == "fail" for key in policy.structural_gates):
        structural = "contradicted"
    elif all(gates[key].status == "pass" for key in policy.structural_gates):
        structural = "verified"
    else:
        structural = "unresolved"
    # Internal mapping materialization can leave departing atoms unmapped. Keep
    # original structure-backed observations for the same query used by condition
    # retrieval; supplied mappings remain authoritative. Disclose any difference
    # from the assessor's materialized projection rather than silently replacing it.
    evidence_reaction = assessment.normalized_mapped_reaction_smiles if structural == "verified" else None
    query = reaction_evidence_core(evidence_reaction) if evidence_reaction else None
    corpus_reaction = (
        evidence_reaction if proposal.mapped_reaction_smiles else
        f"{assessment.canonical_precursor_smiles}>>{assessment.canonical_target_smiles}"
    ) if structural == "verified" and condition_index is not None else None
    corpus_query = reaction_evidence_core(corpus_reaction) if corpus_reaction else None
    query_disagreement = bool(query and corpus_query and query.levels != corpus_query.levels)
    evidence = []
    skipped = 0
    for precedent in assessment.precedent_matches:
        projection = reaction_evidence_core(precedent.mapped_reaction_smiles) if query else None
        match = qualify_reaction_precedent(
            query, projection, reaction_id=precedent.reaction_id,
            reference_id=precedent.reference_id,
            reaction_smiles=precedent.mapped_reaction_smiles,
            source="operator_library", template_id=precedent.template_id,
        ) if query and projection else None
        if match:
            evidence.append(match)
        else:
            skipped += 1
    operator_support = summarize_reaction_support(
        tuple(evidence), match_limit=matches_limit,
        candidate_count=len(assessment.precedent_matches),
        candidate_truncated=(len(assessment.precedent_matches)
                             >= load_external_proposal_admission_policy().limits.maximum_precedent_matches),
        exclusions=(("NO_QUALIFIED_SHARED_CORE", skipped),) if skipped else (),
        warnings=("COUNTS_DESCRIBE_SAVED_OPERATOR_MATCHES_NOT_FULL_LIBRARY",),
        definition_hash=query.definition_hash if query else None,
        scope="canonical_external_assessment_saved_matches",
    )
    corpus = None
    if condition_index is not None and shared_core_index is not None and corpus_query is not None:
        corpus = assess_reaction_precedent_support(
            corpus_reaction, condition_index, shared_core_index,
            candidate_limit=candidates, match_limit=matches_limit,
        )
    levels = {level: (grade, rank) for level, grade, rank in policy.evidence_levels}
    supported = [item.strongest_level for item in (operator_support, corpus) if item and item.strongest_level]
    best = max(supported, key=lambda level: levels[level][1]) if supported else None
    if best:
        grade, rank = levels[best]
    elif structural == "verified" and any(
        match.match_level == "exact_operator_signature" for match in assessment.operator_matches
    ):
        grade, rank = levels["operator_only"]
    else:
        grade, rank = "unresolved", None
    recipe = (assess_reaction_recipe(
        assessment.normalized_mapped_reaction_smiles
        or f"{proposal.precursor_smiles}>>{proposal.target_smiles}", proposal.proposed_conditions,
    ) if proposal.proposed_conditions is not None else None)
    cautions = set()
    if query_disagreement:
        cautions.add("STRUCTURAL_MAPPING_AND_CORPUS_QUERY_CORE_DIFFER_REVIEW_REQUIRED")
    unresolved = set()
    for key in policy.structural_gates:
        if gates[key].status in {"unresolved", "not_run", "out_of_scope"}:
            unresolved.add(key)
    compatibility = assessment.compatibility_evidence
    hard_conflict = False
    if compatibility:
        for name, item in (("precursor_compatibility", compatibility.precursor),
                           ("reaction_compatibility", compatibility.reaction)):
            if item.disposition != "pass":
                cautions.add(f"{name}:{item.disposition}")
            hard_conflict |= item.disposition == "reject"
        if compatibility.functional_group_competition:
            cautions.add(compatibility.functional_group_competition.code)
    else:
        unresolved.add("compatibility")
    for support in (operator_support, corpus):
        if support:
            for match in support.matches:
                cautions.update(match.differences)
    if recipe is None:
        unresolved.add("supplied_recipe_not_assessed")
    else:
        hard_conflict |= recipe.status == "conflict" or bool(recipe.hard_conflicts)
        cautions.update(recipe.penalty_ids)
        cautions.update(recipe.hard_conflicts)
        unresolved.update(recipe.unresolved_requirements)
        if recipe.status in {"unknown", "invalid_input"}:
            unresolved.add(f"recipe:{recipe.status}")
    if condition_index is None:
        unresolved.add("corpus_index_not_supplied")
    elif corpus_query is None:
        unresolved.add("corpus_lookup_requires_verified_graph")
    elif corpus and corpus.status == "unresolved":
        unresolved.add("corpus_projection_unresolved")
    forward_status = forward_evidence.validity if forward_evidence else forward_execution_status
    if forward_status in {"not_run", "inconclusive", "out_of_scope", "error", "timed_out", "cancelled"}:
        unresolved.add(f"forward:{forward_status}")
    if forward_status == "structurally_supported_with_competition":
        cautions.add("COMPETING_FORWARD_PRODUCTS")
    if forward_evidence and forward_evidence.intended_match in {"stereo_relaxed", "connectivity_only"}:
        cautions.add("FORWARD_MATCH_DOES_NOT_ESTABLISH_EXACT_STEREOCHEMISTRY")
    if structural == "invalid_input":
        status = "invalid_input"
    elif structural == "contradicted" or hard_conflict or forward_status == "contradicted":
        status = "contradicted"
    elif structural != "verified" or rank is None or rank == 0:
        status = "insufficient_evidence"
    elif cautions:
        status = "supported_with_cautions"
    else:
        status = "supported"
    warnings = set(assessment.warnings) | {
        "EVIDENCE_RANK_IS_NOT_EXPERIMENTAL_SUCCESS_PROBABILITY",
        "COMPATIBILITY_RULE_COVERAGE_IS_INCOMPLETE",
        "FORWARD_AND_PRECEDENT_EVIDENCE_MAY_SHARE_SOURCE_DATA",
        "REPORTED_OUTCOMES_AND_RECIPES_REQUIRE_SOURCE_INSPECTION",
    }
    for support in (operator_support, corpus):
        if support:
            warnings.update(support.warnings)
    if forward_evidence:
        warnings.update(forward_evidence.warnings)
    payload = dict(
        status=status, structural_status=structural, precedent_grade=grade, evidence_rank=rank,
        suggested_action=dict(policy.suggested_actions)[status],
        precedent_query_reaction_smiles=evidence_reaction,
        precedent_query_evidence=("validated_supplied_mapping" if proposal.mapped_reaction_smiles
                                  else "canonical_assessor_materialized_mapping"),
        corpus_query_reaction_smiles=corpus_reaction,
        structural_assessment=assessment, operator_precedent_support=operator_support,
        corpus_precedent_support=corpus, recipe_assessment=recipe,
        forward_evidence=forward_evidence, forward_execution_status=forward_execution_status,
        forward_status=forward_status,
        cautions=tuple(sorted(cautions)), unresolved_checks=tuple(sorted(unresolved)),
        warnings=tuple(sorted(warnings)), definition_id=policy.definition_id,
        definition_version=policy.schema_version,
    )
    identity = json.dumps(asdict(RetroValidityAssessment(assessment_id="", **payload)),
                          sort_keys=True, separators=(",", ":"))
    return RetroValidityAssessment(assessment_id=digest("RETROVAL1", identity), **payload)
