"""Validity axes remain distinct for positive, ambiguous and conflicting steps."""

from dataclasses import replace
import json

import pytest

from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.shared_core_index import build_shared_core_index, load_shared_core_index
from core_retrosynthesis import (
    RetroForwardEvidence, assess_retro_validity, load_retro_validity_policy, validate_retro_validity_policy,
)
from core_retrosynthesis.external_proposal_assessment import ExternalRetrosynthesisProposal
from core_retrosynthesis.generic_library import build_generic_library
from core_retrosynthesis.retro_validity import POLICY_PATH
from tests.condition_recommender.test_shared_core_retrieval import record
from tests.core_retrosynthesis_tests.test_external_proposal_admission import FIRST_REACTION, _row


@pytest.fixture(scope="module")
def library():
    return build_generic_library((_row(FIRST_REACTION, 1),), levels=("L0", "L1", "L2"))


def proposal(reaction=FIRST_REACTION, **kwargs):
    precursors, target = reaction.split(">>")
    return ExternalRetrosynthesisProposal(target_smiles=target, precursor_smiles=precursors, **kwargs)


def forward(validity="structurally_supported", **kwargs):
    return RetroForwardEvidence(
        starting_materials="CC=O.N", intended_product="CCN", validity=validity,
        targeted_replay_status="structurally_reproduced", intended_match="exact",
        best_competitor_product=None, warnings=(), **kwargs,
    )


def test_exact_source_grade_deduplicates_template_levels_and_is_not_probability(library):
    result = assess_retro_validity(proposal(), library)
    assert result.structural_status == "verified" and result.status == "supported"
    assert result.precedent_grade == "exact_reaction" and result.evidence_rank == 4
    assert result.operator_precedent_support.qualified_count == 1
    assert result.operator_precedent_support.distinct_reference_count == 1
    assert result.forward_status == "not_run"
    assert "supplied_recipe_not_assessed" in result.unresolved_checks
    assert result.experimental_feasibility == "not_established"
    assert result.ranking_influence == "none_advisory_only"
    assert result.score_semantics == "ordinal_evidence_not_success_probability"


def test_deterministic_identity_ignores_partner_order_and_claimed_name(library):
    normal = assess_retro_validity(proposal(), library)
    swapped = assess_retro_validity(ExternalRetrosynthesisProposal(
        target_smiles="NCC", precursor_smiles="N.CC=O", claimed_transformation="untrusted name",
    ), library)
    assert normal.assessment_id == swapped.assessment_id


@pytest.mark.parametrize("value,structural", [
    (ExternalRetrosynthesisProposal(target_smiles="invalid", precursor_smiles="CCO"), "invalid_input"),
    (ExternalRetrosynthesisProposal(target_smiles="CCO", precursor_smiles="CC=O.N",
       mapped_reaction_smiles="[CH3:1][CH:2]=[O:3].[NH3:4]>>[CH3:1][CH2:2][NH2:4]"), "contradicted"),
    (ExternalRetrosynthesisProposal(target_smiles="CCN", precursor_smiles="CC.CN"), "unresolved"),
])
def test_invalid_conflicting_and_ambiguous_mapping_cannot_receive_priority(library, value, structural):
    result = assess_retro_validity(value, library)
    assert result.structural_status == structural
    assert result.status in {"invalid_input", "contradicted", "insufficient_evidence"}
    assert result.evidence_rank is None


@pytest.mark.parametrize("validity,status,action", [
    ("structurally_supported", "supported", "prioritize_evidence_inspection"),
    ("structurally_supported_with_competition", "supported_with_cautions", "investigate_cautions"),
    ("contradicted", "contradicted", "reject_or_revise_realization"),
    ("inconclusive", "supported", "prioritize_evidence_inspection"),
    ("out_of_scope", "supported", "prioritize_evidence_inspection"),
])
def test_forward_result_does_not_rewrite_precedent_grade(library, validity, status, action):
    result = assess_retro_validity(proposal(), library, forward_evidence=forward(validity))
    assert result.precedent_grade == "exact_reaction" and result.evidence_rank == 4
    assert result.structural_status == "verified" and result.status == status
    assert result.suggested_action == action and result.forward_execution_status == "completed"
    if validity in {"inconclusive", "out_of_scope"}:
        assert f"forward:{validity}" in result.unresolved_checks


@pytest.mark.parametrize("execution", ["timed_out", "error", "cancelled"])
def test_incomplete_forward_execution_remains_unresolved(library, execution):
    result = assess_retro_validity(proposal(), library, forward_execution_status=execution)
    assert result.status == "supported" and result.forward_status == execution
    assert f"forward:{execution}" in result.unresolved_checks


def test_foreign_structure_or_recipe_forward_evidence_is_rejected(library):
    for evidence in (replace(forward(), intended_product="CCO"),
                     replace(forward(), starting_materials="CCC=O.N"),
                     forward(audited_recipe={"temperature_c": 100})):
        with pytest.raises(ValueError, match="match this realization"):
            assess_retro_validity(proposal(), library, forward_evidence=evidence)


def test_real_recipe_conflict_overrides_exact_precedent(library):
    reaction = "BrCCc1ccc(C=O)cc1.N>>NCCc1ccc(C=O)cc1"
    recipe = {"oxidants": [{"identity_status": "resolved", "role_status": "assigned",
                           "primary_role": "oxidant", "roles": [{"role_id": "oxidant"}]}]}
    result = assess_retro_validity(proposal(reaction, proposed_conditions=recipe), library)
    assert result.recipe_assessment.status == "conflict"
    assert result.status == "contradicted"


def test_no_operator_match_can_still_have_corpus_backed_support(tmp_path, library):
    reaction = "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1"
    index = build_generic_index([record(1, reaction)])
    path = tmp_path / "shared.sqlite"
    build_shared_core_index(index, path)
    result = assess_retro_validity(proposal(reaction), library, condition_index=index,
                                   shared_core_index=load_shared_core_index(path, index))
    assert not result.structural_assessment.operator_matches
    assert result.precedent_grade == "exact_reaction" and result.evidence_rank == 4
    assert result.corpus_precedent_support.matches[0].yield_pct == 70


def test_operator_only_support_is_exploratory(library):
    bare = replace(library, templates=tuple(replace(template, precedents=()) for template in library.templates))
    result = assess_retro_validity(proposal(), bare)
    assert result.structural_status == "verified"
    assert result.precedent_grade == "operator_only" and result.evidence_rank == 0
    assert result.status == "insufficient_evidence" and result.suggested_action == "seek_evidence"


@pytest.mark.parametrize("mutation", [
    lambda data: data.update(ranking_influence="rerank"),
    lambda data: data["evidence_levels"][0].update(rank=True),
    lambda data: data.update(structural_gates=[]),
    lambda data: data["limits"].update(candidate_limit=0),
])
def test_versioned_policy_rejects_weakened_gates_and_invalid_ranks(mutation):
    assert load_retro_validity_policy().definition_id == "retro_validity.v1"
    data = json.loads(POLICY_PATH.read_text("utf-8"))
    mutation(data)
    with pytest.raises(ValueError):
        validate_retro_validity_policy(data)
