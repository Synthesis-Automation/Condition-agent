"""Scientific summaries retain domain meaning and expose their preview limits."""

from copy import deepcopy
from dataclasses import asdict, replace
import json

import pytest

from chem_coworker.scientific_workspace.call_summaries import summarize_call
from condition_recommender import assess_reaction_recipe
from condition_recommender.models import GenericRecommendationResult
from core_retrosynthesis.external_route_admission import ExternalRouteProposal, assess_external_route_proposal
from core_retrosynthesis.generic_library import build_generic_library
from reactive_taxonomy import audit_target, featurize_reaction
from tests.condition_recommender.test_condition_constraints import _recommendation
from tests.core_retrosynthesis_tests.test_external_proposal_admission import _row, _route_value, FIRST_REACTION, SECOND_REACTION


def _call(operation, result):
    return {"operation": operation, "execution_status": "completed", "result": result,
            "duration_seconds": 1.234, "result_bytes": 9876}


@pytest.fixture(scope="module")
def route_record():
    library = build_generic_library((_row(FIRST_REACTION, 1), _row(SECOND_REACTION, 2)), levels=("L0",))
    assessment = assess_external_route_proposal(ExternalRouteProposal.from_dict(_route_value()), library)
    return {
        "schema_version": "route_investigation.v1", "origin": "agent_proposal",
        "review_status": "unreviewed", "assessment": assessment.to_dict(),
        "experimental_feasibility": "not_established",
        "assessment_options": {"include_conditions": False, "include_forward": False},
        "material_constraints": {"status": "unknown", "stock_availability": "not_assessed", "blocked_leaf_smiles": []},
        "proposed_recipe_assessments": {"step-1": None, "step-2": None},
        "limitations": ["A structurally assessed proposal is not experimental evidence."],
    }


def test_molecule_summary_exposes_identity_and_stereo_without_full_graph():
    result = audit_target("N[C@@H](C)C(=O)O").to_dict()
    summary = summarize_call(_call("analyze_molecule", result))
    projected = summary["result_summary"]
    assert projected["canonical_smiles"] == result["canonical_smiles"]
    assert projected["valid"] is True
    assert projected["stereocenters"][0]["assignment"] == result["stereocenters"][0]["assignment"]
    assert projected["warnings"] == list(result["warnings"])
    assert summary["duration_seconds"] == 1.234
    assert summary["result_bytes"] == 9876
    assert "atom_indices" not in json.dumps(projected)


def test_reaction_summary_preserves_independent_observation_and_interpretation():
    result = featurize_reaction(FIRST_REACTION).to_dict()
    summary = summarize_call(_call("analyze_reaction", result))["result_summary"]
    assert summary["evidence_quality"] == result["evidence_quality"]
    assert summary["observation"]["evidence_quality"] == result["observation"]["evidence_quality"]
    assert summary["interpretation"]["warnings"] == list(result["interpretation"]["warnings"])
    assert summary["reaction_completeness"]["status"] == result["reaction_completeness"]["status"]
    assert "reactants" not in summary


@pytest.mark.parametrize("reaction", ["CC>>CCN", "invalid>>CCN"])
def test_recipe_unknown_or_invalid_is_not_relabelled_as_incompatible(reaction):
    result = asdict(assess_reaction_recipe(reaction, {}))
    assert result["status"] in {"unknown", "invalid_input"}
    summary = summarize_call(_call("assess_recipe", result))
    assert summary["execution_status"] == "completed"
    assert summary["result_summary"]["status"] == result["status"]
    assert summary["result_summary"]["compatible"] is False
    assert summary["result_summary"]["hard_conflicts"] == []
    assert summary["result_summary"]["unresolved_requirements"] == ["VERIFIED_REACTION_SIGNATURE_REQUIRED"]


def test_recommendation_view_exposes_precedent_scope_cautions_and_preview_counts():
    recommendation = replace(_recommendation(1, "cas:7440-50-8"),
                             precedent_reaction_ids=("observed-reaction",),
                             precedent_reference_ids=("patent-example",),
                             cautions=("Analogue conditions do not prove transfer.",))
    result = GenericRecommendationResult(
        query_reaction_smiles=SECOND_REACTION, valid=True, retrieval_level="generic_signature",
        candidate_count=200, compatible_candidate_count=14,
        recommendations=tuple(replace(recommendation, rank=index + 1) for index in range(8)),
    ).to_dict()
    summary = summarize_call(_call("recommend_conditions", result))
    projected = summary["result_summary"]
    assert projected["candidate_count"] == 200
    assert projected["recommendations"][0]["precedent_reference_ids"] == ["patent-example"]
    assert projected["recommendations"][0]["cautions"] == list(recommendation.cautions)
    assert len(projected["recommendations"]) == 3
    metadata = next(row for row in summary["inspection"]["collections"] if row["path"] == "$.result.recommendations")
    assert metadata == {"path": "$.result.recommendations", "total": 8, "shown": 3, "omitted": 5, "count_scope": "saved_result"}
    assert "resolved_recipe" not in projected["recommendations"][0]


def test_route_step_gates_keep_unknown_not_run_and_out_of_scope(route_record):
    assessment = deepcopy(route_record["assessment"]["step_assessments"][0]["assessment"])
    assessment["gates"] = [
        {"gate_id": "mapped_gate", "status": status, "summary": status, "warnings": ["Review source evidence."]}
        for status in ("UNKNOWN", "not_run", "out_of_scope", "unresolved")
    ]
    summary = summarize_call(_call("assess_route_step", {
        "assessment": assessment, "experimental_feasibility": "not_established",
        "proposed_recipe_assessment": None,
    }))["result_summary"]
    assert [gate["status"] for gate in summary["assessment"]["gates"]] == ["UNKNOWN", "not_run", "out_of_scope", "unresolved"]
    assert summary["assessment"]["status"] == assessment["status"]
    assert summary["assessment"]["gates"][0]["warnings"] == ["Review source evidence."]
    assert summary["proposed_recipe_assessment"] is None


def test_revision_preserves_limits_material_uncertainty_and_step_identities(route_record):
    result = deepcopy(route_record)
    result["revision"] = {"improvement_status": "not_automatically_established", "reassessment_scope": "all_steps_and_route_topology",
                          "preserved_step_ids": ["step-2"], "risks": ["Unproven operating conditions."]}
    summary = summarize_call(_call("revise_route_branch", result))["result_summary"]
    assert summary["experimental_feasibility"] == "not_established"
    assert summary["assessment"]["actionable"] is False
    assert summary["assessment"]["status"] == result["assessment"]["status"]
    assert summary["assessment"]["step_assessments"][0]["external_step_id"] == result["assessment"]["step_assessments"][0]["external_step_id"]
    assert summary["material_constraints"]["status"] == "unknown"
    assert summary["material_constraints"]["stock_availability"] == "not_assessed"
    assert summary["revision"] == result["revision"]
    assert summary["proposed_recipe_assessments"] == {"step-1": None, "step-2": None}
    assert "admitted_route_tree" not in summary["assessment"]


def test_comparison_does_not_choose_a_winner_or_merge_different_statuses():
    result = {"ranking": "not_performed", "experimental_feasibility": "not_established", "alternatives": [
        {"source_ref": "before", "route_id": "r1", "status": "invalid", "step_count": 2,
         "material_constraints": {"status": "unknown", "stock_availability": "not_assessed"}},
        {"source_ref": "after", "route_id": "r2", "status": "partially_supported", "step_count": 3,
         "material_constraints": {"status": "satisfied_for_declared_constraints", "stock_availability": "not_assessed"}},
    ]}
    summary = summarize_call(_call("compare_route_proposals", result))["result_summary"]
    assert summary["ranking"] == "not_performed"
    assert [item["status"] for item in summary["alternatives"]] == ["invalid", "partially_supported"]
    assert "winner" not in summary


def test_pagination_distinguishes_saved_rows_from_total_dataset_records():
    result = {"total": 81, "offset": 20, "next_offset": 40, "records": [
        {"reaction_id": f"r{i}", "reference_id": "source", "huge_graph": {"atoms": [0] * 100}}
        for i in range(20)
    ]}
    summary = summarize_call(_call("get_precedents", result))
    assert summary["result_summary"]["total"] == 81
    assert summary["result_summary"]["next_offset"] == 40
    assert summary["inspection"]["collections"][0]["total"] == 20
    assert "huge_graph" not in json.dumps(summary)


@pytest.mark.parametrize("result", [None, 1, "unexpected output", ["unexpected list"], {"status": "UNKNOWN", "novel": {"graph": 1}}])
def test_unexpected_results_are_inspectable_without_fabricated_success(result):
    summary = summarize_call(_call("future_operation", result))
    assert summary["inspection"]["result_present"] is True
    assert summary["inspection"]["result_type"] == type(result).__name__
    assert summary["inspection"]["result_path"] == "$.result"
    assert summary["inspection"]["projection_only"] is True
    if isinstance(result, dict):
        assert summary["result_summary"]["status"] == "UNKNOWN"
        assert "novel" in summary["inspection"]["result_fields"]


def test_failed_calls_preserve_error_without_inventing_scientific_result():
    payload = {"operation": "recommend_conditions", "execution_status": "error",
               "error": {"type": "FileNotFoundError", "message": "Missing configured condition_index"}}
    summary = summarize_call(payload)
    assert summary["error"] == payload["error"]
    assert summary["execution_status"] == "error"
    assert summary["result_summary"] == {}
    assert summary["inspection"]["result_present"] is False
    assert "duration_seconds" not in summary


def test_large_or_malformed_fields_are_bounded_and_do_not_modify_source():
    result = {"warnings": ["warning" * 10000] * 100, "valid": False,
              "recommendations": [{"rank": index, "cautions": ["risk" * 10000] * 100,
                                   "explanation": {"unexpected": {"nested": ["x"] * 100, "oversize" * 1000: "x"}}}
                                  for index in range(1000)],
              "unknown_corpus": ["blob" * 1000] * 100}
    payload = _call("recommend_conditions", result)
    original = deepcopy(payload)
    summary = summarize_call(payload)
    assert payload == original
    assert len(json.dumps(summary)) < 26000
    assert summary["inspection"]["truncation_count"] > 0
    assert summary["result_summary"]["valid"] is False
    assert "unknown_corpus" not in summary["result_summary"]
    assert summarize_call(payload) == summary


def test_inspected_condition_evidence_reports_missing_information_and_compatibility():
    result = {"observation_count": 1, "distinct_reference_count": 0,
              "missing_reference_observation_ids": ["o1"], "procedure_catalog": {"availability": "catalog_unavailable"},
              "precedents": [{"observation": {"observation_id": "o1", "reaction_id": "r1", "reference_id": None},
                              "compatibility": {"status": "unknown", "compatible": False, "hard_conflicts": []},
                              "missing_operating_fields": ["temperature_c", "time_h"]}]}
    summary = summarize_call(_call("inspect_condition_precedents", result))["result_summary"]
    assert summary["procedure_catalog"]["availability"] == "catalog_unavailable"
    assert summary["precedents"][0]["compatibility"]["status"] == "unknown"
    assert summary["precedents"][0]["missing_operating_fields"] == ["temperature_c", "time_h"]
