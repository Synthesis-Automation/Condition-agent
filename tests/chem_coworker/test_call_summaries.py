"""Scientific summaries retain domain meaning and expose their preview limits."""

from copy import deepcopy
from dataclasses import asdict, replace
import json

import pytest

from chem_coworker.scientific_workspace.views.call_summaries import summarize_call
from chem_coworker.scientific_workspace.views.brief_summaries import summarize_call_brief
from condition_recommender import assess_reaction_recipe
from condition_recommender.models import GenericRecommendationResult
from core_retrosynthesis.external_route_admission import (
    ExternalRouteProposal,
    assess_external_route_proposal,
)
from core_retrosynthesis.generic_library import build_generic_library
from reactive_taxonomy import audit_target, featurize_reaction
from tests.condition_recommender.test_condition_constraints import _recommendation
from tests.core_retrosynthesis_tests.test_external_proposal_admission import (
    _row,
    _route_value,
    FIRST_REACTION,
    SECOND_REACTION,
)


def _call(operation, result):
    return {"operation": operation, "execution_status": "completed", "result": result,
            "duration_seconds": 1.234, "result_bytes": 9876}


@pytest.mark.parametrize("summarize", [summarize_call, summarize_call_brief])
def test_disconnection_summary_retains_precursor_burden_without_reranking(summarize):
    complexity = {
        "definition_id": "strategic_complexity.v1", "evidence": "molecular_graph_only",
        "graph_complexity_reduction_fraction": -0.4, "extra_precursor_heavy_atom_count": 12,
        "warnings": ["RETROSYNTHETIC_COMPLEXITY_INCREASE", "STRATEGIC_COMPLEXITY_MAPPING_UNAVAILABLE"],
    }
    result = {"strategies": [{"strategy_id": "first", "strategy_rank": 1, "representative": {
        "precursor_smiles": "CCO", "strategic_class": "complexity_increasing",
        "strategic_complexity_score": 0.0, "strategic_candidate": False,
        "strategic_complexity": complexity,
    }}, {"strategy_id": "second", "strategy_rank": 2, "representative": {"precursor_smiles": "CO"}}]}
    original = deepcopy(result)
    summary = summarize(_call("disconnect_target", result))["result_summary"]
    assert [s["strategy_id"] for s in summary["strategies"]] == ["first", "second"]
    representative = summary["strategies"][0]["representative"]
    assert representative["strategic_complexity"] == complexity
    assert representative["strategic_class"] == "complexity_increasing"
    assert "strategic_complexity" not in summary["strategies"][1]["representative"]
    assert result == original


def test_recipe_summary_exposes_coverage_without_reading_complete_artifact() -> None:
    assessment = asdict(assess_reaction_recipe("CCBr.N>>CCN", {
        "solvents": [{"identity_status": "resolved", "substance_id": "cas:64-17-5"}],
    }))
    call = _call("assess_recipe", assessment)
    full = summarize_call(call)["result_summary"]["coverage"]
    brief = summarize_call_brief(call)["result_summary"]["coverage"]
    for coverage in (full, brief):
        assert coverage["capability_status"] == "not_covered"
        assert coverage["condition_identity_status"] == "resolved"
        assert coverage["score_meaning"] == "absence_of_known_conflicts_not_success_probability"
    assert full["evaluated_hard_conflict_rule_ids"]


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


def test_comparison_brief_keeps_ambiguity_and_timeout_without_atom_tables():
    from reactive_taxonomy import compare_molecules

    result = replace(compare_molecules('Cc1ccccc1', 'Clc1ccccc1'),
                     status='partial_timeout', search_timed_out=True).to_dict()
    brief = summarize_call_brief(_call('compare_molecules', result))
    summary = brief['result_summary']
    assert summary['status'] == 'partial_timeout' and summary['search_timed_out']
    assert summary['alignment_ambiguous'] and summary['alignments_truncated']
    assert summary['stereo_relationship'] == 'not_compared_different_graphs'
    assert summary['alignment_count_observed'] == 12
    assert 'atoms' not in summary['left']
    assert len(json.dumps(brief)) < 4000


def test_focused_brief_keeps_alternative_sites_and_descriptor_provenance():
    from reactive_taxonomy import inspect_reactive_sites

    result = inspect_reactive_sites('BrCCCCCCO')
    bromine = next(a.atom_id for a in result.molecule.atoms if a.element == 'Br')
    focused = inspect_reactive_sites('BrCCCCCCO', [bromine], radius=0)
    brief = summarize_call_brief(_call('inspect_reactive_sites', focused.to_dict()))
    summary = brief['result_summary']
    assert summary['sites_count'] == len(focused.sites)
    assert summary['other_sites_count'] == len(focused.other_sites) > 0
    assert summary['other_sites']
    assert summary['site_profiles'][0]['steric']['evidence']['method']
    assert summary['site_profiles'][0]['status'] == 'derived'
    assert 'atoms' not in summary['molecule']
    assert len(json.dumps(brief)) < 6000


def test_inspection_brief_preserves_unspecified_stereo_warning():
    from reactive_taxonomy import inspect_reactive_sites

    result = inspect_reactive_sites('CC=CC').to_dict()
    brief = summarize_call_brief(_call('inspect_reactive_sites', result))
    assert 'UNSPECIFIED_STEREOCHEMISTRY' in brief['result_summary']['warnings']


def test_brief_selects_failure_after_long_passing_gate_prefix():
    result = {"assessment": {"status": "unresolved", "actionable": False,
              "gates": [{"gate_id": f"g{i}", "status": "pass"} for i in range(25)] + [
                  {"gate_id": "missing_atom_donor", "status": "failed", "warnings": ["Missing oxygen source"]}
              ]}}
    brief = summarize_call_brief(_call("assess_route_step", result))["result_summary"]
    assert brief["assessment"]["non_pass_gate_count"] == 1
    assert brief["assessment"]["non_pass_gates"][0]["gate_id"] == "missing_atom_donor"
    assert "Missing oxygen source" in brief["assessment"]["non_pass_gates"][0]["warnings"]


def test_brief_counts_saved_warnings_and_steps_before_previewing():
    result = {"warnings": [f"warning {i}" for i in range(23)], "assessment": {
        "status": "unresolved", "step_assessments": [
            {"external_step_id": f"s{i}", "assessment": {"status": "supported", "actionable": True}}
            for i in range(12)
        ] + [{"external_step_id": "late-failure", "assessment": {"status": "unresolved", "actionable": False}}],
    }}
    brief = summarize_call_brief(_call("assess_route_proposal", result))["result_summary"]
    assert len(brief["warnings"]) == 3 and brief["warnings_count"] == 23
    assert brief["assessment"]["step_count"] == 13
    assert brief["assessment"]["steps"][0]["step_id"] == "late-failure"


def test_brief_recommendation_contains_usable_recipe_and_source_identity():
    result = {"valid": True, "recommendations": [{
        "rank": 1, "recipe_id": "recipe-1", "reference_support": 1, "observation_support": 25,
        "precedent_reaction_ids": ["r1"], "resolved_recipe": {
            "catalysts": [{"canonical_name": "Pd(OAc)2", "primary_role": "metal_catalyst", "identity_status": "resolved"}],
            "bases": [{"canonical_name": "Potassium carbonate", "identity_status": "resolved"}],
            "temperature_c": 80, "time_h": None,
        },
    }]}
    original = deepcopy(result)
    brief = summarize_call_brief(_call("recommend_conditions", result))["result_summary"]
    record = brief["recommendations"][0]
    assert record["recipe"]["components"][0]["name"] == "Pd(OAc)2"
    assert record["recipe"]["temperature_c"] == 80 and "time_h" not in record["recipe"]
    assert record["reference_support"] == 1 and record["observation_support"] == 25
    assert record["precedent_reaction_ids"] == ["r1"]
    assert result == original


def test_unresolved_recipe_component_remains_visible_after_preview_limit():
    recipe = {"recipe_id": "recipe", "solvents": [
        {"canonical_name": f"solvent-{i}", "identity_status": "resolved"} for i in range(9)
    ], "other_components": [{"raw_identifier": "ambiguous raw name", "identity_status": "ambiguous",
                             "warnings": ["CONDITION_IDENTITY_UNCERTAINTY"]}], "stages": [{}, {}]}
    result = summarize_call_brief(_call("resolve_recipe", recipe))["result_summary"]["recipe"]
    assert result["component_count"] == 10 and result["omitted_component_count"] == 4
    assert result["components"][0]["identity_status"] == "ambiguous"
    assert result["components"][0]["saved_path"] == ["other_components", 0]
    assert result["stage_count"] == 2 and result["stage_details_omitted"]


def test_procedure_brief_reports_text_availability_without_dumping_text():
    source = {"records": [
        {"reaction_id": "r1", "observation_id": "o1", "procedure_text": "reported text " * 5000},
        {"reaction_id": "r1", "observation_id": "o2", "procedure_text": None},
    ], "missing_reaction_ids": ["missing"]}
    result = summarize_call_brief(_call("get_procedures", source))["result_summary"]
    assert result["saved_record_count"] == 2
    assert result["records"][0]["procedure_text_present"]
    assert result["records"][0]["procedure_characters"] == len(source["records"][0]["procedure_text"])
    assert not result["records"][1]["procedure_text_present"]
    assert result["missing_reaction_ids"] == ["missing"]
    assert "reported text" not in json.dumps(result)


def test_precedent_brief_keeps_pagination_and_record_status():
    source = {"records": [{"reaction_id": f"r{i}", "observation_id": f"o{i}",
                           "condition_status": "unresolved", "huge_graph": [0] * 10000} for i in range(10)],
              "total": 30, "offset": 10, "next_offset": 20, "missing_reaction_ids": ["missing"]}
    result = summarize_call_brief(_call("get_precedents", source))["result_summary"]
    assert result["total"] == 30 and result["next_offset"] == 20 and result["saved_record_count"] == 10
    assert result["records"][0]["observation_id"] == "o0"
    assert result["records"][0]["condition_status"] == "unresolved"
    assert "huge_graph" not in json.dumps(result)


def test_brief_preserves_long_smiles_and_flags_extreme_string_truncation():
    result = {"valid": True, "canonical_smiles": "C" * 240}
    assert summarize_call_brief(_call("analyze_molecule", result))["result_summary"]["canonical_smiles"] == result["canonical_smiles"]
    result["canonical_smiles"] = "C" * 5000
    summary = summarize_call_brief(_call("analyze_molecule", result))["result_summary"]
    assert summary["canonical_smiles_truncated"] and len(summary["canonical_smiles"]) == 2001


def test_brief_nested_values_stay_bounded_and_source_untouched():
    source = {"valid": False, "warnings": [{"one": {"two": ["x" * 100000] * 50}}] * 25,
              "recommendations": [{"rank": 1, "cautions": ["risk" * 10000] * 20}] * 100}
    original = deepcopy(source)
    summary = summarize_call_brief(_call("recommend_conditions", source))
    assert source == original and summary["result_summary"]["valid"] is False
    assert summary["result_summary"]["warnings_count"] == 25
    assert len(json.dumps(summary)) < 8000


def test_failed_brief_does_not_invent_zero_candidates():
    summary = summarize_call_brief({"operation": "recommend_conditions", "execution_status": "error",
                                    "error": {"type": "FileNotFoundError", "message": "Missing index"}})
    assert summary["result_summary"] == {}
    assert summary["error"]["type"] == "FileNotFoundError"


def test_brief_keeps_empty_root_warning_contract_and_unknown_edit_counts():
    result = {"valid": True, "warnings": [], "reaction_signature": {"formed_bond_types": []}}
    summary = summarize_call_brief(_call("analyze_reaction", result))["result_summary"]
    assert summary["warnings"] == []
    assert "hydrogen_changes_count" not in summary["bond_changes"]


def test_brief_forward_execution_shows_latest_stage():
    stages = [{"stage": f"phase-{i}", "status": "completed"} for i in range(10)]
    stages.append({"stage": "predict", "status": "running"})
    summary = summarize_call_brief(_call("assess_route_step_forward", {"execution": {"stages": stages}}))
    assert summary["result_summary"]["execution"]["stages"][-1]["stage"] == "predict"
    assert summary["result_summary"]["execution"]["stage_count"] == 11


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


def test_six_step_overview_keeps_each_status_without_repeating_structures(route_record):
    result = deepcopy(route_record)
    template = result["assessment"]["step_assessments"][0]
    result["assessment"]["step_assessments"] = []
    for index in range(6):
        step = deepcopy(template)
        step["external_step_id"] = f"s{index + 1}"
        step["assessment"].update({"status": "unresolved" if index == 5 else "supported",
                                   "warnings": ["Inspect the missing atom donor"] if index == 5 else [],
                                   "canonical_precursor_smiles": "C" * 500})
        result["assessment"]["step_assessments"].append(step)
    summary = summarize_call(_call("assess_route_proposal", result))
    steps = summary["result_summary"]["assessment"]["step_assessments"]
    assert len(steps) == 6
    assert steps[-1]["assessment"]["status"] == "unresolved"
    assert steps[-1]["assessment"]["warnings"] == ["Inspect the missing atom donor"]
    assert all("canonical_precursor_smiles" not in step["assessment"] for step in steps)
    assert summary["result_summary"]["assessment"]["topology_gates"]
    assert not any(row["path"].endswith(".warnings") for row in summary["inspection"]["collections"])


def test_brief_route_view_keeps_failed_decisions_and_step_statuses(route_record):
    result = deepcopy(route_record)
    result["assessment"]["topology_gates"] = [
        {"gate_id": "parse", "status": "pass", "summary": "Parsed"},
        {"gate_id": "admission", "status": "unresolved", "summary": "No exact operator"},
    ]
    result["assessment"]["step_assessments"][0]["assessment"]["warnings"] = [
        "Exact precedent is missing."
    ]
    original = deepcopy(result)
    brief = summarize_call_brief(_call("assess_route_proposal", result))
    assessment = brief["result_summary"]["assessment"]
    assert result == original
    assert assessment["status"] == result["assessment"]["status"]
    assert assessment["actionable"] is False
    assert assessment["non_pass_gates"][0]["gate_id"] == "admission"
    assert assessment["steps"][0]["warnings"] == ["Exact precedent is missing."]
    assert brief["inspection"]["detail_available"] is True
    assert "admitted_route_tree" not in json.dumps(brief)
    assert len(json.dumps(brief)) < 5000


def test_brief_disconnection_retains_selectivity_warning_without_large_nested_graph():
    warning = {"code": "POSSIBLE_FUNCTIONAL_GROUP_COMPETITION",
               "message": "N and O could compete.", "conditions_evaluated": False,
               "competing_outcomes": [{"atoms": list(range(1000))}]}
    result = {"valid": True, "strategies": [{"strategy_id": "s1", "strategy_rank": 1,
              "representative": {"precursor_smiles": "CC.N", "forward_validation_status": "verified_signature",
                                 "selectivity_warnings": [warning]}}]}
    brief = summarize_call_brief(_call("disconnect_target", result))
    representative = brief["result_summary"]["strategies"][0]["representative"]
    assert representative["selectivity_warnings"][0]["code"] == warning["code"]
    assert representative["selectivity_warnings"][0]["conditions_evaluated"] is False
    assert "competing_outcomes" not in json.dumps(brief)
    assert len(json.dumps(brief)) < 3000


def test_brief_fragment_search_discloses_incomplete_source_coverage():
    brief = summarize_call_brief(_call("search_fragment_precedents", {
        "search_status": "partial", "source_coverage_complete": False,
        "ranking_scope": "returned_hits_only", "output_truncated": True,
        "hits": [{"hit_id": "h1", "product_smiles": "COC",
                  "procedure_availability": "unknown"}],
    }))
    result = brief["result_summary"]
    assert result["search_status"] == "partial"
    assert result["source_coverage_complete"] is False
    assert result["ranking_scope"] == "returned_hits_only"
    assert result["output_truncated"] is True


def test_forward_timeout_keeps_question_stage_and_unresolved_result():
    summary = summarize_call({
        "operation": "assess_route_step_forward", "execution_status": "timed_out",
        "error": {"type": "TimeoutError", "message": "Forward worker exceeded 30 seconds"},
        "result": {"schema_version": "route_step_forward_investigation.v1", "execution_status": "timed_out",
                   "source_ref": "saved-route", "step_id": "s4", "question": "Could regioselectivity change the route?",
                   "assessment": None, "experimental_feasibility": "not_established",
                   "execution": {"timings": {"elapsed_seconds": 30.1, "timeout_seconds": 30},
                                 "stages": [{"stage": "load_forward_library", "status": "running", "elapsed_seconds": 0.2}],
                                 "diagnostics": {"stderr.log": "diagnostics/forward_checks/run/stderr.log"}}},
    })
    assert summary["execution_status"] == "timed_out"
    projected = summary["result_summary"]
    assert projected["assessment"] is None
    assert projected["source_ref"] == "saved-route" and projected["step_id"] == "s4"
    assert projected["question"] == "Could regioselectivity change the route?"
    assert projected["execution"]["stages"][0]["status"] == "running"
    assert projected["execution"]["timings"]["timeout_seconds"] == 30
    assert projected["execution"]["diagnostics"]["stderr.log"] == "diagnostics/forward_checks/run/stderr.log"
    assert projected["experimental_feasibility"] == "not_established"


def test_brief_forward_timeout_does_not_hide_execution_failure():
    brief = summarize_call_brief({
        "operation": "assess_route_step_forward", "execution_status": "timed_out",
        "error": {"type": "TimeoutError", "message": "Worker timed out"},
        "result": {"execution_status": "timed_out", "step_id": "s2",
                   "assessment": None, "execution": {"stages": [
                       {"stage": "load_library", "status": "running"}]}}})
    assert brief["execution_status"] == "timed_out"
    assert brief["error"]["type"] == "TimeoutError"
    assert brief["result_summary"]["step_id"] == "s2"
    assert brief["result_summary"]["execution"]["stages"][0]["status"] == "running"


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


def test_brief_comparison_keeps_distinct_alternative_statuses():
    brief = summarize_call_brief(_call("compare_route_proposals", {
        "ranking": "not_performed", "alternatives": [
            {"route_id": "r1", "status": "invalid", "step_count": 2},
            {"route_id": "r2", "status": "partially_supported", "step_count": 3},
        ]}))
    assert brief["result_summary"]["ranking"] == "not_performed"
    assert [item["status"] for item in brief["result_summary"]["alternatives"]] == [
        "invalid", "partially_supported",
    ]


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
