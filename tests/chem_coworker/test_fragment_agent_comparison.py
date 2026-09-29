"""Matched agent runs keep provenance, failures and evaluation limitations explicit."""

import json
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.agent_runtime import AgentResult, AgentStopped
from chem_coworker.scientific_workspace.baseline import code_manifest, environment_versions
from examples.ai_native.fragment_agent_comparison import (
    ComparisonCase, collect_metrics, comparison_question, load_cases, run_comparison,
)


ROOT = Path(__file__).resolve().parents[2]
CASE = ComparisonCase("core", "CC(=O)c1ccc2c(c1)COc1ccccc1-2", "Inspect the core-construction question.")


@pytest.fixture
def baseline():
    return {"repository": str(ROOT), "code_files": code_manifest(ROOT),
            "environment": environment_versions(), "artifacts": {},
            "evaluation_partition": "development_only_not_an_untouched_evaluation"}


class ComparisonRuntime:
    """Stub only the agent; run real recorded chemistry and answer validation."""

    def __init__(self, *, fail=False, contaminate=False, fabricate=False):
        self.fail, self.contaminate, self.fabricate = fail, contaminate, fabricate
        self.requests = []

    def describe(self):
        return {"runtime": "comparison_test_double", "model": "none"}

    def run(self, **kwargs):
        self.requests.append(kwargs)
        assert kwargs["thread_id"] is None
        workspace = ScientificWorkspace(kwargs["workspace"])
        assert not workspace.recall_lessons("retrosynthesis")["enabled"]
        assert not workspace.recall_lessons("retrosynthesis")["lessons"]
        assert "construction witness" in workspace.task_guide("retrosynthesis")["text"]
        is_control = "For this comparison arm, do not use" in kwargs["prompt"]
        operation = "analyze_molecule" if is_control and not self.contaminate else "suggest_search_fragments"
        arguments = {"smiles": CASE.target_smiles} if operation == "analyze_molecule" else {"target_smiles": CASE.target_smiles}
        event = workspace.run(operation, arguments)
        if self.fail:
            raise AgentStopped("timed_out")
        return AgentResult({
            "schema_version": "scientific_answer.v2", "sources": [], "molecules": [],
            "target_molecule_ids": [], "steps": [], "routes": [], "claims": [],
            "answer_markdown": "Structural query recorded; synthesis remains unresolved.",
            "evidence_refs": ["sha256:" + "0" * 64 if self.fabricate else event.artifact_ref],
            "uncertainties": ["No inspected construction precedent"], "needs_user_input": False,
        }, f"thread-{len(self.requests)}", {"input_tokens": 10})


def test_matched_arms_use_fresh_threads_same_baseline_and_no_memory(tmp_path, baseline):
    runtime = ComparisonRuntime()
    cases = (CASE, ComparisonCase("second", CASE.target_smiles, CASE.question))
    report = run_comparison(tmp_path / "comparison", cases, baseline, runtime, 30)
    assert [r["arm"] for r in report["trials"]] == ["agent_only", "fragment_assisted", "fragment_assisted", "agent_only"]
    assert all(r["status"] == "completed" for r in report["trials"])
    assert len({r["thread_id"] for r in report["trials"]}) == 4
    assert all(pair["same_scientific_baseline"] and pair["same_requested_runtime"] for pair in report["pairs"])
    assert not any(pair["recorded_constraint_violation"] for pair in report["pairs"])
    assert report["trials"][0]["metrics"]["fragment_call_count"] == 0
    assert report["trials"][1]["metrics"]["fragment_call_count"] == 1
    assert all(r["metrics"]["manual_review"]["unsupported_route_steps"] is None for r in report["trials"])
    assert baseline["evaluation_partition"] == "development_only_not_an_untouched_evaluation"
    assert "learning_context" not in baseline  # Caller baseline was not mutated.
    for row in report["trials"]:
        workspace = ScientificWorkspace(row["directory"])
        assert workspace.store.read_artifact(row["answer_ref"])["review_status"] == "unreviewed"
    saved = json.loads((tmp_path / "comparison" / "comparison.json").read_text())
    assert saved["pairs"] == report["pairs"]
    with pytest.raises(FileExistsError):
        run_comparison(tmp_path / "comparison", cases, baseline, runtime, 30)


@pytest.mark.parametrize("option,expected", [("fail", "timed_out"), ("fabricate", "failed")])
def test_failed_or_invalid_answers_remain_in_comparison(tmp_path, baseline, option, expected):
    runtime = ComparisonRuntime(**{option: True})
    report = run_comparison(tmp_path / "comparison", (CASE,), baseline, runtime, 30)
    assert len(report["trials"]) == 2
    assert all(row["status"] == expected and row["error"] and row["answer_ref"] is None for row in report["trials"])
    assert all(row["metrics"]["scientific_call_count"] == 1 for row in report["trials"])
    assert not report["pairs"][0]["both_completed"]
    assert report["chemistry_review_status"] == "pending"


def test_control_using_fragment_tool_is_flagged_not_silently_scored(tmp_path, baseline):
    report = run_comparison(tmp_path / "comparison", (CASE,), baseline, ComparisonRuntime(contaminate=True), 30)
    assert report["pairs"][0]["recorded_constraint_violation"]
    assert report["trials"][0]["metrics"]["recorded_arm_violations"][0]["operation"] == "suggest_search_fragments"


def test_returned_construction_is_not_scored_as_inspected_useful_evidence(tmp_path, baseline):
    report = run_comparison(tmp_path / "comparison", (CASE,), baseline, ComparisonRuntime(), 30)
    workspace = ScientificWorkspace(report["trials"][1]["directory"])
    # Report-metric fixture: duplicate construction observations and a retained-only hit.
    for _ in range(2):
        workspace.store.append("call", {"operation": "search_fragment_precedents", "execution_status": "completed",
                                        "result": {"hits": [{"observation_id": "one", "relationships": ["constructed", "unresolved"]},
                                                            {"observation_id": "two", "relationships": ["carried_through"]}]}})
    metrics = collect_metrics(workspace, "fragment_assisted", None)
    assert metrics["returned_construction_observation_ids"] == ["one"]
    assert metrics["manual_review"]["useful_inspected_construction_precedents"] is None


def test_cases_and_prompts_preserve_common_task_budget_and_optional_guidance(tmp_path):
    cases = load_cases(ROOT / "examples/ai_native/fragment_agent_cases.json")
    assert len(cases) == 2
    for arm in ("agent_only", "fragment_assisted"):
        prompt = comparison_question(CASE, arm, 240)
        assert CASE.target_smiles in prompt and CASE.question in prompt
        assert "at most four recorded scientific" in prompt and "240 seconds" in prompt
        assert "Do not read other" in prompt
    invalid = tmp_path / "cases.json"
    invalid.write_text(json.dumps({"schema_version": "fragment_agent_cases.v1", "cases": [
        {"case_id": "../escape", "target_smiles": "CO", "question": "Inspect core"}]}))
    with pytest.raises(ValueError, match="safe lowercase"):
        load_cases(invalid)
    with pytest.raises(ValueError, match="Unknown"):
        comparison_question(CASE, "unknown", 240)
