"""Development pilot prompts preserve agent choice and observable-only metrics."""

from pathlib import Path

from chem_coworker.scientific_workspace import ScientificWorkspace
from examples.ai_native.fragment_agent_comparison import load_cases
from examples.ai_native.fragment_loop_pilot import loop_metrics, pilot_question


def test_live_loop_cases_validate_and_do_not_supply_fragments():
    path = Path(__file__).resolve().parents[1] / "examples/ai_native/fragment_loop_cases.json"
    cases = load_cases(path)
    assert len(cases) == 2
    for case in cases:
        question = pilot_question(case, 600)
        assert case.target_smiles in question
        assert "Choose your own connected key fragment" in question
        assert "No minimum query count" in question
        assert "six fragment" in question
        assert "twelve scientific calls" in question


def test_loop_metrics_keep_decision_links_and_do_not_score_inspections(tmp_path):
    root = Path(__file__).resolve().parents[1]
    w = ScientificWorkspace.create(tmp_path / "pilot", objective="Metric regression", repository=root)
    first = w.store.append("call", {"operation": "search_fragment_precedents", "arguments": {"query": "CO"},
                                   "execution_status": "completed", "result": {"hits": []}})
    w.store.note("decision", "Narrow after a broad match", evidence_refs=(first.artifact_ref,))
    w.store.append("console_inspection", {"source_ref": first.artifact_ref,
                                         "review_status": "opened_not_adjudicated"},
                   evidence_refs=(first.artifact_ref,))
    metrics = loop_metrics(w)
    assert metrics["notes"][0]["evidence_refs"] == (first.artifact_ref,)
    assert metrics["console_inspections"][0]["review_status"] == "opened_not_adjudicated"
    assert "pending" in metrics["manual_chemistry_review"]
    assert not metrics["search_budget_exceeded"]
