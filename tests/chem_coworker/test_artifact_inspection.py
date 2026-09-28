"""Targeted saved-result previews remain bounded, literal, and scientifically neutral."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace


@pytest.fixture
def workspace(tmp_path: Path) -> ScientificWorkspace:
    # Inspection is read-only and must remain usable when an old baseline can no
    # longer execute. No datasets, domain engine or model are needed for this view.
    InvestigationStore.create(tmp_path / "run", objective="Inspect saved results", baseline={})
    return ScientificWorkspace(tmp_path / "run")


def _save(workspace: ScientificWorkspace, payload):
    return workspace.store.append("call", payload).artifact_ref


def test_list_pages_preserve_stored_indices_counts_and_scientific_context(workspace):
    payload = {
        "operation": "recommend_conditions", "execution_status": "completed",
        "result": {"valid": False, "status": "unknown", "warnings": ["Mapping remains unresolved"],
                   "experimental_feasibility": "not_established", "candidate_count": 900,
                   "recommendations": [{"rank": index + 1, "recipe_id": f"r{index}",
                                         "resolved_recipe": {"components": ["NaOH"]}} for index in range(8)]},
    }
    reference = _save(workspace, payload)
    view = workspace.inspect_artifact(reference, path=("result", "recommendations"), offset=3, limit=2)
    assert view["preview"] == payload["result"]["recommendations"][3:5]
    assert view["path"] == ["result", "recommendations"]
    assert view["json_path"] == '$["result"]["recommendations"]'
    assert view["artifact_ref"] == reference
    assert view["page"] == {
        "total": 8, "shown": 2, "offset": 3, "limit": 2, "next_offset": 5,
        "omitted_before": 3, "omitted_after": 3, "count_scope": "saved_artifact", "value_type": "list",
    }
    assert view["context"][0]["fields"]["execution_status"] == "completed"
    science = view["context"][1]["fields"]
    assert science["valid"] is False and science["status"] == "unknown"
    assert science["warnings"] == payload["result"]["warnings"]
    assert science["experimental_feasibility"] == "not_established"
    assert view["inspection"]["truncated"] is True
    assert any(item["reason"] == "collection_page" for item in view["inspection"]["truncations"])
    nested = workspace.inspect_artifact(reference, ("result", "recommendations", 3, "resolved_recipe"))
    assert nested["preview"] == {"components": ["NaOH"]}
    assert nested["json_path"] == '$["result"]["recommendations"][3]["resolved_recipe"]'


def test_mapping_pages_follow_saved_keys_and_literal_punctuation_is_not_evaluated(workspace):
    reference = _save(workspace, {"result": {"a.b[0]": {"alpha": 1, "beta": 2, "gamma": 3}}})
    view = workspace.inspect_artifact(reference, ["result", "a.b[0]"], offset=1, limit=1)
    assert view["preview"] == {"beta": 2}
    assert view["json_path"] == '$["result"]["a.b[0]"]'
    assert view["page"]["total"] == 3 and view["page"]["next_offset"] == 2
    with pytest.raises(KeyError):
        workspace.inspect_artifact(reference, ("result", "__class__"))


def test_failed_call_keeps_error_without_a_fabricated_result(workspace):
    payload = {"operation": "disconnect_target", "execution_status": "error",
               "error": {"type": "FileNotFoundError", "message": "retro_library is unavailable"}}
    reference = _save(workspace, payload)
    view = workspace.inspect_artifact(reference, ("error",))
    assert view["preview"] == payload["error"]
    assert view["context"][0]["fields"]["execution_status"] == "error"
    with pytest.raises(KeyError):
        workspace.inspect_artifact(reference, ("result",))


def test_preview_discloses_nested_and_text_limits_without_mutating_evidence(workspace, monkeypatch):
    payload = {"operation": "future_operation", "execution_status": "completed", "result": {
        "warnings": ["Warning " * 3000], "status": "unresolved",
        "rows": [{"message": "x" * 10000, "nested": {"deeper": {"deeper": {"secret": "not in preview"}}}}
                 for _ in range(80)],
    }}
    reference = _save(workspace, payload)
    before = deepcopy(workspace.store.read_artifact(reference))
    events_before = workspace.store.events()
    monkeypatch.setattr(workspace.operations, "invoke", lambda *_: pytest.fail("Inspection must not execute chemistry"))
    view = workspace.inspect_artifact(reference, ("result", "rows"))
    assert len(json.dumps(view).encode()) < 24000
    assert view["page"]["total"] == 80 and view["page"]["shown"] == 5
    assert view["inspection"]["truncation_count"] > 0
    assert any(item["reason"] == "text_preview" for item in view["inspection"]["truncations"])
    assert any(item["reason"] == "nested_detail" for item in view["inspection"]["truncations"])
    assert workspace.inspect_artifact(reference, ("result", "rows")) == view
    assert workspace.store.read_artifact(reference) == before
    assert workspace.store.events() == events_before
    warning = workspace.inspect_artifact(reference, ("result", "warnings", 0))
    assert warning["preview"].startswith("Warning ")
    assert warning["inspection"]["truncated"] is True


def test_serialized_cap_reports_omission_even_for_giant_mapping_keys(workspace):
    reference = _save(workspace, {"result": {"huge_key" * 10000: "small value"}})
    view = workspace.inspect_artifact(reference, ("result",))
    assert len(json.dumps(view, ensure_ascii=False).encode()) < 24000
    assert view["preview"] == {"summary_omitted": True}
    assert view["page"]["shown"] == 0 and view["page"]["selected"] == 1
    assert view["inspection"]["metadata_omitted"] is True
    assert view["inspection"]["truncations"][0]["reason"] == "serialized_preview_budget"


def test_scalar_inspection_keeps_compact_ancestor_cautions_without_repeating_trees(workspace):
    payload = {"operation": "assess_route_proposal", "execution_status": "completed", "result": {
        "status": "unresolved", "warnings": [{"full_tree": {"rows": ["large" * 500] * 100}}] * 20,
        "evidence_refs": [f"sha256:{index:064x}" for index in range(20)],
        "experimental_feasibility": "not_established",
        "assessment": {"status": "unknown", "warnings": ["Missing source contributor"],
                       "actionable": False, "detail": "selected value"},
    }}
    reference = _save(workspace, payload)
    view = workspace.inspect_artifact(reference, ("result", "assessment", "detail"))
    assert view["preview"] == "selected value"
    assert len(json.dumps(view)) < 6000
    assert "full_tree" not in json.dumps(view)
    assert view["context"][1]["fields"]["status"] == "unresolved"
    assert view["context"][2]["fields"]["actionable"] is False
    assert view["context"][2]["fields"]["warnings"] == ["Missing source contributor"]
    assert view["context"][1]["fields"]["evidence_refs"] == payload["result"]["evidence_refs"][:2]
    assert any(item["reason"] == "ancestor_context_detail" for item in view["inspection"]["truncations"])
    assert view["inspection"]["truncated"] is True
    assert workspace.store.read_artifact(reference) == payload


def test_requested_text_budget_is_independent_of_large_ancestor_warnings(workspace):
    reference = _save(workspace, {"result": {"warnings": ["Surrounding warning " * 1000] * 40,
                                            "finding": "Requested finding " * 100}})
    view = workspace.inspect_artifact(reference, ("result", "finding"))
    assert view["preview"].startswith("Requested finding ")
    assert len(view["preview"]) == 401
    assert view["inspection"]["context_text_budget_characters"] == 1200


@pytest.mark.parametrize(("path", "options"), [
    ("result.rows", {}), ((True,), {}), ((-1,), {}), (("x" * 81,), {}),
    (("result",) * 13, {}), (("result",), {"offset": -1}),
    (("result",), {"offset": True}), (("result",), {"limit": False}),
    (("result",), {"limit": 0}), (("result",), {"limit": 21}),
    (("result", "rows", "0"), {}), (("result", "scalar"), {"offset": 1}),
])
def test_invalid_paths_and_page_options_are_rejected(workspace, path, options):
    reference = _save(workspace, {"result": {"rows": ["item"], "scalar": "text"}})
    with pytest.raises(ValueError):
        workspace.inspect_artifact(reference, path, **options)


def test_empty_and_past_end_pages_distinguish_missing_paths_from_empty_results(workspace):
    reference = _save(workspace, {"result": {"rows": [], "null": None}})
    view = workspace.inspect_artifact(reference, ("result", "rows"), offset=10)
    assert view["preview"] == []
    assert view["page"]["total"] == 0 and view["page"]["shown"] == 0
    assert view["page"]["next_offset"] is None
    assert workspace.inspect_artifact(reference, ("result", "null"))["preview"] is None
    with pytest.raises(IndexError):
        workspace.inspect_artifact(reference, ("result", "rows", 0))


def test_artifact_reference_and_checksum_validation_cannot_be_bypassed(workspace):
    with pytest.raises(ValueError, match="Invalid artifact reference"):
        workspace.inspect_artifact("../outside.json")
    reference = _save(workspace, {"result": {"status": "unknown"}})
    artifact = workspace.store.root / "artifacts" / (reference.removeprefix("sha256:") + ".json")
    artifact.write_text('{"result":{"status":"verified"}}', encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        workspace.inspect_artifact(reference, ("result",))
