"""Lossless route storage, cross-call reuse, legacy access, and integrity failures."""

from copy import deepcopy
import hashlib
import json

from fastapi.testclient import TestClient
import pytest

from app.web_api.main import create_app
from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.artifact_storage import STORAGE_VERSION, SECTION_VERSION
from chem_coworker.scientific_workspace.call_summaries import summarize_call_brief
from chem_coworker.scientific_workspace.conversation import ConversationService
from chem_coworker.scientific_workspace.store import canonical_bytes, _write_json


@pytest.fixture
def store(tmp_path):
    return InvestigationStore.create(tmp_path / ("a" * 32), objective="Compact storage regression", baseline={})


def payload():
    """Representative assessment shape with bulky library membership and evidence."""
    return {
        "operation": "assess_route_proposal", "execution_status": "completed",
        "arguments": {"proposal": {"target_smiles": "CCN"}},
        "result": {"schema_version": "route_investigation.v1", "proposal": {"target_smiles": "CCN"},
                   "assessment": {"status": "partially_supported", "warnings": ["A real unresolved gate"],
                       "step_assessments": [{"external_step_id": "s1", "assessment": {
                           "operator_matches": [{"operator_id": "op", "template_ids": [f"TPL:{i:064x}" for i in range(500)]}],
                           "reaction_signature": {"signature_id": "sig", "events": [{"atom_index": i} for i in range(300)]},
                           "gates": [{"status": "unknown", "summary": "Mapping remains unresolved"}],
                           "precedent_matches": [{"reaction_id": "source", "precursor_smiles": "CC=O.N", "product_smiles": "CCN"}],
                       }}]}},
    }


def raw(store, reference):
    return store.root / "artifacts" / (reference.split(":")[1] + ".json")


def test_compact_record_restores_exact_result_and_deduplicates_between_calls(store):
    original = payload()
    before = deepcopy(original)
    first = store.append("call", original)
    stored = store.read_artifact(first.artifact_ref, expanded=False)
    assert stored["storage_schema_version"] == STORAGE_VERSION
    assert len(canonical_bytes(stored)) < len(canonical_bytes(original)) / 5
    assert original == before
    assert canonical_bytes(store.read_artifact(first.artifact_ref)) == canonical_bytes(original)
    page = ScientificWorkspace(store.root).inspect_artifact(first.artifact_ref, path=[
        "result", "assessment", "step_assessments", 0, "assessment", "operator_matches", 0, "template_ids",
    ], offset=20, limit=2)
    assert page["preview"] == original["result"]["assessment"]["step_assessments"][0]["assessment"]["operator_matches"][0]["template_ids"][20:22]
    assert page["page"]["total"] == 500
    assert summarize_call_brief(store.read_artifact(first.artifact_ref)) == summarize_call_brief(original)
    projected = stored["payload"]["result"]["assessment"]["step_assessments"][0]["assessment"]
    assert projected["gates"] == original["result"]["assessment"]["step_assessments"][0]["assessment"]["gates"]
    assert projected["precedent_matches"][0]["reaction_id"] == "source"
    assert projected["operator_matches"][0]["template_ids"]["item_count"] == 500
    details_before = {p.name for p in (store.root / "artifacts").glob("*.json") if json.loads(p.read_text("utf-8")).get("schema_version") == SECTION_VERSION}
    assert len(details_before) == 2
    second_value = deepcopy(original)
    second_value["arguments"]["label"] = "Other route with the same operator evidence"
    second = store.append("call", second_value)
    assert first.artifact_ref != second.artifact_ref
    assert len(list((store.root / "artifacts").glob("*.json"))) == 4
    reopened = InvestigationStore(store.root)
    assert reopened.read_artifact(second.artifact_ref) == second_value
    replay = store.append("replay", {"source_ref": first.artifact_ref, "matches": True, "result": original["result"]})
    assert store.read_artifact(replay.artifact_ref, expanded=False)["storage_schema_version"] == STORAGE_VERSION
    assert store.read_artifact(replay.artifact_ref)["result"] == original["result"]


@pytest.mark.parametrize("damage", ["missing", "corrupt"])
def test_linked_evidence_failure_rejects_full_compact_and_history_reads(store, damage):
    event = store.append("call", payload())
    stored = store.read_artifact(event.artifact_ref, expanded=False)
    marker = stored["payload"]["result"]["assessment"]["step_assessments"][0]["assessment"]["operator_matches"][0]["template_ids"]
    path = raw(store, marker["artifact_ref"])
    if damage == "missing":
        path.unlink()
    else:
        path.write_text("{}", "utf-8")
    for expanded in (False, True):
        with pytest.raises((FileNotFoundError, ValueError)):
            store.read_artifact(event.artifact_ref, expanded=expanded)
    with pytest.raises((FileNotFoundError, ValueError)):
        store.events()


@pytest.mark.parametrize("change", ["version", "count", "path", "digest"])
def test_inconsistent_envelopes_are_rejected_even_with_valid_outer_hash(store, change):
    event = store.append("call", payload())
    value = store.read_artifact(event.artifact_ref, expanded=False)
    if change == "version":
        value["storage_schema_version"] = "unknown"
    elif change == "count":
        value["payload"]["result"]["assessment"]["step_assessments"][0]["assessment"]["operator_matches"][0]["template_ids"]["item_count"] = 1
    elif change == "path":
        value["section_paths"].append(["result", "absent"])
    else:
        value["expanded_sha256"] = "0" * 64
    reference = "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()
    _write_json(raw(store, reference), value)
    with pytest.raises(ValueError):
        store.read_artifact(reference)


def test_plain_legacy_records_and_non_route_payloads_remain_readable(store):
    value = payload()
    reference = "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()
    _write_json(raw(store, reference), value)
    assert store.read_artifact(reference) == store.read_artifact(reference, expanded=False) == value
    ordinary = store.append("derived_file", value)
    assert ordinary.artifact_ref == reference
    assert json.loads(raw(store, reference).read_text("utf-8")) == value


def test_failed_root_write_does_not_publish_a_partial_event(store, monkeypatch):
    from chem_coworker.scientific_workspace import store as module

    write = module._write_json

    def fail_root(path, value):
        if value.get("storage_schema_version") == STORAGE_VERSION:
            raise OSError("Injected root write failure")
        write(path, value)

    monkeypatch.setattr(module, "_write_json", fail_root)
    with pytest.raises(OSError, match="Injected"):
        store.append("call", payload())
    assert store.events() == ()
    assert not (store.root / ".writer.lock").exists()
    assert all(json.loads(p.read_text("utf-8"))["schema_version"] == SECTION_VERSION
               for p in (store.root / "artifacts").glob("*.json"))


def test_step_inspection_audits_use_literal_keys_and_small_records_stay_plain(store):
    value = {"operation": "inspect_route_step", "execution_status": "completed", "result": {
        "schema_version": "route_step_inspection.v1", "molecule_audits": {
            "C[C@H](N)C(=O)O": {"atoms": [{"index": i, "warnings": []} for i in range(200)]},
        }, "assessment": {"warnings": ["Keep this inline"]},
    }}
    event = store.append("call", value)
    assert store.read_artifact(event.artifact_ref) == value
    assert store.read_artifact(event.artifact_ref, expanded=False)["section_paths"] == [
        ["result", "molecule_audits", "C[C@H](N)C(=O)O"],
    ]
    value["result"]["molecule_audits"] = {}
    small = store.append("call", value)
    assert store.read_artifact(small.artifact_ref, expanded=False) == value


def test_web_artifact_defaults_to_compact_and_offers_explicit_full_expansion(store):
    event = store.append("call", payload())
    service = ConversationService(store.root.parent, runtime=object())
    try:
        with TestClient(create_app(runtime=object(), scientific_service=service, recommendation_only=False), base_url="http://127.0.0.1") as client:
            url = f"/api/v1/scientific/conversations/{store.root.name}/artifacts/{event.artifact_ref}"
            response = client.get(url)
            assert response.status_code == 200
            assert response.json()["storage_schema_version"] == STORAGE_VERSION
            assert client.get(url + "?expanded=true").json() == payload()
            marker = response.json()["payload"]["result"]["assessment"]["step_assessments"][0]["assessment"]["operator_matches"][0]["template_ids"]
            section_url = url.rsplit("/", 1)[0] + "/" + marker["artifact_ref"]
            assert len(client.get(section_url).json()["value"]) == 500
    finally:
        service.close()
