"""Web fragment search uses the standalone index and retains evidence semantics."""

from dataclasses import asdict
import json

from fastapi.testclient import TestClient
import pytest

from app.web_api.contracts import FragmentSearchRequest
from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime
from condition_recommender.fragment_index import build_fragment_index
from reactive_taxonomy import featurize_reaction


@pytest.fixture
def runtime(tmp_path):
    reaction = "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]"
    source = tmp_path / "records.jsonl"
    source.write_text(json.dumps({
        "observation_id": "obs-1", "reaction_id": "rxn-1", "reference_id": "ref-1",
        "reaction_smiles": reaction, "admission_tier": "review",
        "reaction_observation": asdict(featurize_reaction(reaction).observation),
    }) + "\n", encoding="utf-8")
    index = tmp_path / "fragment.sqlite"
    build_fragment_index(source, index)
    return LocalRecommendationRuntime(fragment_index_path=index)


def test_workbench_search_and_empty_result(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.get("/api/v1/capabilities").json()["data"]["fragment_search"]
    response = client.post("/api/v1/fragments/search", json={"query": "COC"})
    assert response.status_code == 200
    data = response.json()["data"]
    assert data["counts"]["observations"] == {"value": 1, "precision": "exact"}
    assert "constructed" in data["hits"][0]["relationships"]
    assert data["hits"][0]["reference_id"] == "ref-1"
    assert data["hits"][0]["admission_tier"] == "review"
    empty = client.post("/api/v1/fragments/search", json={"query": "P(=O)(O)O"})
    assert empty.json()["data"]["search_status"] == "complete"
    assert empty.json()["data"]["hits"] == []


def test_target_only_discovery_and_request_bounds(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post("/api/v1/fragments/discover", json={"target_smiles": "COC"})
    assert response.status_code == 200
    assert response.json()["data"]["schema_version"] == "synthesis_precedent_discovery.v1"
    assert response.json()["data"]["hits"]
    assert response.json()["data"]["execution"]["library_loads"] == 1
    assert client.post("/api/v1/fragments/discover", json={"target_smiles": "bad"}).status_code == 422
    assert client.post("/api/v1/fragments/discover", json={"target_smiles": "CO", "timeout_seconds": 121}).status_code == 422
    focused = TestClient(create_app(runtime=runtime))
    assert focused.post("/api/v1/fragments/discover", json={"target_smiles": "COC"}).status_code == 404


def test_workspace_target_only_discovery_records_evidence(runtime, tmp_path):
    from chem_coworker.scientific_workspace import ScientificWorkspace
    from pathlib import Path

    workspace = ScientificWorkspace.create(tmp_path / "investigation", objective="Discovery fixture",
                                           repository=Path.cwd(), artifacts={"fragment_index": runtime.fragment_index_path})
    event = workspace.run("find_synthesis_precedents", {"target_smiles": "COC", "timeout_seconds": 30})
    payload = workspace.store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "completed"
    assert payload["result"]["hits"][0]["observation_id"] == "obs-1"
    assert payload["result"]["execution"]["diagnostics"]
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]


@pytest.mark.parametrize("payload", [
    {"query": "bad smiles"}, {"query": ""}, {"query": "CO", "limit": 11},
    {"query": "CO", "timeout_seconds": 31}, {"query": "CO", "query_format": "fuzzy"},
    {"query": "CO", "index_path": "arbitrary.sqlite"},
])
def test_invalid_queries_are_422(runtime, payload):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.post("/api/v1/fragments/search", json=payload).status_code == 422


@pytest.mark.parametrize("field,value,label", [
    ("query", "bad smiles", "Core fragment (query)"),
    ("target_smiles", "n1c2c(cncc2)c2C(CCCc12)=O", "Target molecule (target_smiles)"),
])
def test_search_errors_identify_field_in_response_and_log(runtime, caplog, field, value, label):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    payload = {"query": "CO", field: value}
    response = client.post("/api/v1/fragments/search", json=payload)
    assert response.status_code == 422
    assert label in response.json()["detail"]["message"]
    assert label in caplog.text
    assert "/api/v1/fragments/search -> 422" in caplog.text


def test_search_request_validation_preserves_details_and_logs_field(runtime, caplog):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post("/api/v1/fragments/search", json={"query": "CO", "timeout_seconds": 31})
    assert response.status_code == 422
    assert response.json()["detail"][0]["loc"] == ["body", "timeout_seconds"]
    assert "body.timeout_seconds" in caplog.text
    assert "less than or equal to 30" in caplog.text


def test_drawn_explicit_hydrogen_core_is_accepted(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post("/api/v1/fragments/search", json={
        "query": "[n]1([H])c2c(C(=O)CCC2)c2cnccc12",
        "target_smiles": "CC1(C)CC(C)(C)c2[nH]c3ccncc3c2C1=O",
    })
    assert response.status_code == 200
    assert response.json()["data"]["search_status"] == "complete"


def test_explicit_smarts_and_partial_status_are_preserved(runtime, monkeypatch):
    import condition_recommender.fragment_search as search

    original = search.search_fragment_precedents
    def partial(*args, **kwargs):
        assert kwargs["query_format"] == "smarts"
        assert kwargs["topology"] == "subgraph"
        result = original(*args, **kwargs)
        result.update(search_status="partial", stop_reason="deadline")
        result["counts"]["observations"]["precision"] = "at_least"
        return result
    monkeypatch.setattr(search, "search_fragment_precedents", partial)
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post("/api/v1/fragments/search", json={
        "query": "C[O,N]C", "query_format": "smarts", "topology": "subgraph",
    })
    assert response.status_code == 200
    data = response.json()["data"]
    assert data["search_status"] == "partial" and data["hits"]
    assert data["counts"]["observations"]["precision"] == "at_least"


def test_missing_index_and_concurrent_load_are_explicit(runtime, tmp_path):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    runtime._fragment_search_lock.acquire()
    try:
        response = client.post("/api/v1/fragments/search", json={"query": "CO"})
        assert response.status_code == 503
        assert "already running" in response.json()["detail"]["message"]
    finally:
        runtime._fragment_search_lock.release()
    runtime.fragment_index_path = tmp_path / "missing.sqlite"
    assert client.post("/api/v1/fragments/search", json={"query": "CO"}).status_code == 503
    assert not client.get("/api/v1/capabilities").json()["data"]["fragment_search"]


def test_failed_query_releases_search_lock(runtime):
    with pytest.raises(ValueError):
        runtime.search_fragments(FragmentSearchRequest(query="bad smiles"))
    assert runtime.search_fragments(FragmentSearchRequest(query="COC"))["hits"]


def test_focused_profile_stays_recommendation_only(runtime):
    client = TestClient(create_app(runtime=runtime))
    assert client.post("/api/v1/fragments/search", json={"query": "CO"}).status_code == 404
    assert "fragment_search" not in client.get("/api/v1/capabilities").json()["data"]


def test_fragment_index_environment_override(tmp_path, monkeypatch):
    path = tmp_path / "configured.sqlite"
    monkeypatch.setenv("FRAGMENT_PRECEDENT_INDEX", str(path))
    assert LocalRecommendationRuntime().fragment_index_path == path


def test_suggestions_work_without_index_and_then_feed_search(runtime, monkeypatch, tmp_path):
    import condition_recommender.fragment_search as search

    def forbidden(*args, **kwargs):
        pytest.fail("Suggestion must not search the corpus")

    monkeypatch.setattr(search, "search_fragment_precedents", forbidden)
    runtime.fragment_index_path = tmp_path / "missing.sqlite"
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.get("/api/v1/capabilities").json()["data"]["fragment_suggestions"]
    response = client.post("/api/v1/fragments/suggest", json={
        "target_smiles": "CC(=O)c1ccc2c(c1)COc1ccccc1-2", "limit": 2,
    })
    assert response.status_code == 200
    result = response.json()["data"]
    candidate = result["candidates"][0]
    assert candidate["query"] == "c1ccc2c(c1)COc1ccccc1-2"
    assert "<svg" in candidate["target_highlight_svg"]
    assert "ellipse" in candidate["target_highlight_svg"]
    selected = client.post("/api/v1/fragments/suggest", json={
        "target_smiles": result["target_smiles"], "selected_atom_ids": candidate["target_atom_ids"],
    })
    assert selected.json()["data"]["candidates"][0]["query"] == candidate["query"]


@pytest.mark.parametrize("payload", [
    {"target_smiles": "C.O"}, {"target_smiles": "bad"}, {"target_smiles": "CO", "limit": True},
    {"target_smiles": "CO", "selected_atom_ids": [True]},
    {"target_smiles": "CO", "selected_atom_ids": [0.5]},
    {"target_smiles": "CO", "selected_atom_ids": ["0"]},
    {"target_smiles": "CO", "selected_atom_ids": [90]},
    {"target_smiles": "CO", "auto_search": True},
])
def test_invalid_suggestion_requests_are_422(runtime, payload):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.post("/api/v1/fragments/suggest", json=payload).status_code == 422


def test_suggestion_route_respects_focused_profile(runtime):
    client = TestClient(create_app(runtime=runtime))
    assert client.post("/api/v1/fragments/suggest", json={"target_smiles": "CO"}).status_code == 404
    assert "fragment_suggestions" not in client.get("/api/v1/capabilities").json()["data"]
