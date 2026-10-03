"""Research API preserves discovery uncertainty and canonical chemistry gates."""

from dataclasses import asdict
import json

from fastapi.testclient import TestClient
import pytest

from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime
from condition_recommender.fragment_index import build_fragment_index
from core_retrosynthesis import build_generic_library, save_generic_library
from core_retrosynthesis.fragment_guidance import evaluate_transfers
from reactive_taxonomy import featurize_reaction


ENDPOINT = "/api/v1/retrosynthesis/fragment-guided"
TRANSFER = "/api/v1/retrosynthesis/fragment-transfer"
RING_REACTION = (
    "[CH2:1]=[CH2:2].[CH2:3]=[CH2:4]>>[CH2:1]1[CH2:2][CH2:3][CH2:4]1"
)


def make_runtime(tmp_path, reaction=RING_REACTION, library_reaction=None):
    """Prepare discovery and executable libraries through their canonical builders."""
    row = {
        "observation_id": "obs-1", "reaction_id": "source-1", "reference_id": "ref-1",
        "reaction_smiles": reaction, "admission_tier": "review",
        "reaction_observation": asdict(featurize_reaction(reaction).observation),
    }
    source, index = tmp_path / "source.jsonl", tmp_path / "fragments.sqlite"
    source.write_text(json.dumps(row) + "\n", encoding="utf-8")
    build_fragment_index(source, index)
    library_row = {**row, "reaction_smiles": library_reaction or reaction}
    library = build_generic_library([library_row], admission_mode="data_driven")
    root = tmp_path / "operators"
    path = root / "compact" / "operator_library_v3.json.gz"
    path.parent.mkdir(parents=True)
    save_generic_library(library, path)
    return LocalRecommendationRuntime(
        fragment_index_path=index, retrosynthesis_library_root=root,
    )


@pytest.fixture
def runtime(tmp_path):
    return make_runtime(tmp_path)


def test_single_step_endpoint_uses_shared_core_and_does_not_write_library(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.get("/api/v1/capabilities").json()["data"]["fragment_guided_retrosynthesis"]
    response = client.post(ENDPOINT, json={"target_smiles": "C1CCC1"})
    assert response.status_code == 200
    data = response.json()["data"]
    assert data["schema_version"] == "fragment_guided_retrosynthesis.v1"
    assert data["guidance"]["focus_bonds"]
    assert data["suggestions"]["candidates"][0]["target_highlight_svg"].startswith("<?xml")
    assert data["searches"][0]["result"]["target_validation"]["matches_target"]
    shared = evaluate_transfers(
        {"guidance": data["guidance"], "policy": data["policy"]},
        library=runtime._get_retrosynthesis_library("compact"), repeat=False,
        source_library_path=None,
    )
    assert data["transfers"] == json.loads(json.dumps(shared))
    assert shared["repeat_scientific_results_identical"] is None
    assert shared["source_library_file"] is None
    assert not (tmp_path / "selected_source_operators.json").exists()
    proposals = shared["guided"][0]["direct_source_transfer"]["candidates"]
    assert proposals and all(c["forward_validation_status"] == "verified_signature" for c in proposals)
    assert all(c["bond_focus_check"]["status"] == "verified" for c in proposals)
    assert data["resources"]["fragment_index"]["size_bytes"] > 0
    assert not runtime._fragment_search_lock.locked()


def test_local_construction_does_not_override_whole_source_rejection(tmp_path):
    reaction = "[CH3:1][CH2:2][Br:5].[NH2:3][CH3:4]>>[CH3:1][CH2:2][NH:3][CH3:4]"
    runtime = make_runtime(tmp_path, reaction, reaction.replace("[Br:5]", "Br"))
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    data = client.post(ENDPOINT, json={"target_smiles": "CCNC"}).json()["data"]
    assert data["guidance"]["focus_bonds"]
    assert data["transfers"]["compiled_source_template_count"] == 0
    assert data["transfers"]["source_admissions"][0]["reason"] == "materialized_core_not_verified"
    assert data["transfers"]["guided"][0]["direct_source_transfer"]["candidates"] == []
    assert data["transfers"]["guided"][0]["witness_directed_library"]["candidates"]


@pytest.mark.parametrize("status", ["partial", "too_broad"])
def test_incomplete_search_preserves_results_but_cannot_seed_guidance(runtime, monkeypatch, status):
    import condition_recommender.fragment_search as search
    original = search.search_fragment_precedents

    def incomplete(*args, **kwargs):
        result = original(*args, **kwargs)
        result.update(search_status=status, stop_reason="deadline")
        return result

    monkeypatch.setattr(search, "search_fragment_precedents", incomplete)
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    data = client.post(ENDPOINT, json={"target_smiles": "C1CCC1"}).json()["data"]
    assert data["searches"][0]["result"]["search_status"] == status
    assert data["searches"][0]["result"]["hits"]
    assert not data["guidance"]["focus_bonds"]
    assert data["guidance"]["exclusions"][0]["reason"] == "incomplete_or_unvalidated_search"
    assert data["transfers"]["baseline"]["candidates"]


def test_no_matches_is_complete_empty_guidance(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post(ENDPOINT, json={"target_smiles": "P(=O)(O)(O)O"})
    assert response.status_code == 200
    data = response.json()["data"]
    assert all(s["result"]["search_status"] == "complete" for s in data["searches"])
    assert not data["guidance"]["focus_bonds"]
    assert not data["transfers"]["guided"]
    assert data["transfers"]["comparison"]["additional_guided_precursor_sets"] == []


@pytest.mark.parametrize("payload", [
    {"target_smiles": "bad smiles"}, {"target_smiles": ""},
    {"target_smiles": "C1CCC1", "query_limit": 6},
    {"target_smiles": "C1CCC1", "query_limit": True},
    {"target_smiles": "C1CCC1", "max_focus_bonds": 0},
    {"target_smiles": "C1CCC1", "top_k": 11},
    {"target_smiles": "C1CCC1", "library_mode": "arbitrary"},
    {"target_smiles": "C1CCC1", "require_complete_search": False},
    {"target_smiles": "C1CCC1", "index_path": "arbitrary.sqlite"},
])
def test_strict_bounded_requests_and_lock_release(runtime, payload):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.post(ENDPOINT, json=payload).status_code == 422
    assert not runtime._fragment_search_lock.locked()


def test_concurrent_fragment_search_rejected(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    runtime._fragment_search_lock.acquire()
    try:
        assert client.post(ENDPOINT, json={"target_smiles": "C1CCC1"}).status_code == 503
    finally:
        runtime._fragment_search_lock.release()


def test_missing_library_and_index_are_unavailable(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.post(ENDPOINT, json={"target_smiles": "C1CCC1", "library_mode": "full"}).status_code == 503
    assert not runtime._fragment_search_lock.locked()
    runtime.fragment_index_path.unlink()
    assert not client.get("/api/v1/capabilities").json()["data"]["fragment_guided_retrosynthesis"]
    assert client.post(ENDPOINT, json={"target_smiles": "C1CCC1"}).status_code == 503


def test_search_error_releases_lock(runtime, monkeypatch):
    import condition_recommender.fragment_search as search

    def fail(*args, **kwargs):
        raise ValueError("Conflicting target evidence")

    monkeypatch.setattr(search, "search_fragment_precedents", fail)
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post(ENDPOINT, json={"target_smiles": "C1CCC1"})
    assert response.status_code == 422
    assert "Conflicting target" in response.text
    assert not runtime._fragment_search_lock.locked()


def test_focused_deployment_hides_research_endpoint(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=True))
    assert client.post(ENDPOINT, json={"target_smiles": "C1CCC1"}).status_code == 404
    assert "fragment_guided_retrosynthesis" not in client.get("/api/v1/capabilities").json()["data"]


def transfer_request(**overrides):
    return {"target_smiles": "C1CCC1", "query": "C1CCC1",
            "selected_observation_ids": ["obs-1"], **overrides}


def test_manual_transfer_selects_indexed_source_and_preserves_target_alignment_alternatives(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post(TRANSFER, json=transfer_request())
    assert response.status_code == 200, response.text
    data = response.json()["data"]
    assert data["guidance"]["selected_observation_ids"] == ["obs-1"]
    assert data["guidance"]["target_alignment_count"] > 1
    assert data["guidance"]["focus_bonds"]
    assert data["transfers"]["baseline"]["status"] == "not_requested"
    assert data["transfers"]["comparison"]["baseline_requested"] is False
    assert data["transfers"]["comparison"]["additional_guided_precursor_sets"] == []
    assert data["source_comparisons"][0]["comparison"]["same_constitution"]
    assert data["query_search"]["target_validation"]["matches_target"]
    assert data["transfers"]["guided"][0]["direct_source_transfer"]["candidates"]
    assert not runtime._fragment_search_lock.locked()


def test_manual_transfer_optional_baseline_is_explicit(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    data = client.post(TRANSFER, json=transfer_request(include_baseline=True)).json()["data"]
    assert data["transfers"]["baseline"]["status"] == "completed"
    assert data["transfers"]["baseline"]["candidates"]
    assert data["transfers"]["comparison"]["baseline_requested"]


def test_selected_source_rejection_does_not_erase_local_observations(tmp_path):
    reaction = "[CH3:1][CH2:2][Br:5].[NH2:3][CH3:4]>>[CH3:1][CH2:2][NH:3][CH3:4]"
    runtime = make_runtime(tmp_path, reaction, reaction.replace("[Br:5]", "Br"))
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    data = client.post(TRANSFER, json=transfer_request(target_smiles="CCNC", query="CCNC")).json()["data"]
    assert data["guidance"]["focus_bonds"]
    assert data["transfers"]["source_admissions"][0]["reason"] == "materialized_core_not_verified"
    assert data["transfers"]["compiled_source_template_count"] == 0
    assert data["transfers"]["guided"][0]["witness_directed_library"]["candidates"]


@pytest.mark.parametrize("changes", [
    {"selected_observation_ids": []}, {"selected_observation_ids": ["missing"]},
    {"selected_observation_ids": ["obs-1", "obs-1"]},
    {"query": "COC"}, {"include_baseline": "false"},
    {"source_records": [{"reaction_smiles": "fake>>fake"}]},
])
def test_manual_transfer_rejects_unbound_queries_and_untrusted_sources(runtime, changes):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    assert client.post(TRANSFER, json=transfer_request(**changes)).status_code == 422
    assert not runtime._fragment_search_lock.locked()


def test_incomplete_selected_search_cannot_seed_transfer(runtime, monkeypatch):
    import condition_recommender.fragment_search as search
    original = search.search_fragment_precedents

    def partial(*args, **kwargs):
        result = original(*args, **kwargs)
        result["search_status"] = "partial"
        return result

    monkeypatch.setattr(search, "search_fragment_precedents", partial)
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    data = client.post(TRANSFER, json=transfer_request()).json()["data"]
    assert data["query_search"]["search_status"] == "partial"
    assert not data["guidance"]["focus_bonds"]
    assert data["guidance"]["exclusions"]


def test_target_bound_fragment_search_checks_membership(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    response = client.post("/api/v1/fragments/search", json={"target_smiles": "C1CCC1", "query": "C1CCC1"})
    assert response.status_code == 200
    assert response.json()["data"]["target_validation"]["matches_target"]
    assert client.post("/api/v1/fragments/search", json={"target_smiles": "C1CCC1", "query": "COC"}).status_code == 422


def test_focused_deployment_hides_transfer_endpoint(runtime):
    client = TestClient(create_app(runtime=runtime, recommendation_only=True))
    assert client.post(TRANSFER, json=transfer_request()).status_code == 404
