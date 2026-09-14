"""Workbench context route composes planners without affecting recommendations."""

from fastapi.testclient import TestClient

from app.web_api.contracts import ReactionContextRequest
from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime


def test_context_route_validates_request_and_preserves_result():
    class Runtime:
        def reaction_context(self, request):
            assert request.library_mode == "compact"
            return {
                "advisory_only": True,
                "query_reaction_smiles": request.reaction_smiles,
            }

    client = TestClient(create_app(runtime=Runtime(), recommendation_only=False))
    response = client.post(
        "/api/v1/recommendations/context",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "library_mode": "compact",
        },
    )
    assert response.status_code == 200
    assert response.json()["data"]["advisory_only"]
    assert (
        client.post(
            "/api/v1/recommendations/context",
            json={
                "reaction_smiles": "CCBr.N>>CCN",
                "library_mode": "invalid",
            },
        ).status_code
        == 422
    )


def test_context_route_reports_invalid_structures():
    client = TestClient(create_app(recommendation_only=False))
    response = client.post(
        "/api/v1/recommendations/context",
        json={
            "reaction_smiles": "invalid",
        },
    )
    assert response.status_code == 422


def test_missing_artifacts_return_partial_status_without_building(
    tmp_path, monkeypatch
):
    runtime = LocalRecommendationRuntime(retrosynthesis_library_root=tmp_path)

    def must_not_build(*args, **kwargs):
        raise AssertionError("must not rebuild in an interactive context request")

    monkeypatch.setattr(runtime, "_get_forward_library", must_not_build)
    result = runtime.reaction_context(
        ReactionContextRequest(reaction_smiles="CCBr.N>>CCN")
    )
    assert result["reactant_analysis"]["status"] == "unavailable"
    assert result["product_analysis"]["status"] == "unavailable"


def test_focused_deployment_keeps_research_context_private():
    client = TestClient(create_app(recommendation_only=True))
    assert (
        client.post(
            "/api/v1/recommendations/context",
            json={
                "reaction_smiles": "CCBr.N>>CCN",
            },
        ).status_code
        == 404
    )
