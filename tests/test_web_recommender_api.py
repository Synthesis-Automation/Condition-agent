"""Contract tests for the local web application boundary."""

from __future__ import annotations

from typing import Any, Dict

from fastapi.testclient import TestClient

from app.web_api.main import create_app


class FakeRuntime:
    def capabilities(self) -> Dict[str, Any]:
        return {
            "index_available": True,
            "rxnmapper_available": False,
            "local_only": True,
        }

    def ranking_profiles(self) -> tuple[Dict[str, Any], ...]:
        return (
            {
                "profile_id": "default",
                "label": "Balanced",
                "description": "Balanced evidence",
                "weights": {"similarity": 1.0},
            },
        )

    def prepare_reaction(self, reaction_smiles: str) -> Dict[str, Any]:
        if reaction_smiles == "invalid":
            raise ValueError("INVALID_REACTION")
        return {
            "valid": True,
            "input_reaction_smiles": reaction_smiles,
            "completion_proposal": {
                "status": "not_required",
                "requirements": [],
            },
            "warnings": [],
        }

    def recommend(self, request: Any) -> Dict[str, Any]:
        return {
            "query_reaction_smiles": request.reaction_smiles,
            "valid": True,
            "recommendations": [
                {
                    "rank": 1,
                    "recipe_id": "recipe:1",
                    "score": 0.9,
                }
            ],
            "schema_version": "test",
        }

    def discover(self, request: Any) -> Dict[str, Any]:
        return {
            "query_reaction_smiles": request.reaction_smiles,
            "valid": True,
            "hits": [{"rank": 1, "reaction_id": "rxn:1"}],
            "schema_version": "test",
        }

    def render_reaction(
        self, reaction_smiles: str, *, width: int, height: int
    ) -> bytes:
        assert reaction_smiles
        assert width == 760
        assert height == 220
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"


def client() -> TestClient:
    return TestClient(create_app(runtime=FakeRuntime(), frontend_dist="missing"))


def test_health_and_capabilities_are_versioned() -> None:
    web = client()
    health = web.get("/api/v1/health")
    capabilities = web.get("/api/v1/capabilities")

    assert health.status_code == 200
    assert health.json()["api_schema_version"] == "1.0"
    assert capabilities.json()["data"]["index_available"] is True


def test_prepare_reaction_returns_typed_error() -> None:
    response = client().post(
        "/api/v1/reactions/prepare",
        json={"reaction_smiles": "invalid"},
    )

    assert response.status_code == 422
    assert response.json()["detail"] == {
        "code": "INVALID_REACTION",
        "message": "INVALID_REACTION",
    }


def test_recommendation_contract_forwards_validated_options() -> None:
    response = client().post(
        "/api/v1/recommendations",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "top_k": 4,
            "use_rxnmapper": False,
            "ranking_preferences": {
                "profile_id": "default",
                "weights": {},
            },
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["api_schema_version"] == "1.0"
    assert payload["data"]["query_reaction_smiles"] == "CCBr.N>>CCN"
    assert payload["data"]["recommendations"][0]["recipe_id"] == "recipe:1"


def test_discovery_and_svg_rendering_contracts() -> None:
    web = client()
    discovery = web.post(
        "/api/v1/discovery",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "view": "closest_chemistry",
            "use_rxnmapper": False,
        },
    )
    drawing = web.post(
        "/api/v1/render/reaction",
        json={"reaction_smiles": "CCBr.N>>CCN"},
    )

    assert discovery.status_code == 200
    assert discovery.json()["data"]["hits"][0]["reaction_id"] == "rxn:1"
    assert drawing.status_code == 200
    assert drawing.headers["content-type"].startswith("image/svg+xml")
    assert drawing.content.startswith(b"<svg")


def test_request_contract_rejects_unknown_fields() -> None:
    response = client().post(
        "/api/v1/recommendations",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "unknown_option": True,
        },
    )

    assert response.status_code == 422

