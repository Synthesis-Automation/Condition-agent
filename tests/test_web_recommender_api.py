"""Contract tests for the local web application boundary."""

from __future__ import annotations

from typing import Any, Dict

from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime


class FakeRuntime:
    def capabilities(self) -> Dict[str, Any]:
        return {
            "index_available": True,
            "rxnmapper_available": False,
            "featurization": True,
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
            "library_mode": request.library_mode,
            "valid": True,
            "recommendations": [
                {
                    "rank": 1,
                    "recipe_id": "recipe:1",
                    "score": 0.9,
                    "synthesis_protocol": {
                        "protocol_id": "CP1:test",
                        "materials": [
                            {
                                "category": "condition",
                                "canonical_name": "Potassium carbonate",
                                "cas": "584-08-7",
                            }
                        ],
                        "execution_readiness": "review_required",
                    },
                }
            ],
            "schema_version": "test",
        }

    def discover(self, request: Any) -> Dict[str, Any]:
        return {
            "query_reaction_smiles": request.reaction_smiles,
            "library_mode": request.library_mode,
            "valid": True,
            "hits": [{"rank": 1, "reaction_id": "rxn:1"}],
            "schema_version": "test",
        }

    def analyze_features(self, request: Any) -> Dict[str, Any]:
        input_kind = "reaction" if ">" in request.input_smiles else "molecule"
        return {
            "input_kind": input_kind,
            "input_smiles": request.input_smiles,
            "valid": True,
            "overview": {"atom_count": 6},
            "motifs": [],
            "reactive_sites": [],
            "reaction_core": None,
            "partners": [],
            "mapping": None,
            "core_graphic_svg": None,
            "core_projection": None,
            "warnings": [],
            "error": None,
            "analysis": {"schema_version": "test"},
        }

    def render_reaction(
        self, reaction_smiles: str, *, width: int, height: int
    ) -> bytes:
        assert reaction_smiles
        assert width == 760
        assert height == 220
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"

    def render_molecule(
        self, molecule_smiles: str, *, width: int, height: int
    ) -> bytes:
        assert molecule_smiles
        assert width == 760
        assert height == 220
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"


def client() -> TestClient:
    return TestClient(create_app(runtime=FakeRuntime(), frontend_dist="missing"))


def test_local_runtime_reports_isolated_full_and_compact_indexes(tmp_path) -> None:
    for mode in ("full", "compact"):
        mode_dir = tmp_path / mode
        mode_dir.mkdir()
        (mode_dir / "generic_index.sqlite").touch()

    capabilities = LocalRecommendationRuntime(
        library_root=tmp_path
    ).capabilities()

    assert capabilities["default_library_mode"] == "full"
    assert capabilities["library_modes"]["full"]["index_available"] is True
    assert capabilities["library_modes"]["compact"]["index_available"] is True


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
            "library_mode": "compact",
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
    assert payload["data"]["library_mode"] == "compact"
    assert payload["data"]["recommendations"][0]["recipe_id"] == "recipe:1"
    assert (
        payload["data"]["recommendations"][0]["synthesis_protocol"]["materials"][0]["cas"]
        == "584-08-7"
    )


def test_discovery_and_svg_rendering_contracts() -> None:
    web = client()
    discovery = web.post(
        "/api/v1/discovery",
        json={
            "reaction_smiles": "CCBr.N>>CCN",
            "library_mode": "compact",
            "view": "closest_chemistry",
            "use_rxnmapper": False,
        },
    )
    drawing = web.post(
        "/api/v1/render/reaction",
        json={"reaction_smiles": "CCBr.N>>CCN"},
    )

    assert discovery.status_code == 200
    assert discovery.json()["data"]["library_mode"] == "compact"
    assert discovery.json()["data"]["hits"][0]["reaction_id"] == "rxn:1"
    assert drawing.status_code == 200
    assert drawing.headers["content-type"].startswith("image/svg+xml")
    assert drawing.content.startswith(b"<svg")


def test_feature_analysis_and_molecule_rendering_contracts() -> None:
    web = client()
    feature_result = web.post(
        "/api/v1/features/analyze",
        json={"input_smiles": "Brc1ccccc1", "use_rxnmapper": False},
    )
    drawing = web.post(
        "/api/v1/render/molecule",
        json={"molecule_smiles": "Brc1ccccc1"},
    )

    assert feature_result.status_code == 200
    assert feature_result.json()["data"]["input_kind"] == "molecule"
    assert feature_result.json()["data"]["overview"]["atom_count"] == 6
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
