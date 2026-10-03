"""Workbench composite endpoint uses the shared catalogue and chemistry gates."""

import json

from fastapi.testclient import TestClient

from app.web_api.main import create_app
from app.web_api.runtime import LocalRecommendationRuntime
from core_retrosynthesis.cli import main
from core_retrosynthesis.composite_actions import (
    build_composite_strategy_catalog,
    save_composite_strategy_catalog,
)
from core_retrosynthesis.generic_library import save_generic_library
from tests.core_retrosynthesis_tests import test_composite_actions as composite_fixtures
from tests.core_retrosynthesis_tests.test_composite_actions import (
    ACTIVATION,
    SUBSTITUTION,
    candidate,
    strategy,
)


# Explicitly expose the shared chemistry fixture to pytest.
composite_library = composite_fixtures.composite_library


def test_real_workbench_and_cli_share_composite_actions(
    tmp_path, composite_library, capsys
):
    library = tmp_path / "operators.json.gz"
    catalog_path = tmp_path / "catalog.json"
    save_generic_library(composite_library, library)
    catalog = build_composite_strategy_catalog(
        (strategy(candidate(ACTIVATION), candidate(SUBSTITUTION)),)
    )
    save_composite_strategy_catalog(catalog, catalog_path)
    runtime = LocalRecommendationRuntime(
        coupled_strategy_library_path=library,
        coupled_strategy_catalog_path=catalog_path,
    )
    with TestClient(create_app(runtime=runtime, recommendation_only=False)) as client:
        response = client.post(
            "/api/v1/retrosynthesis/coupled-strategies",
            json={"target_smiles": "CCN", "top_k": 1},
        )
        assert response.status_code == 200, response.text
        result = response.json()["data"]
        assert result["catalog_id"] == catalog.catalog_id and "panel_id" not in result
        assert result["schema_version"] == "1.1"
        action = result["actions"][0]
        assert (
            action["physical_step_cost"] == 2
            and action["dependency"]["status"] == "verified"
        )
        assert len(action["physical_steps"]) == 2
        assert client.get("/api/openapi.json").status_code == 200
    assert (
        main(
            [
                "disconnect-composite",
                str(library),
                str(catalog_path),
                "NCC",
                "--top-k",
                "1",
            ]
        )
        == 0
    )
    cli = json.loads(capsys.readouterr().out)
    assert cli["catalog_id"] == catalog.catalog_id
    assert cli["actions"] == result["actions"]


def test_evaluation_panel_is_not_a_runtime_catalogue(tmp_path, composite_library):
    library = tmp_path / "operators.json.gz"
    panel = tmp_path / "panel.json"
    save_generic_library(composite_library, library)
    panel.write_text(
        '{"artifact_type":"v1_coupled_strategy_frozen_panel","strategies":[]}',
        encoding="utf-8",
    )
    runtime = LocalRecommendationRuntime(
        coupled_strategy_library_path=library, coupled_strategy_catalog_path=panel
    )
    with TestClient(create_app(runtime=runtime, recommendation_only=False)) as client:
        response = client.post(
            "/api/v1/retrosynthesis/coupled-strategies", json={"target_smiles": "CCN"}
        )
        assert response.status_code == 422
        assert "export" in response.text
