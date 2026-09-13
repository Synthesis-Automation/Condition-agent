"""Focused HTTP surface and partial-failure integration contracts."""

from fastapi.testclient import TestClient

from app.web_api.main import create_app


REACTION = "CCBr.N>>CCN"


class ConditionsRuntime:
    def __init__(self):
        self.requests = []
        self.fail_generic = False

    def capabilities(self):
        return {
            "recommendation": True,
            "weak_label_recommendation": True,
            "rxnmapper_available": False,
            "retrosynthesis": True,
        }

    def recommend(self, request):
        self.requests.append(request)
        if self.fail_generic and request.recommendation_mode == "generic":
            raise FileNotFoundError("/private/path/index.sqlite")
        return {
            "valid": True,
            "query_reaction_smiles": request.reaction_smiles,
            "recommendation_mode": "verified_signature"
            if request.recommendation_mode == "generic"
            else request.recommendation_mode,
            "recommendations": [
                {
                    "resolved_recipe": {
                        "recipe_id": "RCR2:test",
                        "temperature_c": 80,
                        "time_h": 1,
                    },
                    "cautions": [],
                    "support": 1,
                }
            ],
            "warnings": [],
        }


def test_combined_api_queries_both_sources_and_exports_each_recipe(tmp_path):
    runtime = ConditionsRuntime()
    client = TestClient(create_app(runtime=runtime, frontend_dist=tmp_path))
    response = client.post(
        "/api/v1/conditions/recommend", json={"reaction_smiles": REACTION}
    )
    assert response.status_code == 200
    data = response.json()["data"]
    assert [request.recommendation_mode for request in runtime.requests] == [
        "generic",
        "weak_label_fallback",
    ]
    assert not any(request.use_rxnmapper for request in runtime.requests)
    assert len(data["recommendations"]) == 1
    option = data["recommendations"][0]
    export = data["automation_exports"][option["option_id"]]
    assert export["execution_ready"] is False
    assert export["experiments"][0]["protocol"]["reaction_smiles"] == REACTION


def test_focused_deployment_has_no_research_routes(tmp_path):
    client = TestClient(
        create_app(
            runtime=ConditionsRuntime(),
            frontend_dist=tmp_path / "absent",
            recommendation_only=True,
        )
    )
    paths = client.get("/api/openapi.json").json()["paths"]
    assert "/api/v1/conditions/recommend" in paths
    assert "/api/v1/render/reaction" in paths
    assert "/api/v1/reactions/prepare" in paths
    assert "retrosynthesis" not in client.get("/api/v1/capabilities").json()["data"]
    for path in (
        "/api/v1/retrosynthesis",
        "/api/v1/forward-synthesis",
        "/api/v1/features/analyze",
        "/api/v1/experimental/assistance",
        "/api/v1/recommendations",
    ):
        assert path not in paths
        assert client.post(path, json={}).status_code == 404


def test_engine_failure_is_reported_without_losing_weak_label_options(tmp_path):
    runtime = ConditionsRuntime()
    runtime.fail_generic = True
    client = TestClient(create_app(runtime=runtime, frontend_dist=tmp_path))
    response = client.post(
        "/api/v1/conditions/recommend", json={"reaction_smiles": REACTION}
    )
    data = response.json()["data"]
    assert data["valid"]
    assert data["sources"][0]["status"] == "unavailable"
    assert data["recommendations"][0]["evidence_kind"] == "weak_label"
    assert "/private/path" not in response.text


def test_completion_choices_are_not_silently_dropped_for_weak_labels(tmp_path):
    runtime = ConditionsRuntime()
    client = TestClient(create_app(runtime=runtime, frontend_dist=tmp_path))
    response = client.post(
        "/api/v1/conditions/recommend",
        json={
            "reaction_smiles": REACTION,
            "completion_choices": [
                {"requirement_id": "fragment", "option_id": "source"}
            ],
        },
    )
    assert response.json()["data"]["sources"][1]["status"] == "skipped"
    assert len(runtime.requests) == 1
    assert runtime.requests[0].completion_choices[0].option_id == "source"


def test_api_rejects_empty_queries_and_user_supplied_library_paths(tmp_path):
    client = TestClient(create_app(runtime=ConditionsRuntime(), frontend_dist=tmp_path))
    assert (
        client.post(
            "/api/v1/conditions/recommend", json={"reaction_smiles": ""}
        ).status_code
        == 422
    )
    assert (
        client.post(
            "/api/v1/conditions/recommend", json={"reaction_smiles": "   "}
        ).status_code
        == 422
    )
    assert (
        client.post(
            "/api/v1/conditions/recommend",
            json={"reaction_smiles": REACTION, "library_path": "/tmp/other"},
        ).status_code
        == 422
    )


def test_outdated_index_reports_rebuild_requirement_without_private_paths(tmp_path):
    class OutdatedRuntime(ConditionsRuntime):
        def recommend(self, request):
            if request.recommendation_mode == "generic":
                raise ValueError(
                    "Unsupported index schema; rebuild /private/index.sqlite"
                )
            return super().recommend(request)

    client = TestClient(create_app(runtime=OutdatedRuntime(), frontend_dist=tmp_path))
    response = client.post(
        "/api/v1/conditions/recommend", json={"reaction_smiles": REACTION}
    )
    source = response.json()["data"]["sources"][0]
    assert source["error_code"] == "INDEX_REBUILD_REQUIRED"
    assert "rebuilt" in source["message"]
    assert "/private" not in response.text


def test_current_structure_index_produces_a_real_protocol_through_focused_api(
    tmp_path, monkeypatch
):
    from app.web_api.runtime import LocalRecommendationRuntime
    from condition_recommender.conversion.generic import convert_record
    from condition_recommender.conversion.input_schema import adapt_row
    from condition_recommender.generic_indexing import build_generic_index, load_generic_index
    from condition_recommender.shared_core_index import build_shared_core_index
    from condition_recommender.sqlite_indexing import save_sqlite_generic_index

    reaction = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    records = []
    for number in range(6):
        raw = adapt_row(
            {
                "reaction_id": f"test-{number}",
                "reaction_smiles": reaction,
                "yield_pct": "80",
                "catalyst_cas": "14221-01-3",
                "reagent_cas": "584-08-7",
                "solvent_cas": "108-88-3",
                "reference": f"TEST FIXTURE reference {number}",
                "stages": "1",
            },
            source_dataset="TEST FIXTURE",
            source_path="test.csv",
            source_row_number=number + 2,
        )
        records.append(convert_record(raw).to_dict())
    index_path = tmp_path / "generic_index.sqlite"
    save_sqlite_generic_index(build_generic_index(records), index_path)
    build_shared_core_index(
        load_generic_index(index_path), index_path.with_suffix(".shared_core.sqlite")
    )
    monkeypatch.delenv("CONDITION_SHARED_CORE_EXPERIMENTAL", raising=False)
    runtime = LocalRecommendationRuntime(
        index_path, weak_label_records_path=tmp_path / "absent.csv"
    )
    monkeypatch.setenv("CONDITION_RECOMMENDER_USE_RXNMAPPER", "false")
    client = TestClient(
        create_app(
            runtime=runtime, recommendation_only=True, frontend_dist=tmp_path / "absent"
        )
    )
    response = client.post(
        "/api/v1/conditions/recommend", json={"reaction_smiles": reaction}
    )
    assert response.status_code == 200
    result = response.json()["data"]
    assert result["sources"][0]["status"] == "ok", result["sources"]
    assert result["sources"][0]["result"]["recommendation_mode"] == "experimental_shared_core"
    # Selecting v2 by default does not erase its pending-review provenance.
    assert result["recommendations"][0]["evidence_kind"] == "structure_review"
    protocol = result["recommendations"][0]["synthesis_protocol"]
    assert (
        sum(item["category"] == "reaction_input" for item in protocol["materials"]) == 2
    )
    assert any(item["cas"] == "584-08-7" for item in protocol["materials"])
