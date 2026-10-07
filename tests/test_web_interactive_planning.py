"""HTTP planner tests using real session edits and a bounded fake search engine."""

from fastapi.testclient import TestClient

from app.web_api.main import create_app
from tests.core_retrosynthesis_tests.test_interactive_planning import search_result


class PlanningRuntime:
    """Serve predictable candidate sets through the existing runtime interface."""

    def __init__(self):
        self.requests = []

    def retrosynthesize(self, request):
        self.requests.append(request)
        return search_result(request.target_smiles, "CCBr.O")

    def retrosynthesis_conditions(self, request):
        self.requests.append(request)
        return {
            "status": "insufficient_evidence",
            "recommendations": [],
            "warnings": ["NO_COMPATIBLE_PRECEDENTS"],
        }


def test_api_edits_use_existing_runtime_and_preserve_condition_evidence():
    runtime = PlanningRuntime()
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    url = "/api/v1/retrosynthesis/planner"
    started = client.post(url, json={"action": "start", "target_smiles": "OCC"})
    assert started.status_code == 200
    session = started.json()["data"]["session"]
    found = client.post(url, json={"action": "search", "session": session}).json()[
        "data"
    ]
    assert runtime.requests[-1].target_smiles == "CCO"
    session = found["session"]
    chosen = client.post(
        url,
        json={
            "action": "select",
            "session": session,
            "search_id": session["searches"][0]["search_id"],
        },
    )
    assert chosen.status_code == 200
    data = chosen.json()["data"]
    assert data["summary"]["reaction_count"] == 1
    conditions = client.post(
        url, json={"action": "conditions", "session": data["session"]}
    )
    assert conditions.status_code == 200
    data = conditions.json()["data"]
    assert (
        list(data["session"]["conditions"].values())[0]["status"]
        == "insufficient_evidence"
    )
    assert runtime.requests[-1].starting_materials == "CCBr.O"
    assert runtime.requests[-1].intended_product == "CCO"
    assert (
        client.post(
            url, json={"action": "restore", "session": data["session"]}
        ).status_code
        == 200
    )


def test_invalid_requests_are_actionable_and_profile_boundary_is_preserved():
    client = TestClient(
        create_app(runtime=PlanningRuntime(), recommendation_only=False)
    )
    url = "/api/v1/retrosynthesis/planner"
    assert (
        client.post(url, json={"action": "start", "target_smiles": "bad"}).status_code
        == 422
    )
    assert client.post(url, json={"action": "select"}).status_code == 422
    assert (
        client.post(
            url, json={"action": "restore", "session": {"schema_version": "unknown"}}
        ).status_code
        == 422
    )
    session = client.post(url, json={"action": "start", "target_smiles": "CCO"}).json()[
        "data"
    ]["session"]
    assert (
        client.post(
            url, json={"action": "stop", "session": session, "node_id": "missing"}
        ).status_code
        == 422
    )
    assert (
        client.post(
            url, json={"action": "select", "session": session, "strategy_index": -1}
        ).status_code
        == 422
    )
    restricted = TestClient(
        create_app(runtime=PlanningRuntime(), recommendation_only=True)
    )
    assert (
        restricted.post(
            url, json={"action": "start", "target_smiles": "CCO"}
        ).status_code
        == 404
    )


def test_stock_check_uses_exact_supplier_evidence_without_automatically_stopping(
    tmp_path,
):
    from cas_tools import StockSourceDefinition, build_stock_portfolio
    from app.web_api.runtime import LocalRecommendationRuntime

    source = tmp_path / "stock.smi"
    source.write_text("CCO\tstock-1\n", encoding="utf-8")
    stock_path = tmp_path / "stock.sqlite"
    build_stock_portfolio(
        (
            StockSourceDefinition(
                path=str(source),
                supplier="Test supplier",
                collection="In stock",
                snapshot_date="2026-10-07",
                availability_status="in_stock",
                evidence_level="supplier_in_stock",
                terminal_eligible=True,
                source_url="https://example.test/stock",
                terms_url="https://example.test/terms",
            ),
        ),
        stock_path,
    )
    runtime = LocalRecommendationRuntime(stock_portfolio_path=stock_path)
    assert runtime.planning_stock("OCC")["status"] == "verified_stock_match"
    assert runtime.planning_stock("CCN")["status"] == "no_verified_match"
    unavailable = LocalRecommendationRuntime(
        stock_portfolio_path=tmp_path / "missing.sqlite"
    )
    assert unavailable.planning_stock("CCO")["status"] == "unavailable"
    client = TestClient(create_app(runtime=runtime, recommendation_only=False))
    url = "/api/v1/retrosynthesis/planner"
    session = client.post(url, json={"action": "start", "target_smiles": "OCC"}).json()[
        "data"
    ]["session"]
    checked = client.post(url, json={"action": "stock", "session": session})
    assert checked.status_code == 200
    result = checked.json()["data"]
    assert (
        result["session"]["stock"]["CCO"]["source_records"][0]["supplier"]
        == "Test supplier"
    )
    assert result["summary"]["unresolved_count"] == 1
    assert not result["session"]["root"]["stopped"]
