"""Experimental API wiring for the common assistance service."""

from __future__ import annotations

from fastapi.testclient import TestClient

from app.web_api.main import create_app


class _Service:
    def __init__(self) -> None:
        self.started = None
        self.confirmed = None

    def start(self, **kwargs):
        self.started = kwargs
        return {"state": {"status": "completed"}, "rollout_state": "advisory"}

    def confirm_condition_constraint(self, **kwargs):
        self.confirmed = kwargs
        return {"state": {"status": "completed"}, "rollout_state": "advisory"}


def test_assistance_endpoint_is_disabled_unless_explicitly_configured() -> None:
    client = TestClient(create_app(runtime=object(), frontend_dist="missing"))

    response = client.post(
        "/api/v1/experimental/assistance",
        json={
            "objective": "Explain conditions",
            "mode": "conditions",
            "structure_input": "CCBr.O>>CCO",
        },
    )

    assert response.status_code == 503
    assert response.json()["detail"]["error"] == "ASSISTANCE_NOT_CONFIGURED"


def test_assistance_endpoint_uses_injected_common_service() -> None:
    service = _Service()
    client = TestClient(
        create_app(
            runtime=object(),
            assistance_service=service,
            frontend_dist="missing",
        )
    )

    response = client.post(
        "/api/v1/experimental/assistance",
        json={
            "objective": "Explain conditions",
            "mode": "conditions",
            "structure_input": "CCBr.O>>CCO",
            "provider": {"provider": "openai", "model": "gpt-5.6-terra"},
        },
    )

    assert response.status_code == 200
    assert response.json()["data"]["rollout_state"] == "advisory"
    assert service.started["mode"] == "conditions"


def test_confirmation_endpoint_passes_serialized_state_and_raw_answer() -> None:
    service = _Service()
    client = TestClient(
        create_app(
            runtime=object(),
            assistance_service=service,
            frontend_dist="missing",
        )
    )

    response = client.post(
        "/api/v1/experimental/assistance/confirm-condition",
        json={"state": {"status": "needs_user_input"}, "raw_value": "Pd(PPh3)4"},
    )

    assert response.status_code == 200
    assert service.confirmed["raw_value"] == "Pd(PPh3)4"
