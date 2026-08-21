"""Shared schema-constrained assistance transport tests."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest
from pydantic import ValidationError

from chem_coworker.assistance.transport import (
    AssistanceActionPayload,
    AssistanceProviderSettings,
    OpenAICompatibleAssistanceTransport,
)


def _payload() -> dict:
    return {
        "action_name": "finish",
        "arguments": {
            "terminal_status": "abstained_insufficient_evidence",
            "stopping_reason": "No authoritative evidence is available.",
        },
        "cited_evidence_ids": [],
        "decision_summary": "Abstain instead of inventing evidence.",
        "claims": [],
    }


def test_action_schema_rejects_unknown_actions_and_arguments() -> None:
    invalid_action = {**_payload(), "action_name": "execute_python"}
    invalid_argument = _payload()
    invalid_argument["arguments"]["precursor_smiles"] = "invented"

    with pytest.raises(ValidationError):
        AssistanceActionPayload.model_validate(invalid_action)
    with pytest.raises(ValidationError):
        AssistanceActionPayload.model_validate(invalid_argument)


def test_openai_transport_uses_strict_no_store_schema_and_one_retry(
    monkeypatch,
) -> None:
    monkeypatch.setenv("OPENAI_API_KEY", "fixture-key")
    calls = []
    responses = [
        SimpleNamespace(
            status="incomplete",
            output_text="",
            usage=SimpleNamespace(input_tokens=7, output_tokens=3),
        ),
        SimpleNamespace(
            id="response-2",
            status="completed",
            output_text=json.dumps(_payload()),
            usage=SimpleNamespace(input_tokens=11, output_tokens=5),
        ),
    ]

    class Responses:
        def create(self, **kwargs):
            calls.append(kwargs)
            return responses.pop(0)

    client = SimpleNamespace(responses=Responses())
    transport = OpenAICompatibleAssistanceTransport(
        client_factory=lambda **kwargs: client
    )

    result = transport.complete(
        {"allowed_actions": ["finish"], "evidence": []},
        AssistanceProviderSettings(),
        max_retries=1,
    )

    assert result.attempts == 2
    assert result.input_tokens == 18
    assert result.output_tokens == 8
    assert calls[0]["store"] is False
    assert calls[0]["text"]["format"]["strict"] is True
    assert calls[0]["model"] == "gpt-5.6-terra"


def test_provider_settings_reject_unknown_configuration() -> None:
    with pytest.raises(ValueError, match="unknown provider settings"):
        AssistanceProviderSettings.from_mapping({"temperature": 1.0})
