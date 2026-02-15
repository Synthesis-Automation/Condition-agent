from langchain_core.messages import AIMessage, ToolMessage

from chem_assistant import chemtools_agent as agent_module


def test_extract_tool_payload_parses_python_literal_mapping() -> None:
    messages = [
        AIMessage(
            content="tool call",
            tool_calls=[
                {
                    "id": "call_1",
                    "name": "hte_recommend_conditions",
                    "args": {},
                }
            ],
        ),
        ToolMessage(
            content="{'success': True, 'recommendations': []}",
            tool_call_id="call_1",
        ),
    ]

    payload = agent_module._extract_tool_payload(
        messages,
        "hte_recommend_conditions",
    )

    assert isinstance(payload, dict)
    assert payload.get("success") is True
    assert payload.get("recommendations") == []


def test_build_condition_response_handles_non_mapping_blocks() -> None:
    response = agent_module._build_condition_response(
        query="recommend conditions",
        tool_payload={
            "success": False,
            "error": "tool failed",
            "input": "invalid",
            "diagnostics": "invalid",
            "recommendations": "invalid",
        },
        kb_payload="invalid",  # type: ignore[arg-type]
    )

    assert response.query == "recommend conditions"
    assert isinstance(response.warnings, list)
    assert "hte_recommend_conditions failed." in response.warnings
    assert "tool failed" in response.warnings


class _DummyTool:
    def __init__(self, payload):
        self._payload = payload

    def invoke(self, _args):
        return self._payload


class _DummyBackend:
    def __init__(self, payload):
        self._payload = payload

    def invoke(self, _args, config=None):
        return self._payload


def test_preflight_handles_string_payloads(monkeypatch) -> None:
    from chem_assistant import chemtools_wrapper as wrapper

    monkeypatch.setattr(
        wrapper,
        "unified_featurize_reaction",
        _DummyTool("not-a-dict"),
    )
    monkeypatch.setattr(
        wrapper,
        "analysis_analyze_reaction",
        _DummyTool("not-a-dict"),
    )

    agent = object.__new__(agent_module.ChemToolsAgent)
    agent.auto_check = True

    report = agent._maybe_run_preflight(
        "recommend conditions for reaction: A.B>>C"
    )

    assert isinstance(report, str)
    assert "unexpected payload: str" in report


def test_run_handles_non_mapping_agent_result() -> None:
    agent = object.__new__(agent_module.ChemToolsAgent)
    agent.auto_check = False
    agent.verbose = False
    agent.agent = _DummyBackend("raw-agent-output")

    response = agent.run("hello")

    assert response == "raw-agent-output"
