"""
Phase 6: CriticAgent tests.

All tests use mock LLMs — no real API calls.
"""
import json
import pytest
from unittest.mock import MagicMock, patch


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

def _make_llm_response(text: str) -> MagicMock:
    """Minimal mock of a LangChain AIMessage."""
    msg = MagicMock()
    msg.content = text
    msg.tool_calls = []
    return msg


def _make_critic_response(findings: list, overall: str = "looks reasonable") -> str:
    """Build a valid critic JSON response string."""
    return json.dumps({"findings": findings, "overall": overall})


def _make_agent_with_mock_llm(mock_response_text: str):
    """Create a ChemCoworker instance with a mocked LLM (no API key needed)."""
    from chem_coworker.agent import ChemCoworker
    from chem_coworker.tools import REGISTRY
    agent = object.__new__(ChemCoworker)
    agent.registry = REGISTRY
    agent.verbose = False
    from chem_coworker.event_bus import EventBus
    agent.event_bus = EventBus()
    mock_llm = MagicMock()
    mock_llm.invoke.return_value = _make_llm_response(mock_response_text)
    agent.llm = mock_llm
    return agent


# ---------------------------------------------------------------------------
# Severity tests
# ---------------------------------------------------------------------------

class TestSeverity:
    def test_values(self):
        from chem_coworker.critic import Severity
        assert Severity.INFO.value == "info"
        assert Severity.WARNING.value == "warning"
        assert Severity.CRITICAL.value == "critical"

    def test_ordering_info_lt_warning(self):
        from chem_coworker.critic import Severity
        assert not (Severity.INFO >= Severity.WARNING)

    def test_ordering_warning_lt_critical(self):
        from chem_coworker.critic import Severity
        assert not (Severity.WARNING >= Severity.CRITICAL)

    def test_ordering_equal(self):
        from chem_coworker.critic import Severity
        assert Severity.WARNING >= Severity.WARNING

    def test_ordering_critical_gte_warning(self):
        from chem_coworker.critic import Severity
        assert Severity.CRITICAL >= Severity.WARNING

    def test_ordering_critical_gte_info(self):
        from chem_coworker.critic import Severity
        assert Severity.CRITICAL >= Severity.INFO


# ---------------------------------------------------------------------------
# Finding tests
# ---------------------------------------------------------------------------

class TestFinding:
    def test_str_warning(self):
        from chem_coworker.critic import Finding, Severity
        f = Finding(severity=Severity.WARNING, message="Bad leaving group")
        s = str(f)
        assert "WARNING" in s
        assert "Bad leaving group" in s

    def test_str_critical_with_suggestion(self):
        from chem_coworker.critic import Finding, Severity
        f = Finding(severity=Severity.CRITICAL, message="No stereodirecting group", suggestion="Add chiral auxiliary")
        s = str(f)
        assert "CRITICAL" in s
        assert "No stereodirecting group" in s
        assert "Add chiral auxiliary" in s

    def test_str_info_no_suggestion(self):
        from chem_coworker.critic import Finding, Severity
        f = Finding(severity=Severity.INFO, message="5-step route is long")
        s = str(f)
        assert "INFO" in s
        assert "→" not in s  # no suggestion arrow

    def test_default_suggestion_is_empty(self):
        from chem_coworker.critic import Finding, Severity
        f = Finding(severity=Severity.INFO, message="note")
        assert f.suggestion == ""


# ---------------------------------------------------------------------------
# _parse_critic_response tests
# ---------------------------------------------------------------------------

class TestParseCriticResponse:
    def test_empty_findings(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = json.dumps({"findings": [], "overall": "route is sound"})
        findings, overall = _parse_critic_response(raw, max_findings=5, min_severity=Severity.INFO)
        assert findings == []
        assert overall == "route is sound"

    def test_single_warning_finding(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = json.dumps({
            "findings": [{"severity": "warning", "message": "selectivity risk", "suggestion": "use protecting group"}],
            "overall": "minor issues"
        })
        findings, overall = _parse_critic_response(raw, max_findings=5, min_severity=Severity.INFO)
        assert len(findings) == 1
        assert findings[0].message == "selectivity risk"
        assert findings[0].suggestion == "use protecting group"

    def test_min_severity_filters_info(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = json.dumps({
            "findings": [
                {"severity": "info", "message": "long route"},
                {"severity": "critical", "message": "no selectivity"},
            ],
            "overall": "has issues"
        })
        findings, overall = _parse_critic_response(raw, max_findings=5, min_severity=Severity.WARNING)
        assert len(findings) == 1
        assert findings[0].severity == Severity.CRITICAL

    def test_max_findings_cap(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = json.dumps({
            "findings": [{"severity": "warning", "message": f"issue {i}"} for i in range(10)],
            "overall": "many issues"
        })
        findings, _ = _parse_critic_response(raw, max_findings=3, min_severity=Severity.INFO)
        assert len(findings) <= 3

    def test_invalid_json_returns_fallback_finding(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        findings, overall = _parse_critic_response("not valid json", max_findings=5, min_severity=Severity.INFO)
        assert len(findings) == 1
        assert findings[0].severity == Severity.INFO
        assert overall == "(parse failed)"

    def test_markdown_fence_stripped(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = "```json\n" + json.dumps({"findings": [], "overall": "ok"}) + "\n```"
        findings, overall = _parse_critic_response(raw, max_findings=5, min_severity=Severity.INFO)
        assert overall == "ok"

    def test_unknown_severity_treated_as_info(self):
        from chem_coworker.critic import _parse_critic_response, Severity
        raw = json.dumps({
            "findings": [{"severity": "unknown_level", "message": "msg"}],
            "overall": "test"
        })
        findings, _ = _parse_critic_response(raw, max_findings=5, min_severity=Severity.INFO)
        assert len(findings) == 1
        assert findings[0].severity == Severity.INFO


# ---------------------------------------------------------------------------
# CriticAgent.review() tests (mocked LLM)
# ---------------------------------------------------------------------------

class TestCriticAgentReview:
    def test_review_returns_findings_and_verdict(self):
        from chem_coworker.critic import CriticAgent, Severity
        response_text = _make_critic_response(
            findings=[{"severity": "warning", "message": "protect amine"}],
            overall="route needs adjustment"
        )
        llm = MagicMock()
        llm.invoke.return_value = _make_llm_response(response_text)

        critic = CriticAgent(llm)
        findings, verdict = critic.review(
            query="Synthesize aspirin",
            hypothesis="acetylation of salicylic acid",
            tool_results={},
            answer="Step 1: ...",
        )
        assert len(findings) == 1
        assert "protect amine" in findings[0].message
        assert "adjustment" in verdict

    def test_review_empty_findings_returns_empty_list(self):
        from chem_coworker.critic import CriticAgent, Severity
        response_text = _make_critic_response(findings=[], overall="route is valid")
        llm = MagicMock()
        llm.invoke.return_value = _make_llm_response(response_text)

        critic = CriticAgent(llm)
        findings, verdict = critic.review("query", "hyp", {}, "answer")
        assert findings == []
        assert "valid" in verdict

    def test_review_llm_failure_returns_empty(self):
        from chem_coworker.critic import CriticAgent
        llm = MagicMock()
        llm.invoke.side_effect = RuntimeError("API error")

        critic = CriticAgent(llm)
        findings, verdict = critic.review("q", "h", {}, "a")
        assert findings == []
        assert "unavailable" in verdict.lower()

    def test_review_min_severity_respected(self):
        from chem_coworker.critic import CriticAgent, Severity
        response_text = _make_critic_response(
            findings=[
                {"severity": "info", "message": "step count"},
                {"severity": "critical", "message": "selectivity issue"},
            ],
            overall="serious problem"
        )
        llm = MagicMock()
        llm.invoke.return_value = _make_llm_response(response_text)

        critic = CriticAgent(llm)
        findings, _ = critic.review(
            "q", "h", {}, "a",
            min_severity=Severity.CRITICAL,
        )
        assert len(findings) == 1
        assert findings[0].severity == Severity.CRITICAL

    def test_llm_called_once(self):
        from chem_coworker.critic import CriticAgent
        response_text = _make_critic_response(findings=[], overall="ok")
        llm = MagicMock()
        llm.invoke.return_value = _make_llm_response(response_text)

        critic = CriticAgent(llm)
        critic.review("q", "h", {}, "a")
        llm.invoke.assert_called_once()


# ---------------------------------------------------------------------------
# _run_critic_loop() integration (via ChemCoworker method)
# ---------------------------------------------------------------------------

class TestRunCriticLoop:
    def test_critic_loop_returns_findings_and_call_count(self):
        from chem_coworker.workflow import CriticStep
        response_text = _make_critic_response(
            findings=[{"severity": "warning", "message": "chemo issue"}],
            overall="needs care"
        )
        agent = _make_agent_with_mock_llm(response_text)
        critic_step = CriticStep(enabled=True, max_findings=3, min_severity="warning")

        findings, verdict, calls = agent._run_critic_loop(
            query="target molecule synthesis",
            hypothesis="Wittig route",
            tool_results={},
            answer="your synthesis...",
            critic_step=critic_step,
        )
        assert len(findings) == 1
        assert calls == 1
        assert "needs care" in verdict

    def test_critic_loop_llm_failure_returns_empty_findings(self):
        """When the LLM fails inside CriticAgent.review(), findings are empty.
        The call is still counted (1) because the LLM was invoked."""
        from chem_coworker.workflow import CriticStep
        agent = _make_agent_with_mock_llm("")
        agent.llm.invoke.side_effect = RuntimeError("network error")
        critic_step = CriticStep(enabled=True)

        findings, verdict, calls = agent._run_critic_loop(
            query="q", hypothesis="h", tool_results={}, answer="a",
            critic_step=critic_step,
        )
        assert findings == []
        # CriticAgent.review() gracefully handles the LLM error and returns a verdict
        assert "unavailable" in verdict.lower()
        # The LLM was invoked (and failed), so count is 1
        assert calls == 1

    def test_critic_loop_emits_phase_events(self):
        from chem_coworker.workflow import CriticStep
        from chem_coworker.event_bus import EventBus, ChemEvent
        phases = []
        bus = EventBus()
        bus.subscribe(ChemEvent.PHASE_START, lambda phase, **_: phases.append(f"start:{phase}"))
        bus.subscribe(ChemEvent.PHASE_END,   lambda phase, **_: phases.append(f"end:{phase}"))

        response_text = _make_critic_response(findings=[], overall="ok")
        agent = _make_agent_with_mock_llm(response_text)
        agent.event_bus = bus
        critic_step = CriticStep(enabled=True)

        agent._run_critic_loop("q", "h", {}, "a", critic_step)
        assert "start:critic" in phases
        assert "end:critic" in phases


# ---------------------------------------------------------------------------
# CriticStep dataclass tests
# ---------------------------------------------------------------------------

class TestCriticStep:
    def test_defaults(self):
        from chem_coworker.workflow import CriticStep
        cs = CriticStep()
        assert cs.enabled is True
        assert cs.max_findings == 5
        assert cs.min_severity == "warning"

    def test_custom_values(self):
        from chem_coworker.workflow import CriticStep
        cs = CriticStep(enabled=False, max_findings=10, min_severity="critical")
        assert cs.enabled is False
        assert cs.max_findings == 10
        assert cs.min_severity == "critical"


# ---------------------------------------------------------------------------
# Public exports test
# ---------------------------------------------------------------------------

class TestPublicExports:
    def test_critic_types_exported(self):
        import chem_coworker as cc
        assert hasattr(cc, "CriticAgent")
        assert hasattr(cc, "Finding")
        assert hasattr(cc, "Severity")
        assert hasattr(cc, "CRITIC_SYSTEM_PROMPT")
        assert hasattr(cc, "CriticStep")

    def test_critic_system_prompt_is_str(self):
        from chem_coworker import CRITIC_SYSTEM_PROMPT
        assert isinstance(CRITIC_SYSTEM_PROMPT, str)
        assert len(CRITIC_SYSTEM_PROMPT) > 100
