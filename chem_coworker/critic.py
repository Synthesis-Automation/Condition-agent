"""
CriticAgent — Adversarial reviewer for retrosynthesis proposals.

Runs after the main native tool loop for retrosynthesis tasks, taking the
proposed route and tool results and applying an adversarial lens:
  • Chemoselectivity issues
  • Protecting-group strategy gaps
  • Stereochemical concerns
  • Step-count and atom-economy notes
  • Reagent availability red flags
  • Reaction feasibility under stated conditions

Returns a list of structured Finding objects plus an optional revised note
that can be appended to the main answer.

Designed to be cheap (single focused LLM call) and optional (disabled by
default for forward-synthesis workflows; enabled for retrosynthesis).
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

class Severity(str, Enum):
    """Ordered severity levels for critic findings."""
    INFO     = "info"
    WARNING  = "warning"
    CRITICAL = "critical"

    def __ge__(self, other: "Severity") -> bool:  # type: ignore[override]
        _order = {Severity.INFO: 0, Severity.WARNING: 1, Severity.CRITICAL: 2}
        return _order[self] >= _order[other]


@dataclass
class Finding:
    """A single critic observation about a synthesis proposal."""
    severity: Severity
    message: str
    suggestion: str = ""

    def __str__(self) -> str:
        icon = {"info": "ℹ", "warning": "⚠", "critical": "✘"}.get(self.severity.value, "•")
        parts = [f"{icon} [{self.severity.value.upper()}] {self.message}"]
        if self.suggestion:
            parts.append(f"  → {self.suggestion}")
        return "\n".join(parts)


# ---------------------------------------------------------------------------
# Critic system prompt
# ---------------------------------------------------------------------------

CRITIC_SYSTEM_PROMPT = """\
You are an adversarial chemistry reviewer specialising in retrosynthesis.
Your role is to critically evaluate a proposed synthesis route and identify
genuine problems — not to repeat the route or re-describe it.

Focus on:
1. CHEMOSELECTIVITY — Will the proposed reagents react with unintended functional groups?
2. PROTECTING GROUPS — Are competing nucleophiles or electrophiles left unprotected?
3. STEREOCHEMISTRY — Does the proposed route achieve required stereocontrol?
4. STEP COUNT / ATOM ECONOMY — Obvious inefficiencies or better known alternatives.
5. REAGENT AVAILABILITY — Highly exotic or unavailable reagents.
6. FEASIBILITY — Temperature, pressure, or safety concerns not addressed.

RESPONSE FORMAT: Return a JSON object (no markdown fences):
{
  "findings": [
    {
      "severity": "info"|"warning"|"critical",
      "message": "concise description of the issue",
      "suggestion": "what to do instead (optional)"
    }
  ],
  "overall": "brief one-sentence verdict"
}

Return at most 5 findings. Only list real problems — if the route looks sound, say so in "overall" with an empty findings list.
"""


# ---------------------------------------------------------------------------
# CriticAgent
# ---------------------------------------------------------------------------

class CriticAgent:
    """
    Runs a single focused adversarial LLM call to review a retrosynthesis
    proposal and return structured findings.
    """

    def __init__(self, llm: Any):
        self.llm = llm

    def review(
        self,
        query: str,
        hypothesis: str,
        tool_results: Dict[str, Any],
        answer: str,
        max_findings: int = 5,
        min_severity: Severity = Severity.WARNING,
    ) -> tuple:  # -> (List[Finding], str)  — py3.8 compat
        """
        Review a synthesis proposal.

        Args:
            query:        Original user query (target molecule / task description).
            hypothesis:   The agent's stated hypothesis / synthesis strategy.
            tool_results: Dict of tool name → result dict from the main loop.
            answer:       The agent's final answer to be critiqued.
            max_findings: Cap on number of findings returned.
            min_severity: Minimum severity level to include in returned findings.

        Returns:
            (findings, overall)  where findings is a List[Finding] and
            overall is a one-sentence verdict string.
        """
        from langchain_core.messages import SystemMessage, HumanMessage

        # Build a concise summary of tool evidence for the critic
        tool_summary_parts: List[str] = []
        for tool_name, result in tool_results.items():
            if not isinstance(result, dict):
                continue
            display = {k: v for k, v in result.items()
                       if k not in ("success", "_warnings") and v is not None}
            try:
                summary_str = json.dumps(display, default=str)[:600]
            except Exception:
                summary_str = str(display)[:600]
            tool_summary_parts.append(f"[{tool_name}]\n{summary_str}")

        tool_summary = "\n\n".join(tool_summary_parts) if tool_summary_parts else "(no tool results)"

        critic_input = (
            f"ORIGINAL QUERY:\n{query}\n\n"
            f"PROPOSED STRATEGY:\n{hypothesis or '(not stated)'}\n\n"
            f"TOOL EVIDENCE:\n{tool_summary}\n\n"
            f"AGENT ANSWER (to critique):\n{answer[:2000] if answer else '(empty)'}"
        )

        messages = [
            SystemMessage(content=CRITIC_SYSTEM_PROMPT),
            HumanMessage(content=critic_input),
        ]

        try:
            response = self.llm.invoke(messages)
            raw = _extract_text(response)
        except Exception as exc:
            logger.warning(f"[CriticAgent] LLM call failed: {exc}")
            return [], f"Critic unavailable: {exc}"

        findings, overall = _parse_critic_response(raw, max_findings, min_severity)
        return findings, overall


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _extract_text(llm_response: Any) -> str:
    content = getattr(llm_response, "content", llm_response)
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        return "\n".join(
            str(item.get("text") or item.get("content", "")) if isinstance(item, dict) else str(item)
            for item in content
        )
    return str(content)


def _parse_critic_response(
    raw: str,
    max_findings: int,
    min_severity: Severity,
) -> tuple:  # (List[Finding], str)
    """Parse the JSON critic response, gracefully degrading on bad output."""
    # Strip markdown fences if present
    text = raw.strip()
    if text.startswith("```"):
        lines = text.splitlines()
        text = "\n".join(lines[1:-1] if lines[-1].strip() == "```" else lines[1:])

    try:
        data = json.loads(text)
    except Exception:
        # Fallback: return a single INFO finding with the raw text
        return [Finding(severity=Severity.INFO, message=raw[:300])], "(parse failed)"

    findings: List[Finding] = []
    for item in data.get("findings", [])[:max_findings]:
        try:
            sev = Severity(item.get("severity", "info"))
        except ValueError:
            sev = Severity.INFO
        if sev >= min_severity:
            findings.append(Finding(
                severity=sev,
                message=str(item.get("message", "")),
                suggestion=str(item.get("suggestion", "")),
            ))

    overall = str(data.get("overall", ""))
    return findings, overall
