"""
ExecutionPlan and PlanParser for ChemCoworker.

The LLM (in REASON_PROMPT step) produces a JSON plan at the end of its reasoning text.
PlanParser extracts and validates that plan, producing an ExecutionPlan that
the ToolExecutor can run.

Plan JSON format (expected from LLM):
  {
    "hypothesis": "Suzuki-Miyaura coupling (C-C bond, aryl bromide + boronic acid)",
    "confidence": 0.9,
    "groups": [
      [{"name": "normalize_reaction", "args": {"smiles": "..."}},
       {"name": "detect_reaction_type", "args": {"reaction_smiles": "..."}}],
      [{"name": "analyze_bond_changes", "args": {"reaction_smiles": "..."}}],
      [{"name": "recommend_conditions", "args": {"reaction_smiles": "...", "top_k": 5}}]
    ],
    "rationale": "HIGH confidence Suzuki hypothesis; minimal tools needed."
  }

Groups are sequential; tools within each group run in parallel.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


class PlanRejected(Exception):
    """Raised by a plan_callback to cancel tool execution (A2 — plan approval)."""


@dataclass
class ToolCall:
    """A single tool call specification (name + arguments)."""
    name: str
    args: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        args_str = ", ".join(f"{k}={v!r}" for k, v in self.args.items())
        return f"{self.name}({args_str})"


@dataclass
class ExecutionPlan:
    """
    Structured execution plan produced by the LLM's reasoning phase.

    Attributes:
        hypothesis: Named reaction / chemistry concept identified by the LLM.
        confidence: 0.0–1.0 reflecting LLM certainty about the hypothesis.
        groups: Sequential groups of ToolCalls. Within each group, calls are
                independent and will run concurrently.
        rationale: Brief explanation of why these tools were chosen.
        raw_plan_text: The full LLM text from Phase 1 (reasoning + plan JSON).
    """
    hypothesis: str = ""
    confidence: float = 0.5
    groups: List[List[ToolCall]] = field(default_factory=list)
    rationale: str = ""
    raw_plan_text: str = ""

    @property
    def all_tool_calls(self) -> List[ToolCall]:
        return [tc for grp in self.groups for tc in grp]

    @property
    def tool_names(self) -> List[str]:
        return [tc.name for tc in self.all_tool_calls]

    @property
    def is_empty(self) -> bool:
        """True if no tools are planned (pure LLM reasoning task)."""
        return len(self.all_tool_calls) == 0

    def __repr__(self) -> str:
        groups_repr = " → ".join(
            f"[{', '.join(tc.name for tc in grp)}]" for grp in self.groups
        )
        return (
            f"ExecutionPlan(hypothesis={self.hypothesis!r}, "
            f"confidence={self.confidence:.2f}, groups={groups_repr})"
        )


class PlanParser:
    """
    Extracts and validates an ExecutionPlan from LLM Phase 1 output text.

    Strategies (in order):
    1. JSON inside ```json ... ``` fence
    2. Last {...} block in the text
    3. Entire text as JSON
    4. Fallback: empty plan (LLM will answer from knowledge alone)
    """

    def parse(
        self,
        llm_text: str,
        known_tools: Optional[List[str]] = None,
        smiles_context: Optional[str] = None,
    ) -> ExecutionPlan:
        """
        Parse LLM Phase 1 output into an ExecutionPlan.

        Args:
            llm_text: Full text from the LLM's reasoning phase.
            known_tools: List of registered tool names (for validation).
            smiles_context: Primary SMILES from query (for filling in missing args).

        Returns:
            ExecutionPlan (may be empty if LLM chose to answer from knowledge).
        """
        data = self._extract_json(llm_text)

        if data is None:
            # LLM didn't produce a parseable plan — treat as knowledge-only response
            return ExecutionPlan(
                hypothesis=self._extract_hypothesis_text(llm_text),
                confidence=0.5,
                groups=[],
                rationale="No structured plan found; answering from LLM knowledge.",
                raw_plan_text=llm_text,
            )

        hypothesis = str(data.get("hypothesis", ""))
        confidence = float(data.get("confidence", 0.5))
        confidence = max(0.0, min(1.0, confidence))
        rationale = str(data.get("rationale", ""))

        raw_groups = data.get("groups", [])
        groups: List[List[ToolCall]] = []
        for raw_group in raw_groups:
            if not isinstance(raw_group, list):
                continue
            group: List[ToolCall] = []
            for item in raw_group:
                if not isinstance(item, dict):
                    continue
                name = str(item.get("name", "")).strip()
                if not name:
                    continue
                # Validate against known tools if provided
                if known_tools and name not in known_tools:
                    continue
                args = dict(item.get("args", {}))
                # Fill in primary SMILES if args reference it but it's missing
                if smiles_context:
                    args = self._fill_smiles_args(args, smiles_context)
                group.append(ToolCall(name=name, args=args))
            if group:
                groups.append(group)

        return ExecutionPlan(
            hypothesis=hypothesis,
            confidence=confidence,
            groups=groups,
            rationale=rationale,
            raw_plan_text=llm_text,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _extract_json(self, text: str) -> Optional[Dict[str, Any]]:
        """Try multiple strategies to extract a JSON plan from text."""
        if not text:
            return None

        # Strip control characters (MiniMax issue)
        text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", "", text)

        # Strategy 1: JSON in code fence
        fence = re.search(r"```(?:json)?\s*\n?([\s\S]*?)\n?```", text)
        if fence:
            try:
                data = json.loads(fence.group(1).strip())
                if isinstance(data, dict) and "groups" in data:
                    return data
            except json.JSONDecodeError:
                pass

        # Strategy 2: All {...} spans — greedy match finds outermost blocks
        # Try each starting position; reversed so we try the last (often outermost) first
        blocks = list(re.finditer(r"\{[\s\S]*\}", text))
        for block in reversed(blocks):
            try:
                data = json.loads(block.group(0))
                if isinstance(data, dict) and "groups" in data:
                    return data
            except json.JSONDecodeError:
                pass

        # Strategy 3: Full text as JSON
        try:
            data = json.loads(text.strip())
            if isinstance(data, dict) and "groups" in data:
                return data
        except json.JSONDecodeError:
            pass

        return None

    def _extract_hypothesis_text(self, text: str) -> str:
        """Try to extract a hypothesis statement from free-form LLM text."""
        # Look for "HYPOTHESIS: ..." pattern
        m = re.search(r"hypothesis[:\s]+([^\n.]{10,200})", text, re.IGNORECASE)
        if m:
            return m.group(1).strip()
        # First sentence as fallback
        sentences = text.strip().split(".")
        return sentences[0].strip() if sentences else ""

    def _fill_smiles_args(self, args: Dict[str, Any], smiles: str) -> Dict[str, Any]:
        """Fill in SMILES arguments that the LLM left as placeholders."""
        result = dict(args)
        placeholder_pattern = re.compile(
            r"^(\.\.\.|\<.*?\>|SMILES_HERE|reaction_smiles_here|YOUR_SMILES)$",
            re.IGNORECASE,
        )
        for key, val in result.items():
            if isinstance(val, str) and (placeholder_pattern.match(val) or not val):
                result[key] = smiles
        return result
