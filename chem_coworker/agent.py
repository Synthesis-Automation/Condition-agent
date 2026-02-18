"""
ChemCoworker — General-purpose chemistry AI agent.

"Claude Code for chemistry" — combines LLM chemistry expertise with a
lightweight tool ecosystem. Follows a 6-step workflow:

  1. INTAKE    : classify task, extract SMILES
  2. REASON    : LLM reasons from chemistry knowledge → produces execution plan
  3. PARSE     : extract ExecutionPlan from LLM output
  4. EXECUTE   : run tool groups in parallel batches
  5. VALIDATE  : check if results match hypothesis (optional revision)
  6. SYNTHESIZE: LLM writes final expert answer

Total LLM calls: 2 (reason + synthesize). 3 only if hypothesis needs revision.
"""
from __future__ import annotations

import json
import logging
import os
import time
from typing import Any, Dict, List, Optional, Tuple

from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)

# Local imports (no circular dependency — these modules don't import agent.py)
from .response import ChemResponse  # noqa: E402
from .plan import ExecutionPlan     # noqa: E402


# ---------------------------------------------------------------------------
# LLM client factory (reuses pattern from reasoning_agent.py)
# ---------------------------------------------------------------------------

def _get_llm_client(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0,
) -> Any:
    """Get a LangChain ChatOpenAI client for the given provider/model."""
    from langchain_openai import ChatOpenAI

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "openai")
    if model is None:
        model = os.getenv("LLM_MODEL", "o4-mini")

    if provider == "aliyun":
        api_key = os.getenv("ALIYUN_API_KEY")
        base_url = os.getenv(
            "ALIYUN_BASE_URL",
            "https://dashscope.aliyuncs.com/compatible-mode/v1",
        )
    elif provider == "openai":
        api_key = os.getenv("OPENAI_API_KEY")
        base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")
    else:
        raise ValueError(f"Unsupported provider '{provider}'. Use 'openai' or 'aliyun'.")

    if not api_key:
        raise RuntimeError(f"{provider.upper()}_API_KEY environment variable not set")

    # OpenAI o-series and gpt-5.x models do not accept temperature
    _no_temp_prefixes = ("o1", "o3", "o4", "gpt-5")
    is_reasoning_model = (provider == "openai") and any(
        model.startswith(p) for p in _no_temp_prefixes
    )

    kwargs: Dict[str, Any] = {"model": model, "api_key": api_key, "base_url": base_url}
    if not is_reasoning_model:
        kwargs["temperature"] = temperature

    return ChatOpenAI(**kwargs)


# ---------------------------------------------------------------------------
# ChemCoworker
# ---------------------------------------------------------------------------

class ChemCoworker:
    """
    General-purpose chemistry AI agent.

    Handles any chemistry question: reaction analysis, condition prediction,
    mechanism explanation, reagent lookup, troubleshooting, and more.

    Usage:
        coworker = ChemCoworker(provider="openai", model="o4-mini")
        response = coworker.run("Recommend conditions for Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        print(response.answer)
        print(response.conditions)

    Multi-turn:
        history = []
        resp, history = coworker.chat("What is this reaction?", history)
        resp, history = coworker.chat("Now recommend conditions", history)
    """

    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0,
        verbose: bool = False,
    ):
        self.provider = provider or os.getenv("LLM_PROVIDER", "openai")
        self.model_name = model or os.getenv("LLM_MODEL", "o4-mini")
        self.verbose = verbose

        self.llm = _get_llm_client(provider, model, temperature)

        from .tools import REGISTRY
        from .classifier import TaskClassifier
        from .plan import PlanParser
        from .executor import ToolExecutor

        self.registry = REGISTRY
        self.classifier = TaskClassifier()
        self.parser = PlanParser()
        self.executor = ToolExecutor(verbose=verbose)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, query: str) -> "ChemResponse":
        """
        Single-turn run. Returns a ChemResponse with plan, tool results,
        and final synthesized answer.
        """
        from .response import ChemResponse
        start = time.monotonic()
        llm_calls = 0

        from .prompts import REASON_PROMPT, SYNTHESIZE_PROMPT
        from langchain_core.messages import HumanMessage

        # ── Step 1: Intake ─────────────────────────────────────────────
        clf = self.classifier.classify(query)
        task_type = clf.task_type_str
        smiles_list = clf.all_smiles
        primary_smiles = clf.primary_smiles or ""

        if self.verbose:
            logger.info(f"[ChemCoworker] task={task_type} smiles={smiles_list}")

        reason_text = REASON_PROMPT.format(
            task_type=task_type.upper(),
            query=query,
            smiles_list=", ".join(smiles_list) if smiles_list else "(none)",
            tool_descriptions=self.registry.describe_tools(),
        )

        try:
            reason_response = self.llm.invoke([HumanMessage(content=reason_text)])
            plan_text = self._get_text(reason_response)
            llm_calls += 1
        except Exception as exc:
            logger.error(f"[ChemCoworker] Reasoning LLM call failed: {exc}")
            return ChemResponse(
                query=query,
                task_type=task_type,
                warnings=[f"Reasoning failed: {exc}"],
                model=self.model_name,
                elapsed_s=time.monotonic() - start,
                llm_calls=llm_calls,
            )

        if self.verbose:
            logger.info(f"[ChemCoworker] Plan text length: {len(plan_text)} chars")

        # ── Step 3: Parse plan ─────────────────────────────────────────
        plan = self.parser.parse(
            plan_text,
            known_tools=self.registry.names(),
            smiles_context=primary_smiles,
        )

        if self.verbose:
            logger.info(f"[ChemCoworker] {plan}")

        # ── Step 4: Execute tool groups ────────────────────────────────
        warnings: List[str] = []
        tool_results: Dict[str, Any] = {}

        if not plan.is_empty:
            callables = self.registry.get_callables()
            tool_results = self.executor.run_plan(plan, callables)

            # Collect any tool-level warnings
            for name, result in tool_results.items():
                if isinstance(result, dict) and not result.get("success", True):
                    warnings.append(f"Tool '{name}' failed: {result.get('error', '?')}")

        tools_called = [tc.name for tc in plan.all_tool_calls if tc.name in tool_results]

        # ── Step 5: Validate (optional revision) ──────────────────────
        # Currently lightweight — just log if hypothesis contradicted.
        # Future: add a revision LLM call here if validation fails.
        contradiction = self._check_hypothesis(plan, tool_results)
        if contradiction and self.verbose:
            logger.info(f"[ChemCoworker] Hypothesis may need revision: {contradiction}")

        # ── Step 6: Synthesize (LLM Call 2) ───────────────────────────
        tool_results_text = self._format_tool_results(tool_results)

        synth_text = SYNTHESIZE_PROMPT.format(
            query=query,
            task_type=task_type.upper(),
            hypothesis=plan.hypothesis or "(not yet identified)",
            confidence=plan.confidence,
            tool_results_text=tool_results_text,
        )

        try:
            synth_response = self.llm.invoke([HumanMessage(content=synth_text)])
            answer = self._get_text(synth_response)
            llm_calls += 1
        except Exception as exc:
            logger.error(f"[ChemCoworker] Synthesis LLM call failed: {exc}")
            answer = f"Tool results gathered but synthesis failed: {exc}"
            warnings.append(f"Synthesis failed: {exc}")

        elapsed = time.monotonic() - start

        # ── Build structured outputs ───────────────────────────────────
        structured = self._extract_structured(tool_results)

        return ChemResponse(
            query=query,
            task_type=task_type,
            hypothesis=plan.hypothesis,
            plan_rationale=plan.rationale,
            plan_text=plan_text,
            answer=answer,
            tools_called=tools_called,
            tool_results=tool_results,
            structured=structured,
            confidence=plan.confidence,
            warnings=warnings,
            model=self.model_name,
            provider=self.provider,
            elapsed_s=round(elapsed, 2),
            llm_calls=llm_calls,
        )

    def chat(
        self,
        query: str,
        history: Optional[List[Dict[str, str]]] = None,
    ) -> Tuple["ChemResponse", List[Dict[str, str]]]:
        """
        Multi-turn conversation. Maintains history as a list of
        {"role": "user"/"assistant", "content": "..."} dicts.

        Returns: (ChemResponse, updated_history)
        """
        if history is None:
            history = []

        # Prepend recent history as context for the reasoning step
        context = ""
        if history:
            recent = history[-6:]  # last 3 turns
            context = "\n".join(
                f"[{h['role'].upper()}]: {h['content']}" for h in recent
            )
            query_with_context = f"Previous conversation:\n{context}\n\n---\nNew question: {query}"
        else:
            query_with_context = query

        response = self.run(query_with_context)

        # Update history
        history = list(history) + [
            {"role": "user", "content": query},
            {"role": "assistant", "content": response.answer},
        ]

        return response, history

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get_text(self, llm_response: Any) -> str:
        """Extract text from a LangChain message response."""
        content = getattr(llm_response, "content", llm_response)
        if isinstance(content, str):
            return content
        if isinstance(content, list):
            parts = []
            for item in content:
                if isinstance(item, dict):
                    parts.append(str(item.get("text") or item.get("content", "")))
                elif isinstance(item, str):
                    parts.append(item)
            return "\n".join(parts)
        return str(content)

    def _format_tool_results(self, results: Dict[str, Any]) -> str:
        """Format tool results as readable text for the synthesis prompt."""
        if not results:
            return "(No tools were called — answering from LLM chemistry knowledge)"

        lines = []
        for name, result in results.items():
            lines.append(f"\n--- {name} ---")
            if isinstance(result, dict):
                # Omit the success flag and show rest
                display = {k: v for k, v in result.items() if k != "success"}
                try:
                    lines.append(json.dumps(display, indent=2, default=str)[:2000])
                except Exception:
                    lines.append(str(display)[:2000])
            else:
                lines.append(str(result)[:2000])

        return "\n".join(lines)

    def _extract_structured(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract key machine-readable outputs from tool results."""
        structured: Dict[str, Any] = {}

        if "detect_reaction_type" in results:
            r = results["detect_reaction_type"]
            if isinstance(r, dict) and r.get("success"):
                structured["reaction_type"] = r.get("reaction_type")
                structured["reaction_family"] = r.get("family_label")

        if "analyze_bond_changes" in results:
            r = results["analyze_bond_changes"]
            if isinstance(r, dict) and r.get("success"):
                structured["bonds_formed"] = r.get("bonds_formed", [])
                structured["bonds_broken"] = r.get("bonds_broken", [])
                structured["key_bond_type"] = r.get("key_bond_type")

        if "recommend_conditions" in results:
            r = results["recommend_conditions"]
            if isinstance(r, dict) and r.get("success"):
                structured["conditions"] = r.get("recommendations", [])

        if "get_molecular_descriptors" in results:
            r = results["get_molecular_descriptors"]
            if isinstance(r, dict) and r.get("success"):
                structured["descriptors"] = r.get("descriptors", {})
                structured["is_drug_like"] = r.get("is_drug_like")

        if "search_reaction_types" in results:
            r = results["search_reaction_types"]
            if isinstance(r, dict) and r.get("success"):
                structured["taxonomy_matches"] = r.get("matches", [])

        return structured

    def _check_hypothesis(
        self, plan: "ExecutionPlan", results: Dict[str, Any]
    ) -> Optional[str]:
        """
        Lightweight hypothesis validation.
        Returns a string description of the contradiction, or None if OK.
        """
        if not plan.hypothesis or not results:
            return None

        # If detect_reaction_type ran and found nothing, flag it
        if "detect_reaction_type" in results:
            r = results["detect_reaction_type"]
            if isinstance(r, dict) and r.get("success"):
                rt = r.get("reaction_type")
                if not rt and plan.confidence >= 0.8:
                    return "Deterministic classifier found no reaction type despite HIGH hypothesis confidence"

        return None


def create_coworker(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    verbose: bool = False,
) -> ChemCoworker:
    """Factory function — convenience alternative to ChemCoworker(...)."""
    return ChemCoworker(provider=provider, model=model, verbose=verbose)


