"""
ChemCoworker — General-purpose chemistry AI agent.

"Claude Code for chemistry" — combines LLM chemistry expertise with a
lightweight tool ecosystem. Follows a 3-step workflow:

  1. INTAKE   : classify task, extract SMILES
  2. REASON   : LLM reasons via native tool calls (bind_tools loop)
  3. VALIDATE : collect caveats from validator _warnings; append to answer

Uses native LangChain tool calling throughout — no classic text-plan path.
Typical LLM call count: 1–9 (reasoning iterations) + optional exhaustion-guard.
"""
from __future__ import annotations

import json
import logging
import os
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

from dotenv import load_dotenv
from .event_bus import EventBus, ChemEvent
from .workflow import WORKFLOW_REGISTRY, WorkflowDefinition

load_dotenv()

logger = logging.getLogger(__name__)

# A4 — Conversation compaction: auto-summarize history when it grows too long.
_COMPACT_THRESHOLD  = 20   # total messages before triggering compaction
_COMPACT_KEEP_RECENT = 6   # most-recent messages kept verbatim (= 3 full turns)

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
        # Phase 2 — single unified observer (replaces progress_cb/phase_cb/pre_synth_cb/stream_cb)
        event_bus: Optional[EventBus] = None,
        # A2 — plan approval (kept as explicit gate: it can raise PlanRejected)
        plan_callback: Optional[Callable[["ExecutionPlan"], "ExecutionPlan"]] = None,
        # A1 — pre/post tool hooks
        hooks: Optional[Any] = None,   # HookRegistry
    ):
        self.provider = provider or os.getenv("LLM_PROVIDER", "openai")
        self.model_name = model or os.getenv("LLM_MODEL", "o4-mini")
        self.verbose = verbose
        self.event_bus = event_bus or EventBus()
        self.plan_callback = plan_callback

        self.llm = _get_llm_client(provider, model, temperature)

        from .tools import REGISTRY
        from .classifier import TaskClassifier
        from .executor import ToolExecutor

        self.registry = REGISTRY
        self.classifier = TaskClassifier()
        self.executor = ToolExecutor(
            verbose=verbose,
            event_bus=self.event_bus,
            hooks=hooks,
            registry=self.registry,  # Phase 1: data-contract resolution + validators
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, query: str) -> "ChemResponse":
        """
        Single-turn run. Returns a ChemResponse with tool results and final answer.
        """
        from .response import ChemResponse
        start = time.monotonic()
        llm_calls = 0

        # ── Step 1: Intake ─────────────────────────────────────────────
        clf = self.classifier.classify(query)
        task_type = clf.task_type_str
        smiles_list = clf.all_smiles
        primary_smiles = clf.primary_smiles or ""

        if self.verbose:
            logger.info(f"[ChemCoworker] task={task_type} smiles={smiles_list}")

        # Route to the appropriate workflow
        workflow = WORKFLOW_REGISTRY.get_for_task(task_type)

        # ── Steps 2–4: Native tool calling ─────────────────────────────
        warnings: List[str] = []
        critic_findings: List[Any] = []
        tool_results, hypothesis, confidence, tool_warnings, llm_calls_native, answer, _messages = \
            self._run_native_tool_loop(
                query=query,
                task_type=task_type,
                smiles_list=smiles_list,
                workflow=workflow,
                primary_smiles=primary_smiles,
            )
        warnings.extend(tool_warnings)
        llm_calls += llm_calls_native

        effective_plan = ExecutionPlan(
            hypothesis=hypothesis,
            confidence=confidence,
            groups=[],
            rationale="(native tool calling)",
            raw_plan_text="",
        )

        tools_called = list(tool_results.keys())

        # ── Step 5: Collect caveats from validator _warnings ───────────
        caveats_text = self._collect_caveats(tool_results, warnings)
        if caveats_text and self.verbose:
            logger.info(f"[ChemCoworker] Caveats: {caveats_text[:100]}")

        # Append any critical caveats to the answer
        if answer and caveats_text:
            answer += f"\n\n---\n⚠ **Validation notes**: {caveats_text}"

        # ── Step 5b: Critic pass (retrosynthesis only) ─────────────────
        if workflow.critic_step and workflow.critic_step.enabled:
            critic_findings, critic_verdict, critic_calls = self._run_critic_loop(
                query=query,
                hypothesis=hypothesis,
                tool_results=tool_results,
                answer=answer,
                critic_step=workflow.critic_step,
            )
            llm_calls += critic_calls
            for f in critic_findings:
                warnings.append(f"[critic] {f.message}")
            if critic_findings:
                # ── Step 5c: Revision pass — revise the answer to address findings
                if workflow.critic_step.revision_pass:
                    revised_answer, rev_calls = self._run_revision_pass(
                        query=query,
                        original_answer=answer,
                        findings=critic_findings,
                        critic_verdict=critic_verdict,
                    )
                    llm_calls += rev_calls
                    answer = revised_answer
                # Append the critic's findings for transparency
                finding_lines = "\n".join(str(f) for f in critic_findings)
                answer += f"\n\n---\n🔍 **Critic review** (addressed above): {critic_verdict}\n{finding_lines}"

        # Evidence-based confidence aggregation replaces the placeholder native-loop confidence.
        effective_plan.confidence = self._aggregate_confidence(
            tool_results=tool_results,
            base_confidence=confidence,
            warnings=warnings,
            critic_findings=critic_findings,
        )

        # ── Step 6: Emit events and stream answer ──────────────────────
        self.event_bus.emit(
            ChemEvent.PRE_SYNTH,
            hypothesis=effective_plan.hypothesis or "",
            confidence=effective_plan.confidence,
            rationale=effective_plan.rationale or "",
            tools_called=tools_called,
        )

        streamed = False
        self.event_bus.emit(ChemEvent.PHASE_START, phase="synth")
        if answer:
            for char in answer:
                self.event_bus.emit(ChemEvent.STREAM_TOKEN, token=char)
            streamed = True
        self.event_bus.emit(ChemEvent.PHASE_END, phase="synth")

        elapsed = time.monotonic() - start

        # ── Build structured outputs ───────────────────────────────────
        structured = self._extract_structured(tool_results)

        return ChemResponse(
            query=query,
            task_type=task_type,
            hypothesis=effective_plan.hypothesis,
            plan_rationale=effective_plan.rationale,
            answer=answer,
            tools_called=tools_called,
            tool_results=tool_results,
            structured=structured,
            confidence=effective_plan.confidence,
            warnings=warnings,
            model=self.model_name,
            provider=self.provider,
            elapsed_s=round(elapsed, 2),
            llm_calls=llm_calls,
            streamed=streamed,
        )

    def set_verbose(self, verbose: bool) -> None:
        """Update verbosity across the agent and tool executor."""
        self.verbose = verbose
        self.executor.verbose = verbose

    def set_plan_callback(self, callback: Optional[Callable[["ExecutionPlan"], "ExecutionPlan"]]) -> None:
        """Update the optional plan approval callback at runtime."""
        self.plan_callback = callback

    def compact_history(self, history: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """Public wrapper for conversation history compaction."""
        return self._compact_history(history)

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

        # A4 — Compact history if it has grown too long
        compacted = False
        if len(history) >= _COMPACT_THRESHOLD:
            history = self._compact_history(history)
            compacted = True
            if self.verbose:
                logger.info(f"[ChemCoworker] History compacted to {len(history)} messages")

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
        response.compacted = compacted  # A4: tag so CLI can show notice

        # Update history
        history = list(history) + [
            {"role": "user", "content": query},
            {"role": "assistant", "content": response.answer},
        ]

        return response, history

    def _compact_history(self, history: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """
        A4 — Summarize the oldest conversation turns via one LLM call,
        keeping only the most recent _COMPACT_KEEP_RECENT messages verbatim.
        Falls back to a simple truncation on LLM failure.
        """
        from langchain_core.messages import HumanMessage

        old    = history[:-_COMPACT_KEEP_RECENT]
        recent = history[-_COMPACT_KEEP_RECENT:]

        convo = "\n".join(
            f"{m['role'].upper()}: {m['content'][:400]}" for m in old
        )
        prompt = (
            "Summarize this chemistry conversation concisely in bullet points. "
            "Keep: reaction types discussed, conditions recommended, warnings raised, "
            "conclusions reached. Omit: small talk, failed attempts, repeated information.\n\n"
            + convo
        )

        self.event_bus.emit(ChemEvent.COMPACT_START)
        try:
            resp = self.llm.invoke([HumanMessage(content=prompt)])
            summary_text = self._get_text(resp)
            summary_msg = {
                "role": "assistant",
                "content": f"[Summary of earlier conversation]\n{summary_text}",
            }
            return [summary_msg] + recent
        except Exception as exc:
            logger.warning(f"[ChemCoworker] History compaction LLM call failed: {exc}. Truncating.")
            return recent  # fallback: just keep recent turns, lose old ones
        finally:
            self.event_bus.emit(ChemEvent.COMPACT_END)

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
                structured["reaction_type"] = r.get("reaction_type_id") or r.get("reaction_type")
                structured["reaction_family"] = r.get("family_label")
                if r.get("reaction_type_metadata"):
                    structured["reaction_type_metadata"] = r.get("reaction_type_metadata")

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

    def _collect_available_context_from_tool_results(
        self,
        tool_results: Dict[str, Any],
    ) -> tuple[set[str], set[str]]:
        """
        Return (completed_tool_names, provided_keys) from successful tool results.

        `provided_keys` includes actual dict keys and registered ToolPlugin.provides keys
        present in the returned payload. This is used to enforce runtime `requires`.
        """
        completed_tools: set[str] = set()
        provided_keys: set[str] = set()
        plugins = getattr(self.registry, "_plugins", {})

        for tool_name, result in tool_results.items():
            if not isinstance(result, dict) or not result.get("success", False):
                continue
            completed_tools.add(tool_name)
            provided_keys.update(result.keys())
            plugin = plugins.get(tool_name)
            if plugin is None:
                continue
            for key in getattr(plugin, "provides", []) or []:
                if key in result:
                    provided_keys.add(key)
        return completed_tools, provided_keys

    def _partition_native_tool_calls_by_contracts(
        self,
        response_tool_calls: List[Dict[str, Any]],
        callables: Dict[str, Callable[..., Any]],
        tool_results: Dict[str, Any],
    ) -> tuple[List["ToolCall"], Dict[str, Any], List[str]]:
        """
        Partition tool calls into executable vs blocked/deferred by ToolPlugin contracts.

        Returns:
            (executable_calls, synthetic_results_for_blocked, warnings)
        """
        from .plan import ToolCall

        plugins = getattr(self.registry, "_plugins", {})
        completed_tools, provided_keys = self._collect_available_context_from_tool_results(tool_results)

        requested_names = {
            str(tc.get("name", ""))
            for tc in response_tool_calls
            if isinstance(tc, dict) and tc.get("name")
        }
        providers_by_key: Dict[str, set[str]] = {}
        for name, plugin in plugins.items():
            for key in getattr(plugin, "provides", []) or []:
                providers_by_key.setdefault(key, set()).add(name)

        executable: List[ToolCall] = []
        blocked_results: Dict[str, Any] = {}
        warnings: List[str] = []

        for tc in response_tool_calls:
            tool_name = tc.get("name")
            if tool_name not in callables:
                continue  # preserve existing behavior; caller will produce "not executed" ToolMessage

            plugin = plugins.get(tool_name)
            if plugin is None:
                executable.append(ToolCall(name=tool_name, args=dict(tc.get("args", {}))))
                continue

            missing_prereqs = [p for p in (plugin.prerequisites or []) if p not in completed_tools]
            missing_requires = [k for k in (plugin.requires or []) if k not in provided_keys]
            if not missing_prereqs and not missing_requires:
                executable.append(ToolCall(name=tool_name, args=dict(tc.get("args", {}))))
                continue

            waiting_on_tools: set[str] = set()
            for p in missing_prereqs:
                if p in requested_names:
                    waiting_on_tools.add(p)
            for key in missing_requires:
                for provider_name in providers_by_key.get(key, set()):
                    if provider_name in requested_names:
                        waiting_on_tools.add(provider_name)

            deferred = bool(waiting_on_tools)
            if deferred:
                reason = (
                    f"Deferred '{tool_name}' due to unmet tool contracts; wait for: "
                    f"{', '.join(sorted(waiting_on_tools))}."
                )
            else:
                bits: List[str] = []
                if missing_prereqs:
                    bits.append(f"missing prerequisites={missing_prereqs}")
                if missing_requires:
                    bits.append(f"missing required context keys={missing_requires}")
                reason = f"Blocked '{tool_name}' due to unmet tool contracts ({'; '.join(bits)})."

            warnings.append(reason)
            blocked_results[tool_name] = {
                "success": False,
                "error": reason,
                "tool": tool_name,
                "contract_violation": True,
                "deferred": deferred,
                "missing_prerequisites": missing_prereqs,
                "missing_requires": missing_requires,
            }

        return executable, blocked_results, warnings

    def _describe_resources(self) -> str:
        """
        Return a brief, factual description of the local databases and data
        accessible via the registered tools. Used in SYNTHESIZE_PROMPT so the
        LLM can answer meta-questions like 'what data do you have access to?'
        """
        import json as _json
        import pathlib as _pathlib

        _data_dir = _pathlib.Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data"
        lines: List[str] = []

        # Reaction taxonomy — load directly from JSON
        try:
            rt_data = _json.loads((_data_dir / "reaction_types.v4.0.json").read_text())
            n_types = len(rt_data.get("reaction_types", []))
        except Exception:
            n_types = "~60"
        lines.append(f"  • Reaction taxonomy: {n_types} named reaction types "
                     "(Suzuki, Buchwald-Hartwig, SNAr, Heck, Negishi, Chan-Lam, etc.) "
                     "with scope, mechanism, and conditions metadata")

        # Motif scope index
        try:
            mi_data = _json.loads((_data_dir / "motif_scope_index.v1.json").read_text())
            n_scope = len(mi_data.get("scope_map", {}))
            n_probes = len(mi_data.get("probe_smiles", []))
        except Exception:
            n_scope, n_probes = "~64", "~22"
        lines.append(f"  • Motif scope index: {n_scope} motif-reaction scope entries, "
                     f"{n_probes} structural probe SMILES "
                     "(maps aryl halides, boronic acids, amines, etc. to compatible reactions)")

        # Reagent registry
        try:
            from chemtools.reagent.lookup import _load_all_reagents
            n_reagents = len(_load_all_reagents())
        except Exception:
            n_reagents = "~27,000"
        lines.append(f"  • Reagent registry: {n_reagents:,} entries "
                     "(catalysts, ligands, bases, solvents, oxidants, reductants) "
                     "searchable by name, abbreviation, CAS, or role")

        # HTE conditions database
        try:
            from chemtools.recommend.hte_adapter import get_hte_reaction_keys
            keys = get_hte_reaction_keys()
            hte_info = f"{len(keys)} reaction keys" if keys else "available"
        except Exception:
            hte_info = "available"
        lines.append(f"  • HTE conditions database: {hte_info} "
                     "(high-throughput screening results; used by recommend_conditions tool)")

        # RDKit
        try:
            from rdkit import __version__ as _rdkit_ver
            rdkit_info = f"v{_rdkit_ver}"
        except Exception:
            rdkit_info = "available"
        lines.append(f"  • RDKit {rdkit_info}: molecular descriptors, fingerprints, "
                     "substructure search, SMILES parsing/canonicalization, stereochemistry")

        # Literature folder
        try:
            from .tools.literature import describe_literature_folder
            lines.append(describe_literature_folder())
        except Exception:
            pass

        return "\n".join(lines)

    def _build_native_tools(self) -> List[Any]:
        """Build LangChain StructuredTool objects from all registered ToolPlugins."""
        tools = []
        for name, plugin in self.registry._plugins.items():
            try:
                tools.append(plugin.to_langchain_tool())
            except Exception as exc:
                if self.verbose:
                    logger.warning(f"[ChemCoworker] Skipping tool '{name}' (schema build failed): {exc}")
        return tools

    def _run_native_tool_loop(
        self,
        query: str,
        task_type: str,
        smiles_list: List[str],
        workflow: WorkflowDefinition,
        primary_smiles: str,
    ) -> tuple:
        """
        Run the reasoning+execution phase via native LangChain tool calling.

        Mirrors how Claude Code works:
          - SystemMessage carries the chemistry identity + tool usage rules (static,
            loaded once per session; API can cache it)
          - HumanMessage contains only the specific task query (minimal, per request)
          - Tool descriptions live exclusively in the API schema — NOT duplicated in prompts
          - The model's final turn writes the answer directly, so SYNTHESIZE_PROMPT
            is typically not needed (saves one LLM call)

        Returns:
            (tool_results, hypothesis, confidence, warnings, llm_call_count, final_answer)
            final_answer is non-empty when the model produced a text response after
            finishing all tool calls — caller can use it directly, skipping SYNTHESIZE.
        """
        from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage
        from .plan import ToolCall

        # Use the workflow's pre-selected system prompt (no template vars, API-cacheable)
        native_system = workflow.system_prompt
        max_iterations = workflow.max_iterations

        native_tools = self._build_native_tools()
        if not native_tools:
            logger.warning("[ChemCoworker] No native tools built — falling back to knowledge-only")
            return {}, "", 0.5, [], 0, "", []

        llm_with_tools = self.llm.bind_tools(native_tools)

        # Build message thread:
        #   SystemMessage — chemistry identity + tool rules (static, API-cacheable)
        #   HumanMessage  — just the concrete task (minimal, per request)
        smiles_str = ", ".join(smiles_list) if smiles_list else "(none)"
        user_message = (
            f"Task type: {task_type.upper()}\n"
            f"Query: {query}\n"
            f"SMILES found: {smiles_str}"
        )

        messages: List[Any] = [
            SystemMessage(content=native_system),
            HumanMessage(content=user_message),
        ]
        tool_results: Dict[str, Any] = {}
        hypothesis = ""
        confidence = 0.5
        warnings: List[str] = []
        llm_call_count = 0
        final_answer = ""
        callables = self.registry.get_callables()

        self.event_bus.emit(ChemEvent.PHASE_START, phase="reason")

        for iteration in range(max_iterations):
            try:
                response = llm_with_tools.invoke(messages)
                llm_call_count += 1
            except Exception as exc:
                logger.error(f"[ChemCoworker] Native tool loop iteration {iteration} failed: {exc}")
                warnings.append(f"LLM call failed at iteration {iteration}: {exc}")
                break

            messages.append(response)

            # Extract any reasoning text from this response
            resp_text = self._get_text(response)
            if resp_text and not hypothesis:
                import re as _re
                _m = _re.search(r"hypothesis[:\s]+([^\n.]{10,200})", resp_text, _re.IGNORECASE)
                ph = _m.group(1).strip() if _m else (resp_text.strip().split(".")[0].strip() if resp_text.strip() else "")
                if ph:
                    hypothesis = ph

            # Check for tool calls
            response_tool_calls = getattr(response, "tool_calls", []) or []
            if not response_tool_calls:
                # No more tool calls — model wrote its final answer
                if resp_text:
                    final_answer = resp_text
                    if not hypothesis:
                        hypothesis = resp_text[:300]
                if self.verbose:
                    logger.info(
                        f"[ChemCoworker] Native loop finished after {iteration + 1} iteration(s); "
                        f"tools called: {list(tool_results.keys())}"
                    )
                break

            if self.verbose:
                names = [tc.get("name", "?") for tc in response_tool_calls]
                logger.info(f"[ChemCoworker] Native loop iter {iteration}: {names}")

            # Execute all tool calls in parallel
            tc_list, blocked_results, contract_warnings = self._partition_native_tool_calls_by_contracts(
                response_tool_calls=response_tool_calls,
                callables=callables,
                tool_results=tool_results,
            )
            if blocked_results:
                tool_results.update(blocked_results)
                warnings.extend(contract_warnings)

            if tc_list:
                import time as _t
                t0 = _t.monotonic()
                group_results = self.executor._run_parallel(tc_list, callables)
                tool_results.update(group_results)
                elapsed = _t.monotonic() - t0

                for name, result in group_results.items():
                    if isinstance(result, dict) and not result.get("success", True):
                        warnings.append(f"Tool '{name}' failed: {result.get('error', '?')}")

                if self.verbose:
                    logger.info(
                        f"[ChemCoworker] Native loop iter {iteration}: "
                        f"{list(group_results.keys())} done in {elapsed:.2f}s"
                    )

            # Feed results back as ToolMessages (required for structured tool calling API)
            for tc in response_tool_calls:
                tool_name = tc.get("name", "")
                tool_call_id = tc.get("id", f"{tool_name}_{iteration}")
                result = tool_results.get(tool_name, {"success": False, "error": "not executed"})
                result_str = json.dumps(result, default=str)[:3000]
                messages.append(ToolMessage(
                    content=result_str,
                    tool_call_id=tool_call_id,
                ))

        self.event_bus.emit(ChemEvent.PHASE_END, phase="reason")

        # Exhaustion guard: if the loop hit max_iterations without producing a
        # final_answer (model kept calling tools without writing a conclusion),
        # inject warnings and make one closing call to get a proper answer.
        if not final_answer:
            tool_warnings_text = "\n".join(
                f"• {w}" for r in tool_results.values()
                if isinstance(r, dict)
                for w in r.get("_warnings", [])
            )
            if tool_warnings_text:
                messages.append(HumanMessage(
                    content=f"Validation notes to address in your answer:\n{tool_warnings_text}"
                ))
            messages.append(HumanMessage(
                content="You have gathered sufficient tool evidence. Write your expert answer now."
            ))
            try:
                closing_response = llm_with_tools.invoke(messages)
                llm_call_count += 1
                final_answer = self._get_text(closing_response) or ""
                messages.append(closing_response)
                if not hypothesis:
                    hypothesis = final_answer[:300]
            except Exception as exc:
                logger.error(f"[ChemCoworker] Exhaustion-guard closing call failed: {exc}")
                warnings.append(f"Closing call failed after loop exhaustion: {exc}")

        return tool_results, hypothesis, confidence, warnings, llm_call_count, final_answer, messages

    def _collect_caveats(
        self, tool_results: Dict[str, Any], existing_warnings: List[str]
    ) -> str:
        """
        Collect caveats from tool validator _warnings and existing warnings.
        Returns a deduplicated, newline-joined string (empty string if none).
        Replaces the classic _check_hypothesis() method.
        """
        parts: List[str] = []
        for tool_name, result in tool_results.items():
            if not isinstance(result, dict):
                continue
            for w in result.get("_warnings", []):
                parts.append(f"⚠ {tool_name}: {w}")
        for w in existing_warnings:
            parts.append(f"• {w}")
        # Deduplicate while preserving order
        return "\n".join(dict.fromkeys(parts))

    def _aggregate_confidence(
        self,
        tool_results: Dict[str, Any],
        base_confidence: float = 0.5,
        warnings: Optional[List[str]] = None,
        critic_findings: Optional[List[Any]] = None,
    ) -> float:
        """
        Estimate answer confidence from chemistry evidence, support, and caveats.

        Heuristic (deterministic) aggregator; intended to be calibrated later.
        """
        warnings = warnings or []
        critic_findings = critic_findings or []

        def _clip(x: float) -> float:
            return max(0.05, min(0.99, x))

        def _coerce_num(x: Any) -> Optional[float]:
            try:
                if x is None:
                    return None
                if isinstance(x, (int, float)):
                    return float(x)
                s = str(x).strip()
                if not s:
                    return None
                if s.endswith("%"):
                    val = float(s[:-1]) / 100.0
                else:
                    val = float(s)
                if val > 1.0 and val <= 100.0:
                    val = val / 100.0
                return val
            except Exception:
                return None

        score = float(base_confidence or 0.5)

        # Reaction typing confidence and taxonomy grounding
        rtype = tool_results.get("detect_reaction_type")
        if isinstance(rtype, dict) and rtype.get("success"):
            det_conf = _coerce_num(rtype.get("confidence"))
            if det_conf is not None:
                score = 0.55 * score + 0.45 * det_conf
            rt_id = str(rtype.get("reaction_type_id") or rtype.get("reaction_type") or "").strip().lower()
            meta = rtype.get("reaction_type_metadata") or {}
            if rt_id and rt_id != "unknown":
                score += 0.08
            else:
                score -= 0.18
            if meta.get("category"):
                score += 0.04
            if rtype.get("reacted_motifs") or rtype.get("formed_motifs"):
                score += 0.03

        # Bond-change / mapping support
        bchg = tool_results.get("analyze_bond_changes")
        if isinstance(bchg, dict) and bchg.get("success"):
            map_conf = _coerce_num(bchg.get("mapping_confidence"))
            if map_conf is not None:
                score = 0.75 * score + 0.25 * map_conf
            formed = bchg.get("bonds_formed") or []
            if formed:
                score += 0.05
            else:
                score -= 0.15
            key_bond = str(bchg.get("key_bond_type") or "").strip().lower()
            if key_bond and key_bond != "unknown":
                score += 0.04
            else:
                score -= 0.10

        # Condition recommendation support quality (if relevant)
        cond = tool_results.get("recommend_conditions")
        if isinstance(cond, dict) and cond.get("success"):
            recs = cond.get("recommendations") or []
            if recs:
                score += 0.04
                top = recs[0] if isinstance(recs[0], dict) else {}
                n_exp = int(top.get("num_experiments", 0) or 0)
                top_conf = _coerce_num(top.get("confidence"))
                if n_exp > 0:
                    score += 0.06 if n_exp >= 5 else 0.03
                if top_conf is not None:
                    score = 0.85 * score + 0.15 * top_conf
            else:
                score -= 0.10

        # Penalize contract violations and warnings
        contract_violations = 0
        for result in tool_results.values():
            if isinstance(result, dict) and result.get("contract_violation"):
                contract_violations += 1
        score -= 0.10 * contract_violations

        warning_penalty = 0.0
        for w in warnings:
            text = str(w).lower()
            if "contract" in text:
                warning_penalty += 0.05
            elif "failed" in text or "unknown" in text:
                warning_penalty += 0.03
            else:
                warning_penalty += 0.015
        score -= min(0.25, warning_penalty)

        # Critic findings reduce confidence by severity
        severity_penalty = 0.0
        for f in critic_findings:
            sev = str(getattr(getattr(f, "severity", None), "value", getattr(f, "severity", ""))).lower()
            if sev == "critical":
                severity_penalty += 0.20
            elif sev == "warning":
                severity_penalty += 0.08
            elif sev:
                severity_penalty += 0.03
        score -= min(0.35, severity_penalty)

        return round(_clip(score), 3)

    def _run_critic_loop(
        self,
        query: str,
        hypothesis: str,
        tool_results: Dict[str, Any],
        answer: str,
        critic_step: Any,
    ) -> tuple:  # (List[Finding], str, int)
        """
        Phase 6 — Run the adversarial critic pass for retrosynthesis workflows.

        Calls CriticAgent.review() with the main loop's outputs and returns
        structured findings plus a one-sentence verdict.

        Returns:
            (findings, verdict, llm_call_count)
            findings is a List[Finding]; verdict is a str; llm_call_count is int.
            Returns ([], "", 0) on any failure so the main answer is unaffected.
        """
        from .critic import CriticAgent, Severity

        self.event_bus.emit(ChemEvent.PHASE_START, phase="critic")
        try:
            critic = CriticAgent(self.llm)
            min_sev = Severity(critic_step.min_severity)
            findings, verdict = critic.review(
                query=query,
                hypothesis=hypothesis,
                tool_results=tool_results,
                answer=answer,
                max_findings=critic_step.max_findings,
                min_severity=min_sev,
            )
            if self.verbose:
                logger.info(
                    f"[ChemCoworker] Critic: {len(findings)} finding(s); verdict={verdict[:80]!r}"
                )
            return findings, verdict, 1
        except Exception as exc:
            logger.warning(f"[ChemCoworker] Critic pass failed: {exc}")
            return [], "", 0
        finally:
            self.event_bus.emit(ChemEvent.PHASE_END, phase="critic")


    def _run_revision_pass(
        self,
        query: str,
        original_answer: str,
        findings: List[Any],  # List[Finding]
        critic_verdict: str,
    ) -> tuple:  # (revised_answer: str, llm_call_count: int)
        """
        Phase 6b — Follow-up LLM call that revises the main answer to address
        the critic's findings.

        Returns:
            (revised_answer, llm_call_count)
            On any failure returns (original_answer, 0) so the pipeline is
            unaffected.
        """
        from langchain_core.messages import SystemMessage, HumanMessage

        finding_lines = "\n".join(
            f"- [{f.severity.value.upper()}] {f.message}"
            + (f"\n  Suggestion: {f.suggestion}" if f.suggestion else "")
            for f in findings
        )

        revision_prompt = (
            f"You previously proposed the following synthesis answer:\n\n"
            f"{original_answer}\n\n"
            f"An adversarial chemistry reviewer has identified the following issues "
            f"(verdict: {critic_verdict}):\n\n"
            f"{finding_lines}\n\n"
            f"Please produce a revised answer that directly addresses each issue above. "
            f"Keep the overall structure of your original answer; only modify the parts "
            f"that are affected by the critic's findings. "
            f"Do not repeat the critic's findings verbatim — integrate the corrections "
            f"into the answer naturally."
        )

        self.event_bus.emit(ChemEvent.PHASE_START, phase="revision")
        try:
            messages = [
                SystemMessage(content=(
                    "You are a chemistry expert revising a synthesis proposal based on "
                    "peer-review feedback. Produce a corrected, self-contained answer."
                )),
                HumanMessage(content=revision_prompt),
            ]
            response = self.llm.invoke(messages)
            content = getattr(response, "content", response)
            revised = content if isinstance(content, str) else str(content)
            if self.verbose:
                logger.info(
                    f"[ChemCoworker] Revision pass complete "
                    f"({len(revised)} chars → was {len(original_answer)} chars)"
                )
            return revised, 1
        except Exception as exc:
            logger.warning(f"[ChemCoworker] Revision pass failed: {exc}")
            return original_answer, 0
        finally:
            self.event_bus.emit(ChemEvent.PHASE_END, phase="revision")


def create_coworker(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    verbose: bool = False,
) -> ChemCoworker:
    """Factory function — convenience alternative to ChemCoworker(...)."""
    return ChemCoworker(provider=provider, model=model, verbose=verbose)


