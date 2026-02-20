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
from typing import Any, Callable, Dict, List, Optional, Tuple

from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)

# When initial plan confidence is below this threshold AND the plan has >1 group,
# the observe step fires: Group 0 runs first, then the LLM revises Groups 1+.
_OBSERVE_THRESHOLD = 0.75

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
        # A3 — real-time streaming
        progress_cb: Optional[Callable[[str, str, float], None]] = None,
        phase_cb: Optional[Callable[[str], None]] = None,
        # A5 — answer streaming
        pre_synth_cb: Optional[Callable[[str, float, str, List[str], bool], None]] = None,
        stream_cb: Optional[Callable[[str], None]] = None,
        # A2 — plan approval
        plan_callback: Optional[Callable[["ExecutionPlan"], "ExecutionPlan"]] = None,
        # A1 — pre/post tool hooks
        hooks: Optional[Any] = None,   # HookRegistry
    ):
        self.provider = provider or os.getenv("LLM_PROVIDER", "openai")
        self.model_name = model or os.getenv("LLM_MODEL", "o4-mini")
        self.verbose = verbose
        self.progress_cb = progress_cb
        self.phase_cb = phase_cb
        self.pre_synth_cb = pre_synth_cb
        self.stream_cb = stream_cb
        self.plan_callback = plan_callback

        self.llm = _get_llm_client(provider, model, temperature)

        from .tools import REGISTRY
        from .classifier import TaskClassifier
        from .plan import PlanParser
        from .executor import ToolExecutor

        self.registry = REGISTRY
        self.classifier = TaskClassifier()
        self.parser = PlanParser()
        self.executor = ToolExecutor(
            verbose=verbose,
            progress_cb=progress_cb,
            hooks=hooks,
        )

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

        # Route retrosynthesis queries to specialized prompts
        is_retro = (task_type == "retrosynthesis")
        if is_retro:
            from .retro_prompts import RETRO_REASON_PROMPT, RETRO_SYNTHESIZE_PROMPT
            _reason_prompt = RETRO_REASON_PROMPT
            _synth_prompt = RETRO_SYNTHESIZE_PROMPT
        else:
            _reason_prompt = REASON_PROMPT
            _synth_prompt = SYNTHESIZE_PROMPT

        reason_text = _reason_prompt.format(
            task_type=task_type.upper(),
            query=query,
            smiles_list=", ".join(smiles_list) if smiles_list else "(none)",
            tool_descriptions=self.registry.describe_tools(),
        )

        if self.phase_cb: self.phase_cb("reason_start")  # A3
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
        finally:
            if self.phase_cb: self.phase_cb("reason_done")  # A3

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

        # ── Step 3.5: Plan approval — pause for user confirmation (A2) ─
        if self.plan_callback and not plan.is_empty:
            from .plan import PlanRejected
            try:
                plan = self.plan_callback(plan)
            except PlanRejected as e:
                return ChemResponse(
                    query=query,
                    task_type=task_type,
                    answer=f"Plan cancelled: {e}",
                    hypothesis=plan.hypothesis,
                    confidence=plan.confidence,
                    model=self.model_name,
                    provider=self.provider,
                    elapsed_s=round(time.monotonic() - start, 2),
                    llm_calls=llm_calls,
                )

        # ── Step 4: Execute tool groups (observe-then-plan) ───────────
        warnings: List[str] = []
        tool_results: Dict[str, Any] = {}
        plan_revised = False
        observe_llm_text = ""
        effective_plan = plan

        if not plan.is_empty:
            callables = self.registry.get_callables()

            if plan.confidence < _OBSERVE_THRESHOLD and len(plan.groups) > 1:
                # Uncertain path: run G0 diagnostics first, then let LLM revise G1+
                if self.verbose:
                    logger.info(
                        f"[ChemCoworker] confidence={plan.confidence:.2f} < {_OBSERVE_THRESHOLD} "
                        f"— triggering observe step"
                    )

                # 4a. Execute only Group 0
                g0_plan = ExecutionPlan(
                    hypothesis=plan.hypothesis,
                    confidence=plan.confidence,
                    groups=[plan.groups[0]],
                    rationale="G0 diagnostic run",
                    raw_plan_text="",
                )
                g0_results = self.executor.run_plan(g0_plan, callables)
                tool_results.update(g0_results)

                for name, result in g0_results.items():
                    if isinstance(result, dict) and not result.get("success", True):
                        warnings.append(f"Tool '{name}' failed: {result.get('error', '?')}")

                if self.verbose:
                    logger.info(f"[ChemCoworker] G0 results: {list(g0_results.keys())}")

                # 4b. Observe step — LLM sees G0 results, produces revised plan for G1+
                revised_plan, observe_llm_text, plan_revised = self._run_observe_step(
                    query=query,
                    plan=plan,
                    g0_results=g0_results,
                    primary_smiles=primary_smiles,
                )
                llm_calls += 1

                # 4c. Execute revised Groups 1+
                if not revised_plan.is_empty:
                    g1plus_results = self.executor.run_plan(revised_plan, callables)
                    tool_results.update(g1plus_results)

                    for name, result in g1plus_results.items():
                        if isinstance(result, dict) and not result.get("success", True):
                            warnings.append(f"Tool '{name}' failed: {result.get('error', '?')}")

                effective_plan = revised_plan if not revised_plan.is_empty else plan

            else:
                # Confident path: execute all groups as before (no regression)
                tool_results = self.executor.run_plan(plan, callables)

                for name, result in tool_results.items():
                    if isinstance(result, dict) and not result.get("success", True):
                        warnings.append(f"Tool '{name}' failed: {result.get('error', '?')}")

                effective_plan = plan

        # All executed tool names in insertion order (G0 + G1+ for uncertain path)
        tools_called = list(tool_results.keys())

        # ── Step 5: Validate ───────────────────────────────────────────
        contradiction = self._check_hypothesis(effective_plan, tool_results)
        if contradiction and self.verbose:
            logger.info(f"[ChemCoworker] Hypothesis may need revision: {contradiction}")

        # ── Step 6: Synthesize (final LLM call) ───────────────────────
        tool_results_text = self._format_tool_results(tool_results)

        synth_text = _synth_prompt.format(
            query=query,
            task_type=task_type.upper(),
            hypothesis=effective_plan.hypothesis or "(not yet identified)",
            confidence=effective_plan.confidence,
            tool_results_text=tool_results_text,
            tool_descriptions=self.registry.describe_tools(),
            resource_context=self._describe_resources(),
        )

        # A5 — notify CLI of plan info so it can print hypothesis/tools before streaming
        if self.pre_synth_cb:
            self.pre_synth_cb(
                effective_plan.hypothesis or "",
                effective_plan.confidence,
                effective_plan.rationale or "",
                tools_called,
                plan_revised,
            )

        streamed = False
        if self.phase_cb: self.phase_cb("synth_start")  # A3
        try:
            if self.stream_cb:
                # A5 — stream tokens directly to CLI as they arrive
                answer_chunks: List[str] = []
                for chunk in self.llm.stream([HumanMessage(content=synth_text)]):
                    token = self._get_text(chunk)
                    if token:
                        self.stream_cb(token)
                        answer_chunks.append(token)
                answer = "".join(answer_chunks)
                streamed = True
            else:
                synth_response = self.llm.invoke([HumanMessage(content=synth_text)])
                answer = self._get_text(synth_response)
            llm_calls += 1
        except Exception as exc:
            logger.error(f"[ChemCoworker] Synthesis LLM call failed: {exc}")
            answer = f"Tool results gathered but synthesis failed: {exc}"
            warnings.append(f"Synthesis failed: {exc}")
        finally:
            if self.phase_cb: self.phase_cb("synth_done")  # A3

        elapsed = time.monotonic() - start

        # ── Build structured outputs ───────────────────────────────────
        structured = self._extract_structured(tool_results)

        return ChemResponse(
            query=query,
            task_type=task_type,
            hypothesis=effective_plan.hypothesis,
            plan_rationale=effective_plan.rationale,
            plan_text=plan_text,
            plan_revised=plan_revised,
            observe_text=observe_llm_text,
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

        if self.phase_cb: self.phase_cb("compact_start")
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
            if self.phase_cb: self.phase_cb("compact_done")

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

    def _run_observe_step(
        self,
        query: str,
        plan: "ExecutionPlan",
        g0_results: Dict[str, Any],
        primary_smiles: str,
    ) -> Tuple["ExecutionPlan", str, bool]:
        """
        Mid-pipeline observe step: show Group 0 results to the LLM and get a
        revised plan for Groups 1+.

        Returns:
            (revised_plan, observe_llm_text, revision_happened)
            On LLM failure, gracefully falls back to the original plan.groups[1:].
        """
        from .prompts import OBSERVE_PROMPT
        from langchain_core.messages import HumanMessage

        g0_results_text = self._format_tool_results(g0_results)

        # Exclude G0 tools so PlanParser silently drops any hallucinated duplicates
        _g0_tools = {"normalize_reaction", "detect_reaction_type"}
        remaining_tools = [n for n in self.registry.names() if n not in _g0_tools]

        observe_text = OBSERVE_PROMPT.format(
            query=query,
            hypothesis=plan.hypothesis or "(none)",
            initial_confidence=plan.confidence,
            g0_results_text=g0_results_text,
            tool_descriptions=self.registry.describe_tools(),
        )

        if self.phase_cb: self.phase_cb("observe_start")  # A3
        try:
            observe_response = self.llm.invoke([HumanMessage(content=observe_text)])
            observe_llm_text = self._get_text(observe_response)
        except Exception as exc:
            logger.warning(f"[ChemCoworker] Observe LLM call failed: {exc}. Using original plan.")
            if self.phase_cb: self.phase_cb("observe_done")  # A3
            fallback = ExecutionPlan(
                hypothesis=plan.hypothesis,
                confidence=plan.confidence,
                groups=plan.groups[1:],
                rationale=f"[Observe step failed: {exc}] Falling back to original plan.",
                raw_plan_text="",
            )
            return fallback, "", True

        if self.phase_cb: self.phase_cb("observe_done")  # A3

        revised_plan = self.parser.parse(
            observe_llm_text,
            known_tools=remaining_tools,
            smiles_context=primary_smiles,
        )

        if self.verbose:
            logger.info(f"[ChemCoworker] Revised plan: {revised_plan}")

        return revised_plan, observe_llm_text, True

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


