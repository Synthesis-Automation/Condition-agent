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

import copy
import inspect
import json
import logging
import os
import re
import threading
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

from dotenv import load_dotenv
from .event_bus import EventBus, ChemEvent
from .workflow import WORKFLOW_REGISTRY, WorkflowDefinition
from .skills import (
    SkillRecord,
    build_default_skill_registry,
    format_skill_instruction_block,
)

load_dotenv()

logger = logging.getLogger(__name__)
_TOOL_CALL_RESULTS_META_KEY = "_tool_call_results_by_id"
_REACTION_SMILES_IN_TEXT_RE = re.compile(
    r"([A-Za-z0-9@+\-\[\]\(\)\\\/%=#$:.*]+(?:\.[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$:.*]+)*)\s*>>\s*"
    r"([A-Za-z0-9@+\-\[\]\(\)\\\/%=#$:.*]+(?:\.[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$:.*]+)*)"
)

# A4 — Conversation compaction: auto-summarize history when it grows too long.
_COMPACT_THRESHOLD  = 20   # total messages before triggering compaction
_COMPACT_KEEP_RECENT = 6   # most-recent messages kept verbatim (= 3 full turns)
_MAX_REPEATED_TOOL_SIGNATURES_WITHOUT_PROGRESS = 2


def _conditions_cache_key(
    *,
    top_k: int,
    source_group: str = "",
    reaction_key_only: bool = False,
    use_spectator_groups: bool = True,
    prefer_mixfp_for_similarity: bool = False,
    similarity_mixfp_weight: float = 0.3,
) -> str:
    return (
        f"top_k={int(top_k)}"
        f"|source_group={str(source_group or '').strip().lower()}"
        f"|reaction_key_only={1 if reaction_key_only else 0}"
        f"|use_spectator_groups={1 if use_spectator_groups else 0}"
        f"|prefer_mixfp_for_similarity={1 if prefer_mixfp_for_similarity else 0}"
        f"|similarity_mixfp_weight={float(similarity_mixfp_weight):.4f}"
    )

# Local imports (no circular dependency — these modules don't import agent.py)
from .response import ChemResponse  # noqa: E402
from .plan import ExecutionPlan     # noqa: E402


@dataclass
class ReactionContext:
    """Per-reaction cached chemistry state for a single ChemCoworker run."""
    input_reaction_smiles: str
    normalized_reaction_smiles: str
    normalization_result: Optional[Dict[str, Any]] = None
    featurization_result: Optional[Dict[str, Any]] = None
    bond_change_analysis_raw: Optional[Any] = None
    bond_change_recommended: Optional[Dict[str, Any]] = None
    detect_reaction_type_result: Optional[Dict[str, Any]] = None
    analyze_bond_changes_result: Optional[Dict[str, Any]] = None
    conditions_results_by_top_k: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    reaction_type_candidates: List[Any] = field(default_factory=list)
    motif_evidence: Dict[str, Any] = field(default_factory=dict)
    fg_profile: Dict[str, Any] = field(default_factory=dict)
    _lock: Any = field(default_factory=threading.RLock, repr=False)


@dataclass
class ChemistryRunState:
    """Ephemeral per-query cache shared across tool calls in a single run."""
    reaction_contexts: Dict[str, ReactionContext] = field(default_factory=dict)
    functional_group_results: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    descriptor_results: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    _lock: Any = field(default_factory=threading.RLock, repr=False)


@dataclass
class _ToolRuntimeContext:
    """Formal per-run runtime context injected into tool calls via contextvar."""
    agent: Any
    chemistry_state: ChemistryRunState

    def normalize_reaction(self, smiles: str) -> Optional[Dict[str, Any]]:
        if not self.agent._is_reaction_smiles_input(smiles):
            return None
        rxn_ctx = self.agent._get_or_create_reaction_context(self.chemistry_state, smiles)
        return self.agent._copy_result(self.agent._hydrate_reaction_normalization(rxn_ctx))

    def detect_reaction_type(self, reaction_smiles: str) -> Dict[str, Any]:
        rxn_ctx = self.agent._get_or_create_reaction_context(self.chemistry_state, reaction_smiles)
        return self.agent._detect_reaction_type_from_context(rxn_ctx)

    def analyze_bond_changes(self, reaction_smiles: str) -> Dict[str, Any]:
        rxn_ctx = self.agent._get_or_create_reaction_context(self.chemistry_state, reaction_smiles)
        return self.agent._analyze_bond_changes_from_context(rxn_ctx)

    def get_cached_conditions(
        self,
        reaction_smiles: str,
        top_k: int,
        cache_key: str = "",
    ) -> Optional[Dict[str, Any]]:
        rxn_ctx = self.agent._get_or_create_reaction_context(self.chemistry_state, reaction_smiles)
        key = str(cache_key or "").strip() or _conditions_cache_key(top_k=int(top_k))
        with rxn_ctx._lock:
            cached = rxn_ctx.conditions_results_by_top_k.get(key)
            if cached is None:
                # Legacy fallback for pre-keyed cache entries within the same run.
                cached = rxn_ctx.conditions_results_by_top_k.get(str(int(top_k)))
            if cached is None:
                cached = rxn_ctx.conditions_results_by_top_k.get(int(top_k))
        return self.agent._copy_result(cached) if cached is not None else None

    def set_cached_conditions(
        self,
        reaction_smiles: str,
        top_k: int,
        result: Dict[str, Any],
        cache_key: str = "",
    ) -> None:
        rxn_ctx = self.agent._get_or_create_reaction_context(self.chemistry_state, reaction_smiles)
        key = str(cache_key or "").strip() or _conditions_cache_key(top_k=int(top_k))
        with rxn_ctx._lock:
            rxn_ctx.conditions_results_by_top_k[key] = self.agent._copy_result(result)

    def get_cached_molecule_result(self, cache_name: str, smiles: str) -> Optional[Dict[str, Any]]:
        key = str(smiles or "").strip()
        with self.chemistry_state._lock:
            cache = (
                self.chemistry_state.functional_group_results
                if cache_name == "functional_groups"
                else self.chemistry_state.descriptor_results
            )
            cached = cache.get(key)
        return self.agent._copy_result(cached) if cached is not None else None

    def set_cached_molecule_result(self, cache_name: str, smiles: str, result: Dict[str, Any]) -> None:
        key = str(smiles or "").strip()
        with self.chemistry_state._lock:
            cache = (
                self.chemistry_state.functional_group_results
                if cache_name == "functional_groups"
                else self.chemistry_state.descriptor_results
            )
            cache[key] = self.agent._copy_result(result)


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
        self.skill_registry = build_default_skill_registry()
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

    def _new_chemistry_run_state(self) -> ChemistryRunState:
        """Create an ephemeral per-run chemistry cache shared across tool calls."""
        return ChemistryRunState()

    def _new_tool_runtime_context(self, run_state: ChemistryRunState) -> Any:
        """Create a formal per-run runtime context object for tool calls."""
        return _ToolRuntimeContext(agent=self, chemistry_state=run_state)

    def _get_or_create_reaction_context(
        self,
        run_state: ChemistryRunState,
        reaction_smiles: str,
    ) -> ReactionContext:
        from .tools._helpers import _clean_rxn_smiles

        key = _clean_rxn_smiles(reaction_smiles or "")
        with run_state._lock:
            ctx = run_state.reaction_contexts.get(key)
            if ctx is None:
                ctx = ReactionContext(
                    input_reaction_smiles=reaction_smiles,
                    normalized_reaction_smiles=key,
                )
                run_state.reaction_contexts[key] = ctx
            return ctx

    def _is_reaction_smiles_input(self, value: Any) -> bool:
        if not isinstance(value, str):
            return False
        return ">>" in value or (">" in value and value.count(">") >= 1)

    def _copy_result(self, value: Any) -> Any:
        try:
            return copy.deepcopy(value)
        except Exception:
            return value

    def _hydrate_reaction_normalization(self, ctx: ReactionContext) -> Optional[Dict[str, Any]]:
        from .tools._helpers import _error, _success
        from chemtools.featurizers.analysis.smiles import normalize_reaction

        with ctx._lock:
            if ctx.normalization_result is not None:
                return ctx.normalization_result
            try:
                result = normalize_reaction(ctx.normalized_reaction_smiles)
                if "error" in result:
                    ctx.normalization_result = _error(f"Invalid reaction SMILES: {result['error']}")
                else:
                    ctx.normalization_result = _success({
                        "is_reaction": True,
                        "input_smiles": ctx.normalized_reaction_smiles,
                        "reactants": result.get("reactants", []),
                        "agents": result.get("agents", []),
                        "products": result.get("products", []),
                        "warnings": result.get("warnings", []),
                    })
            except Exception as exc:
                ctx.normalization_result = _error(f"SMILES normalization failed: {exc}")
            return ctx.normalization_result

    def _hydrate_reaction_featurization(self, ctx: ReactionContext) -> Optional[Dict[str, Any]]:
        from chemtools.featurizers.unified import featurize_reaction

        with ctx._lock:
            if ctx.featurization_result is not None:
                return ctx.featurization_result
            result = featurize_reaction(ctx.normalized_reaction_smiles)
            ctx.featurization_result = result if isinstance(result, dict) else None
            detection = ctx.featurization_result.get("detection", {}) if ctx.featurization_result else {}
            evidence = detection.get("evidence", {}) if isinstance(detection, dict) else {}
            candidates = detection.get("candidates") if isinstance(detection, dict) else None
            if isinstance(candidates, list):
                ctx.reaction_type_candidates = candidates
            if isinstance(evidence, dict):
                ctx.motif_evidence = {
                    "reacted_motifs": list(evidence.get("reacted_motifs", []) or []),
                    "formed_motifs": list(evidence.get("formed_motifs", []) or []),
                }
            return ctx.featurization_result

    def _hydrate_bond_change_analysis(self, ctx: ReactionContext) -> Optional[Dict[str, Any]]:
        from chemtools._atom_mapping import analyze_bond_changes_hybrid

        with ctx._lock:
            if ctx.bond_change_recommended is not None:
                return ctx.bond_change_recommended
            result = analyze_bond_changes_hybrid(ctx.normalized_reaction_smiles)
            ctx.bond_change_analysis_raw = result
            if isinstance(result, dict) and "recommended_result" in result:
                ctx.bond_change_recommended = result["recommended_result"]
            elif isinstance(result, dict):
                ctx.bond_change_recommended = result
            else:
                ctx.bond_change_recommended = None
            return ctx.bond_change_recommended

    def _hydrate_reaction_fg_profile(self, ctx: ReactionContext) -> Dict[str, Any]:
        """Lazily compute component-level functional group profile for a reaction."""
        with ctx._lock:
            if ctx.fg_profile:
                return ctx.fg_profile
        norm = self._hydrate_reaction_normalization(ctx)
        if not isinstance(norm, dict) or not norm.get("success"):
            return {}

        try:
            from chemtools.util.functional_groups import get_functional_groups, get_group_categories
        except Exception:
            return {}

        reactants = [r for r in (norm.get("reactants") or []) if isinstance(r, str) and r]
        products = [p for p in (norm.get("products") or []) if isinstance(p, str) and p]

        def _profile_many(items: List[str]) -> List[Dict[str, Any]]:
            prof: List[Dict[str, Any]] = []
            for smi in items:
                try:
                    groups = get_functional_groups(smi)
                    cats = get_group_categories(smi)
                    prof.append({
                        "smiles": smi,
                        "detected_groups": groups,
                        "categories": {k: v for k, v in (cats or {}).items() if v},
                    })
                except Exception:
                    prof.append({"smiles": smi, "detected_groups": [], "categories": {}})
            return prof

        fg = {"reactants": _profile_many(reactants), "products": _profile_many(products)}
        with ctx._lock:
            if not ctx.fg_profile:
                ctx.fg_profile = fg
            return ctx.fg_profile

    def _detect_reaction_type_from_context(self, ctx: ReactionContext) -> Dict[str, Any]:
        from .tools._helpers import _error, _success
        from chemtools.taxonomy.reaction_catalog import get_reaction_type, resolve_reaction_type

        with ctx._lock:
            if ctx.detect_reaction_type_result is not None:
                return self._copy_result(ctx.detect_reaction_type_result)

        try:
            result = self._hydrate_reaction_featurization(ctx)
            if not result:
                out = _error("Reaction featurization returned no result")
            else:
                detection = result.get("detection", {})
                validation = detection.get("validation", {}) if isinstance(detection, dict) else {}
                evidence = detection.get("evidence", {}) if isinstance(detection, dict) else {}

                reaction_type_raw = result.get("reaction_type") or validation.get("validated_detection")
                reaction_type_id = resolve_reaction_type(str(reaction_type_raw)) if reaction_type_raw else None
                reaction_type = reaction_type_id or reaction_type_raw
                confidence = (
                    validation.get("validation_confidence")
                    or result.get("confidence")
                    or 0.0
                )
                rt_def = get_reaction_type(str(reaction_type)) if reaction_type else None
                family_label = (
                    rt_def.name if rt_def else (reaction_type.replace("_", " ") if reaction_type else "")
                )
                taxonomy_metadata = {
                    "id": getattr(rt_def, "id", reaction_type_id or reaction_type or ""),
                    "name": getattr(rt_def, "name", family_label or ""),
                    "category": getattr(rt_def, "category", ""),
                    "aliases": list(getattr(rt_def, "aliases", []) or []),
                    "has_constraints": bool(getattr(rt_def, "constraints", None)),
                } if reaction_type else {}

                reacted_motifs = evidence.get("reacted_motifs", []) if isinstance(evidence, dict) else []
                formed_motifs = evidence.get("formed_motifs", []) if isinstance(evidence, dict) else []

                out = _success({
                    "reaction_smiles": ctx.normalized_reaction_smiles,
                    "reaction_type": reaction_type,
                    "reaction_type_id": reaction_type_id or reaction_type,
                    "family_label": family_label,
                    "confidence": float(confidence),
                    "reacted_motifs": reacted_motifs,
                    "formed_motifs": formed_motifs,
                    "reaction_key": result.get("reaction_key"),
                    "reaction_type_metadata": taxonomy_metadata,
                })
        except Exception as exc:
            out = _error(f"Reaction type detection failed: {exc}")

        with ctx._lock:
            ctx.detect_reaction_type_result = self._copy_result(out)
            # Update context summary fields for downstream consistency
            if isinstance(out, dict) and out.get("success"):
                ctx.motif_evidence = {
                    "reacted_motifs": list(out.get("reacted_motifs", []) or []),
                    "formed_motifs": list(out.get("formed_motifs", []) or []),
                }
        return self._copy_result(out)

    def _analyze_bond_changes_from_context(self, ctx: ReactionContext) -> Dict[str, Any]:
        from .tools._helpers import _error, _success, _to_jsonable
        from .tools.chemistry import _infer_key_bond_type

        with ctx._lock:
            if ctx.analyze_bond_changes_result is not None:
                return self._copy_result(ctx.analyze_bond_changes_result)

        try:
            rec = self._hydrate_bond_change_analysis(ctx)
            if not rec:
                out = _error("Bond change analysis returned no result")
            else:
                broken = rec.get("broken_bonds") or rec.get("bonds_broken", [])
                formed = rec.get("formed_bonds") or rec.get("bonds_formed", [])
                leaving = rec.get("leaving_groups", [])
                key_bond = _infer_key_bond_type(broken, leaving)
                out = _success({
                    "reaction_smiles": ctx.normalized_reaction_smiles,
                    "bonds_broken": _to_jsonable(broken),
                    "bonds_formed": _to_jsonable(formed),
                    "key_bond_type": key_bond,
                    "leaving_groups": _to_jsonable(leaving),
                    "mapping_confidence": rec.get("confidence", rec.get("mapping_confidence", "")),
                })
        except Exception as exc:
            out = _error(f"Bond change analysis failed: {exc}")

        with ctx._lock:
            ctx.analyze_bond_changes_result = self._copy_result(out)
        return self._copy_result(out)

    def _build_context_aware_callables(
        self,
        base_callables: Dict[str, Callable[..., Any]],
        run_state: ChemistryRunState,
    ) -> Dict[str, Callable[..., Any]]:
        """Wrap chemistry-heavy tools so they share per-run cached reaction state."""
        wrapped = dict(base_callables)

        def _normalize_wrapper(*, smiles: str) -> Any:
            if not self._is_reaction_smiles_input(smiles):
                return base_callables["normalize_reaction"](smiles=smiles)
            ctx = self._get_or_create_reaction_context(run_state, smiles)
            result = self._hydrate_reaction_normalization(ctx)
            return self._copy_result(result)

        def _detect_wrapper(*, reaction_smiles: str) -> Any:
            ctx = self._get_or_create_reaction_context(run_state, reaction_smiles)
            return self._detect_reaction_type_from_context(ctx)

        def _bond_wrapper(*, reaction_smiles: str) -> Any:
            ctx = self._get_or_create_reaction_context(run_state, reaction_smiles)
            return self._analyze_bond_changes_from_context(ctx)

        def _conditions_wrapper(
            *,
            reaction_smiles: str,
            top_k: int = 5,
            source_group: str = "",
            reaction_key_only: bool = False,
            use_spectator_groups: bool = True,
            prefer_mixfp_for_similarity: bool = False,
            similarity_mixfp_weight: float = 0.3,
        ) -> Any:
            ctx = self._get_or_create_reaction_context(run_state, reaction_smiles)
            tk = int(top_k)
            cache_key = _conditions_cache_key(
                top_k=tk,
                source_group=source_group,
                reaction_key_only=reaction_key_only,
                use_spectator_groups=use_spectator_groups,
                prefer_mixfp_for_similarity=prefer_mixfp_for_similarity,
                similarity_mixfp_weight=similarity_mixfp_weight,
            )
            with ctx._lock:
                if cache_key in ctx.conditions_results_by_top_k:
                    return self._copy_result(ctx.conditions_results_by_top_k[cache_key])
                legacy_key = str(tk)
                if legacy_key in ctx.conditions_results_by_top_k:
                    return self._copy_result(ctx.conditions_results_by_top_k[legacy_key])
                if tk in ctx.conditions_results_by_top_k:
                    return self._copy_result(ctx.conditions_results_by_top_k[tk])
            call_kwargs = {
                "reaction_smiles": reaction_smiles,
                "top_k": top_k,
                "source_group": source_group,
                "reaction_key_only": reaction_key_only,
                "use_spectator_groups": use_spectator_groups,
                "prefer_mixfp_for_similarity": prefer_mixfp_for_similarity,
                "similarity_mixfp_weight": similarity_mixfp_weight,
            }
            try:
                sig = inspect.signature(base_callables["recommend_conditions"])
                accepts_var_kwargs = any(
                    param.kind == inspect.Parameter.VAR_KEYWORD
                    for param in sig.parameters.values()
                )
                if not accepts_var_kwargs:
                    call_kwargs = {
                        key: value
                        for key, value in call_kwargs.items()
                        if key in sig.parameters
                    }
            except Exception:
                pass
            result = base_callables["recommend_conditions"](**call_kwargs)
            with ctx._lock:
                ctx.conditions_results_by_top_k[cache_key] = self._copy_result(result)
            return self._copy_result(result)

        def _inspect_fg_wrapper(*, smiles: str) -> Any:
            key = str(smiles or "").strip()
            with run_state._lock:
                cached = run_state.functional_group_results.get(key)
            if cached is not None:
                return self._copy_result(cached)
            if self._is_reaction_smiles_input(key):
                ctx = self._get_or_create_reaction_context(run_state, key)
                self._hydrate_reaction_fg_profile(ctx)
            result = base_callables["inspect_functional_groups"](smiles=smiles)
            with run_state._lock:
                run_state.functional_group_results[key] = self._copy_result(result)
            return self._copy_result(result)

        def _descriptor_wrapper(*, smiles: str) -> Any:
            key = str(smiles or "").strip()
            with run_state._lock:
                cached = run_state.descriptor_results.get(key)
            if cached is not None:
                return self._copy_result(cached)
            result = base_callables["get_molecular_descriptors"](smiles=smiles)
            with run_state._lock:
                run_state.descriptor_results[key] = self._copy_result(result)
            return self._copy_result(result)

        if "normalize_reaction" in base_callables:
            wrapped["normalize_reaction"] = _normalize_wrapper
        if "detect_reaction_type" in base_callables:
            wrapped["detect_reaction_type"] = _detect_wrapper
        if "analyze_bond_changes" in base_callables:
            wrapped["analyze_bond_changes"] = _bond_wrapper
        if "recommend_conditions" in base_callables:
            wrapped["recommend_conditions"] = _conditions_wrapper
        if "inspect_functional_groups" in base_callables:
            wrapped["inspect_functional_groups"] = _inspect_fg_wrapper
        if "get_molecular_descriptors" in base_callables:
            wrapped["get_molecular_descriptors"] = _descriptor_wrapper
        return wrapped

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
        active_skill_records = self._resolve_active_skill_records(
            query=query,
            task_type=task_type,
            workflow=workflow,
            smiles_present=bool(smiles_list),
        )

        # ── Steps 2–4: Native tool calling ─────────────────────────────
        warnings: List[str] = []
        critic_findings: List[Any] = []
        chemistry_state = self._new_chemistry_run_state()
        tool_results, hypothesis, confidence, tool_warnings, llm_calls_native, answer, _messages, reason_tokens = \
            self._run_native_tool_loop(
                query=query,
                task_type=task_type,
                smiles_list=smiles_list,
                workflow=workflow,
                primary_smiles=primary_smiles,
                chemistry_state=chemistry_state,
                active_skill_records=active_skill_records,
            )
        warnings.extend(tool_warnings)
        llm_calls += llm_calls_native
        token_sections: List[Dict[str, Any]] = []
        if reason_tokens.get("llm_calls", 0) or reason_tokens.get("total_tokens", 0):
            token_sections.append(reason_tokens)

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
            critic_findings, critic_verdict, critic_calls, critic_tokens = self._run_critic_loop(
                query=query,
                hypothesis=hypothesis,
                tool_results=tool_results,
                answer=answer,
                critic_step=workflow.critic_step,
            )
            llm_calls += critic_calls
            if critic_tokens.get("llm_calls", 0) or critic_tokens.get("total_tokens", 0):
                token_sections.append(critic_tokens)
            for f in critic_findings:
                warnings.append(f"[critic] {f.message}")
            if critic_findings:
                # ── Step 5c: Revision pass — revise the answer to address findings
                if workflow.critic_step.revision_pass:
                    revised_answer, rev_calls, revision_tokens = self._run_revision_pass(
                        query=query,
                        original_answer=answer,
                        findings=critic_findings,
                        critic_verdict=critic_verdict,
                    )
                    llm_calls += rev_calls
                    answer = revised_answer
                    if revision_tokens.get("llm_calls", 0) or revision_tokens.get("total_tokens", 0):
                        token_sections.append(revision_tokens)
                # Append the critic's findings for transparency
                finding_lines = "\n".join(str(f) for f in critic_findings)
                answer += f"\n\n---\n🔍 **Critic review** (addressed above): {critic_verdict}\n{finding_lines}"

        post_verification_penalty = 0.0
        if answer:
            raw_answer_for_repair = str(answer)
            answer, verification_warnings, post_verification_penalty, gate_report = self._apply_output_verification_gate(
                answer=answer,
                tool_results=tool_results,
                task_type=task_type,
                active_skill_records=active_skill_records,
            )
            (
                answer,
                verification_warnings,
                post_verification_penalty,
                repair_warnings,
                repair_calls,
                repair_tokens,
            ) = self._attempt_auto_repair_after_verification(
                query=query,
                raw_answer=raw_answer_for_repair,
                verified_answer=answer,
                verification_warnings=verification_warnings,
                verification_penalty=post_verification_penalty,
                verification_report=gate_report,
                tool_results=tool_results,
                task_type=task_type,
                active_skill_records=active_skill_records,
            )
            warnings.extend(verification_warnings)
            warnings.extend(repair_warnings)
            llm_calls += repair_calls
            if repair_tokens.get("llm_calls", 0) or repair_tokens.get("total_tokens", 0):
                token_sections.append(repair_tokens)

        # Evidence-based confidence aggregation replaces the placeholder native-loop confidence.
        effective_plan.confidence = self._aggregate_confidence(
            tool_results=tool_results,
            base_confidence=confidence,
            warnings=warnings,
            critic_findings=critic_findings,
        )
        if post_verification_penalty > 0.0:
            effective_plan.confidence = round(
                max(0.05, min(0.99, float(effective_plan.confidence) - float(post_verification_penalty))),
                3,
            )
        perf_summary = self._build_performance_summary(tool_results, task_type=task_type)
        if answer and perf_summary:
            answer += f"\n\n---\n{perf_summary}"

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
        prompt_tokens = sum(int(section.get("prompt_tokens", 0) or 0) for section in token_sections)
        completion_tokens = sum(int(section.get("completion_tokens", 0) or 0) for section in token_sections)
        total_tokens = sum(int(section.get("total_tokens", 0) or 0) for section in token_sections)

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
            prompt_tokens=prompt_tokens,
            completion_tokens=completion_tokens,
            total_tokens=total_tokens,
            token_sections=token_sections,
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

    @staticmethod
    def _coerce_token_int(value: Any) -> int:
        try:
            return max(0, int(value or 0))
        except Exception:
            return 0

    def _extract_token_usage(self, llm_response: Any) -> Dict[str, int]:
        """Extract prompt/completion/total token usage from a LangChain response."""
        prompt_tokens = 0
        completion_tokens = 0
        total_tokens = 0

        candidate_dicts: List[Dict[str, Any]] = []
        usage_metadata = getattr(llm_response, "usage_metadata", None)
        if isinstance(usage_metadata, dict):
            candidate_dicts.append(usage_metadata)

        response_metadata = getattr(llm_response, "response_metadata", None)
        if isinstance(response_metadata, dict):
            for key in ("token_usage", "usage"):
                value = response_metadata.get(key)
                if isinstance(value, dict):
                    candidate_dicts.append(value)
            candidate_dicts.append(response_metadata)

        if isinstance(llm_response, dict):
            candidate_dicts.append(llm_response)

        for payload in candidate_dicts:
            prompt_tokens = max(
                prompt_tokens,
                self._coerce_token_int(payload.get("input_tokens")),
                self._coerce_token_int(payload.get("prompt_tokens")),
            )
            completion_tokens = max(
                completion_tokens,
                self._coerce_token_int(payload.get("output_tokens")),
                self._coerce_token_int(payload.get("completion_tokens")),
            )
            total_tokens = max(total_tokens, self._coerce_token_int(payload.get("total_tokens")))

        if total_tokens == 0:
            total_tokens = prompt_tokens + completion_tokens

        return {
            "prompt_tokens": prompt_tokens,
            "completion_tokens": completion_tokens,
            "total_tokens": total_tokens,
        }

    def _new_token_section(self, name: str, label: str) -> Dict[str, Any]:
        return {
            "name": name,
            "label": label,
            "prompt_tokens": 0,
            "completion_tokens": 0,
            "total_tokens": 0,
            "llm_calls": 0,
        }

    def _accumulate_token_usage(
        self,
        section: Dict[str, Any],
        usage_source: Any,
        *,
        llm_calls: int = 0,
    ) -> None:
        usage = self._extract_token_usage(usage_source) if usage_source is not None else {
            "prompt_tokens": 0,
            "completion_tokens": 0,
            "total_tokens": 0,
        }
        section["prompt_tokens"] = int(section.get("prompt_tokens", 0) or 0) + usage["prompt_tokens"]
        section["completion_tokens"] = int(section.get("completion_tokens", 0) or 0) + usage["completion_tokens"]
        section["total_tokens"] = int(section.get("total_tokens", 0) or 0) + usage["total_tokens"]
        section["llm_calls"] = int(section.get("llm_calls", 0) or 0) + int(llm_calls or 0)

    def _format_tool_results(self, results: Dict[str, Any]) -> str:
        """Format tool results as readable text for the synthesis prompt."""
        if not results:
            return "(No tools were called — answering from LLM chemistry knowledge)"

        lines = []
        for name, result in results.items():
            if name == _TOOL_CALL_RESULTS_META_KEY:
                continue
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
        """
        Extract machine-readable outputs from tool results.

        Primary path uses ToolPlugin.structured_projection so tools declare their
        own contribution to `ChemResponse.structured`. A small legacy fallback is
        retained for unannotated tools to avoid breaking older plugins.
        """
        structured: Dict[str, Any] = {}

        plugins = getattr(self.registry, "_plugins", {})
        for tool_name, result in results.items():
            if tool_name == _TOOL_CALL_RESULTS_META_KEY:
                continue
            if not isinstance(result, dict) or not result.get("success"):
                continue
            plugin = plugins.get(tool_name)
            projector = getattr(plugin, "structured_projection", None) if plugin is not None else None
            if callable(projector):
                try:
                    fragment = projector(result) or {}
                    if isinstance(fragment, dict):
                        structured.update(fragment)
                except Exception as exc:
                    if self.verbose:
                        logger.warning(
                            "[ChemCoworker] structured_projection failed for %s: %s",
                            tool_name,
                            exc,
                        )
                continue

            # Legacy fallback for older/unannotated tools (kept intentionally small)
            if tool_name == "detect_reaction_type":
                structured["reaction_type"] = result.get("reaction_type_id") or result.get("reaction_type")
                structured["reaction_family"] = result.get("family_label")
                if result.get("reaction_type_metadata"):
                    structured["reaction_type_metadata"] = result.get("reaction_type_metadata")
            elif tool_name == "analyze_bond_changes":
                structured["bonds_formed"] = result.get("bonds_formed", [])
                structured["bonds_broken"] = result.get("bonds_broken", [])
                structured["key_bond_type"] = result.get("key_bond_type")
            elif tool_name == "recommend_conditions":
                structured["conditions"] = result.get("recommendations", [])
            elif tool_name == "get_molecular_descriptors":
                structured["descriptors"] = result.get("descriptors", {})
                structured["is_drug_like"] = result.get("is_drug_like")
            elif tool_name == "search_reaction_types":
                structured["taxonomy_matches"] = result.get("matches", [])

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
            if tool_name == _TOOL_CALL_RESULTS_META_KEY:
                continue
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

    def _infer_prerequisite_args_from_dependent_call(
        self,
        prereq_name: str,
        dependent_args: Dict[str, Any],
        callables: Dict[str, Callable[..., Any]],
    ) -> Optional[Dict[str, Any]]:
        """
        Best-effort inference of prerequisite args from a dependent tool call.

        Keeps runtime contract enforcement chemistry-safe while avoiding obvious
        dead-ends like `recommend_conditions` without `detect_reaction_type`.
        Returns None if required args cannot be inferred.
        """
        fn = callables.get(prereq_name)
        if fn is None:
            return None

        try:
            sig = inspect.signature(fn)
        except Exception:
            return None

        aliases: Dict[str, List[str]] = {
            "smiles": ["smiles", "reaction_smiles", "target_smiles"],
            "reaction_smiles": ["reaction_smiles", "smiles"],
            "target_smiles": ["target_smiles", "smiles"],
        }
        inferred: Dict[str, Any] = {}

        for param_name, param in sig.parameters.items():
            if param.kind in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD):
                continue

            candidates = aliases.get(param_name, [param_name])
            found = False
            for key in candidates:
                if key in dependent_args:
                    inferred[param_name] = dependent_args[key]
                    found = True
                    break

            if found:
                continue

            if param.default is inspect._empty:
                return None

        return inferred

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
        warning_seen: set[tuple[Any, ...]] = set()
        auto_inserted_this_round: set[str] = set()
        executable_keys: set[tuple[str, str]] = set()

        def _append_executable_once(name: str, args: Dict[str, Any], call_id: str = "") -> Optional[ToolCall]:
            try:
                args_key = json.dumps(args, sort_keys=True, default=str)
            except Exception:
                args_key = str(sorted(args.items()))
            dedupe_key = (name, args_key)
            if dedupe_key in executable_keys:
                return None
            executable_keys.add(dedupe_key)
            tc_obj = ToolCall(name=name, args=dict(args), call_id=call_id or "")
            executable.append(tc_obj)
            return tc_obj

        def _warn_once(msg: str, key: tuple[Any, ...]) -> None:
            if key in warning_seen:
                return
            warning_seen.add(key)
            warnings.append(msg)

        for tc in response_tool_calls:
            tool_name = tc.get("name")
            if tool_name not in callables:
                continue  # preserve existing behavior; caller will produce "not executed" ToolMessage

            plugin = plugins.get(tool_name)
            if plugin is None:
                _append_executable_once(
                    tool_name,
                    dict(tc.get("args", {})),
                    call_id=str(tc.get("id", "") or ""),
                )
                continue

            tc_args = dict(tc.get("args", {}))
            tc_call_id = str(tc.get("id", "") or "")
            missing_prereqs = [p for p in (plugin.prerequisites or []) if p not in completed_tools]
            missing_requires = [k for k in (plugin.requires or []) if k not in provided_keys]
            if not missing_prereqs and not missing_requires:
                _append_executable_once(tool_name, tc_args, call_id=tc_call_id)
                continue

            auto_inserted_prereqs: List[str] = []
            for p in list(missing_prereqs):
                if p in completed_tools or p in requested_names or p in auto_inserted_this_round:
                    continue
                prereq_plugin = plugins.get(p)
                if prereq_plugin is None or p not in callables:
                    continue
                inferred_args = self._infer_prerequisite_args_from_dependent_call(p, tc_args, callables)
                if inferred_args is None:
                    continue
                _append_executable_once(p, inferred_args, call_id=f"auto_{p}_for_{tool_name}")
                auto_inserted_this_round.add(p)
                auto_inserted_prereqs.append(p)

            waiting_on_tools: set[str] = set()
            for p in missing_prereqs:
                if p in requested_names or p in auto_inserted_this_round:
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
                if auto_inserted_prereqs:
                    reason += f" Auto-inserted prerequisite(s): {', '.join(sorted(auto_inserted_prereqs))}."
            else:
                bits: List[str] = []
                if missing_prereqs:
                    bits.append(f"missing prerequisites={missing_prereqs}")
                if missing_requires:
                    bits.append(f"missing required context keys={missing_requires}")
                reason = f"Blocked '{tool_name}' due to unmet tool contracts ({'; '.join(bits)})."

            warning_key = (
                tool_name,
                bool(deferred),
                tuple(sorted(waiting_on_tools)),
                tuple(sorted(missing_prereqs)),
                tuple(sorted(missing_requires)),
            )
            _warn_once(reason, warning_key)
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

    def _build_native_tools(
        self,
        workflow: Optional[WorkflowDefinition] = None,
        active_skill_records: Optional[List[SkillRecord]] = None,
    ) -> List[Any]:
        """Build LangChain StructuredTool objects for the workflow's LLM-visible subset."""
        include_names = self._resolve_native_tool_names(
            workflow=workflow,
            active_skill_records=active_skill_records,
        )
        tools = []
        for name in self.registry.filtered_names(
            llm_exposed_only=False,
            include_names=include_names,
        ):
            plugin = self.registry._plugins[name]
            try:
                tools.append(plugin.to_langchain_tool())
            except Exception as exc:
                if self.verbose:
                    logger.warning(f"[ChemCoworker] Skipping tool '{name}' (schema build failed): {exc}")
        return tools

    def _resolve_native_tool_names(
        self,
        workflow: Optional[WorkflowDefinition] = None,
        active_skill_records: Optional[List[SkillRecord]] = None,
    ) -> List[str]:
        """Resolve the LLM-visible tool surface from workflow policy and active skills."""
        ordered: List[str] = []
        seen: set[str] = set()

        def _append_names(names: List[str]) -> None:
            for name in names:
                if name in seen:
                    continue
                seen.add(name)
                ordered.append(name)

        if workflow and getattr(workflow, "tool_policy", None):
            _append_names(
                self.registry.filtered_names_for_policy(
                    workflow.tool_policy,
                    llm_exposed_only=True,
                )
            )
        if workflow and workflow.llm_visible_tools:
            _append_names(
                self.registry.filtered_names(
                    llm_exposed_only=False,
                    include_names=list(workflow.llm_visible_tools),
                )
            )
        if active_skill_records:
            for record in active_skill_records:
                tool_policy = getattr(record.manifest, "tool_policy", None)
                if tool_policy:
                    _append_names(
                        self.registry.filtered_names_for_policy(
                            tool_policy,
                            llm_exposed_only=False,
                        )
                    )
                _append_names(
                    self.registry.filtered_names(
                        llm_exposed_only=False,
                        include_names=list(record.manifest.tool_allowlist),
                    )
                )
        if not ordered:
            _append_names(self.registry.filtered_names(llm_exposed_only=True))
        return ordered

    def _resolve_active_skill_records(
        self,
        *,
        query: str,
        task_type: str,
        workflow: WorkflowDefinition,
        smiles_present: bool,
    ) -> List[SkillRecord]:
        """Select workflow defaults plus at most one query-matched specialist skill."""
        selected: List[SkillRecord] = []
        seen: set[str] = set()
        skill_registry = getattr(self, "skill_registry", None)
        if skill_registry is None:
            return selected

        def _append_record(record: Optional[SkillRecord]) -> None:
            if record is None:
                return
            skill_id = str(record.manifest.id)
            if skill_id in seen:
                return
            seen.add(skill_id)
            selected.append(record)

        for skill_id in list(getattr(workflow, "default_skill_ids", None) or []):
            _append_record(skill_registry.get_record(skill_id))
        for manifest in skill_registry.match_for_query(
            query=query,
            task_type=task_type,
            smiles_present=smiles_present,
            workflow_name=workflow.name,
        ):
            record = skill_registry.get_record(manifest.id)
            if record is None or record.manifest.id in seen:
                continue
            _append_record(record)
            break
        return selected

    def _build_skill_system_messages(
        self,
        workflow: WorkflowDefinition,
        active_skill_records: List[SkillRecord],
    ) -> List[Any]:
        from langchain_core.messages import SystemMessage

        messages: List[Any] = []
        active_text = format_skill_instruction_block(active_skill_records)
        if active_text:
            messages.append(SystemMessage(content=active_text))
        return messages

    def _tool_progress_marker(self, tool_results: Dict[str, Any]) -> Tuple[int, int, int]:
        completed_tools, provided_keys = self._collect_available_context_from_tool_results(tool_results)
        per_call = tool_results.get(_TOOL_CALL_RESULTS_META_KEY)
        by_call_count = len(per_call) if isinstance(per_call, dict) else 0
        return (len(completed_tools), len(provided_keys), by_call_count)

    def _tool_call_signature(self, response_tool_calls: List[Dict[str, Any]]) -> str:
        normalized: List[Dict[str, Any]] = []
        for tc in response_tool_calls or []:
            normalized.append(
                {
                    "name": str(tc.get("name", "")),
                    "args": dict(tc.get("args", {})),
                }
            )
        try:
            return json.dumps(normalized, sort_keys=True, default=str)
        except Exception:
            return str(normalized)

    def _run_native_tool_loop(
        self,
        query: str,
        task_type: str,
        smiles_list: List[str],
        workflow: WorkflowDefinition,
        primary_smiles: str,
        chemistry_state: Optional[ChemistryRunState] = None,
        active_skill_records: Optional[List[SkillRecord]] = None,
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
            (tool_results, hypothesis, confidence, warnings, llm_call_count, final_answer, messages, token_usage)
            final_answer is non-empty when the model produced a text response after
            finishing all tool calls — caller can use it directly, skipping SYNTHESIZE.
        """
        from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage
        from .plan import ToolCall

        # Use the workflow's pre-selected system prompt (no template vars, API-cacheable)
        native_system = workflow.system_prompt
        max_iterations = workflow.max_iterations

        if active_skill_records is None:
            active_skill_records = self._resolve_active_skill_records(
                query=query,
                task_type=task_type,
                workflow=workflow,
                smiles_present=bool(smiles_list),
            )
        native_tools = self._build_native_tools(workflow=workflow, active_skill_records=active_skill_records)
        if not native_tools:
            logger.warning("[ChemCoworker] No native tools built — falling back to knowledge-only")
            return {}, "", 0.5, [], 0, "", [], self._new_token_section("reason", "Reasoning")

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
        ]
        messages.extend(self._build_skill_system_messages(workflow, active_skill_records))
        messages.append(HumanMessage(content=user_message))
        tool_results: Dict[str, Any] = {}
        tool_results_by_call_id: Dict[str, Any] = {}
        hypothesis = ""
        confidence = 0.5
        warnings: List[str] = []
        llm_call_count = 0
        final_answer = ""
        token_usage = self._new_token_section("reason", "Reasoning")
        callables = self.registry.get_callables()
        runtime_context = self._new_tool_runtime_context(chemistry_state) if chemistry_state is not None else None
        last_tool_signature = ""
        last_progress_marker = (-1, -1, -1)
        repeated_without_progress = 0

        self.event_bus.emit(ChemEvent.PHASE_START, phase="reason")

        for iteration in range(max_iterations):
            try:
                response = llm_with_tools.invoke(messages)
                llm_call_count += 1
                self._accumulate_token_usage(token_usage, response, llm_calls=1)
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

            tool_signature = self._tool_call_signature(response_tool_calls)
            progress_marker = self._tool_progress_marker(tool_results)
            if tool_signature == last_tool_signature and progress_marker == last_progress_marker:
                repeated_without_progress += 1
            else:
                repeated_without_progress = 0
            if repeated_without_progress >= _MAX_REPEATED_TOOL_SIGNATURES_WITHOUT_PROGRESS:
                warnings.append(
                    "Stopped repeated tool loop after identical tool requests produced no new evidence."
                )
                if self.verbose:
                    logger.info(
                        "[ChemCoworker] Breaking repeated tool loop at iteration %d due to no progress.",
                        iteration,
                    )
                break

            # Execute all tool calls in parallel
            tc_list, blocked_results, contract_warnings = self._partition_native_tool_calls_by_contracts(
                response_tool_calls=response_tool_calls,
                callables=callables,
                tool_results=tool_results,
            )
            if blocked_results:
                for _name, _blocked in blocked_results.items():
                    _prev = tool_results.get(_name)
                    # Preserve previously successful evidence; only replace if there
                    # is no prior result or the prior result was unsuccessful.
                    if (
                        isinstance(_prev, dict)
                        and _prev.get("success", True)
                        and isinstance(_blocked, dict)
                        and not _blocked.get("success", True)
                    ):
                        continue
                    tool_results[_name] = _blocked
                warnings.extend(contract_warnings)
                for tc in response_tool_calls:
                    tool_call_id = str(tc.get("id", "") or "")
                    tool_name = str(tc.get("name", "") or "")
                    if tool_call_id and tool_name in blocked_results:
                        tool_results_by_call_id[tool_call_id] = self._copy_result(blocked_results[tool_name])

            if tc_list:
                import time as _t
                t0 = _t.monotonic()
                group_results, group_results_by_call_id = self.executor._run_parallel(
                    tc_list,
                    callables,
                    runtime_context=runtime_context,
                    return_call_results=True,
                )
                for _name, _new in group_results.items():
                    _prev = tool_results.get(_name)
                    if not isinstance(_prev, dict) or not isinstance(_new, dict):
                        tool_results[_name] = _new
                        continue
                    _prev_ok = bool(_prev.get("success", True))
                    _new_ok = bool(_new.get("success", True))
                    # Prefer preserving a successful prior result over a later failure.
                    if _prev_ok and not _new_ok:
                        continue
                    # Prefer newer successful result over older failure.
                    if (not _prev_ok and _new_ok) or (_prev_ok == _new_ok):
                        tool_results[_name] = _new
                tool_results_by_call_id.update(group_results_by_call_id)
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
                result = tool_results_by_call_id.get(tool_call_id)
                if result is None:
                    result = tool_results.get(tool_name, {"success": False, "error": "not executed"})
                result_str = json.dumps(result, default=str)[:3000]
                messages.append(ToolMessage(
                    content=result_str,
                    tool_call_id=tool_call_id,
                ))

            last_tool_signature = tool_signature
            last_progress_marker = self._tool_progress_marker(tool_results)

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
                self._accumulate_token_usage(token_usage, closing_response, llm_calls=1)
                final_answer = self._get_text(closing_response) or ""
                messages.append(closing_response)
                if not hypothesis:
                    hypothesis = final_answer[:300]
            except Exception as exc:
                logger.error(f"[ChemCoworker] Exhaustion-guard closing call failed: {exc}")
                warnings.append(f"Closing call failed after loop exhaustion: {exc}")

        if tool_results_by_call_id:
            tool_results[_TOOL_CALL_RESULTS_META_KEY] = tool_results_by_call_id

        return tool_results, hypothesis, confidence, warnings, llm_call_count, final_answer, messages, token_usage

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
            if tool_name == _TOOL_CALL_RESULTS_META_KEY:
                continue
            if not isinstance(result, dict):
                continue
            for w in result.get("_warnings", []):
                parts.append(f"⚠ {tool_name}: {w}")
        for w in existing_warnings:
            parts.append(f"• {w}")
        # Deduplicate while preserving order
        return "\n".join(dict.fromkeys(parts))

    def _extract_reaction_smiles_candidates_from_text(
        self,
        text: str,
        max_candidates: int = 24,
    ) -> List[str]:
        """Extract likely reaction SMILES candidates from free-form answer text."""
        if not text:
            return []
        candidates: List[str] = []
        for match in _REACTION_SMILES_IN_TEXT_RE.finditer(str(text)):
            left = str(match.group(1) or "").strip().strip("`")
            right = str(match.group(2) or "").strip().strip("`")
            if not left or not right:
                continue
            candidate = f"{left}>>{right}".strip().strip("`'\"")
            candidate = candidate.rstrip(".,;:")
            lower = candidate.lower()
            if "http://" in lower or "https://" in lower:
                continue
            candidates.append(candidate)
            if len(candidates) >= max_candidates:
                break
        return list(dict.fromkeys(candidates))

    def _redact_invalid_reaction_smiles(self, answer: str, invalid_smiles: List[str]) -> str:
        """Redact invalid reaction SMILES so they are not presented as valid chemistry."""
        redacted = str(answer or "")
        for rxn in sorted({str(x) for x in invalid_smiles if str(x)}, key=len, reverse=True):
            redacted = redacted.replace(f"`{rxn}`", "`[INVALID_REACTION_SMILES]`")
            redacted = redacted.replace(rxn, "[INVALID_REACTION_SMILES]")
        return redacted

    def _iter_tool_payloads(self, tool_results: Dict[str, Any]) -> List[Tuple[str, Dict[str, Any]]]:
        """Flatten top-level and per-call tool payloads for evidence checks."""
        payloads: List[Tuple[str, Dict[str, Any]]] = []
        for tool_name, result in (tool_results or {}).items():
            if tool_name == _TOOL_CALL_RESULTS_META_KEY:
                continue
            if isinstance(result, dict):
                payloads.append((str(tool_name), result))
        per_call = tool_results.get(_TOOL_CALL_RESULTS_META_KEY)
        if isinstance(per_call, dict):
            for call_id, result in per_call.items():
                if isinstance(result, dict):
                    payloads.append((f"{_TOOL_CALL_RESULTS_META_KEY}:{call_id}", result))
        return payloads

    @staticmethod
    def _payload_has_condition_recommendations(payload: Dict[str, Any]) -> bool:
        recs = payload.get("recommendations")
        if not isinstance(recs, list) or not recs:
            return False
        key_fields = (
            "catalyst",
            "ligand",
            "base",
            "solvent",
            "secondary_solvent",
            "additive",
            "coupling_reagent",
            "temperature",
            "atmosphere",
        )
        for rec in recs:
            if not isinstance(rec, dict):
                continue
            if any(str(rec.get(k) or "").strip() for k in key_fields):
                return True
        return bool(recs)

    def _has_condition_recommendation_evidence(self, tool_results: Dict[str, Any]) -> bool:
        """Return True if any successful tool payload provides concrete condition recommendations."""
        direct_tools = (
            "recommend_conditions",
            "recommend_reaction_conditions",
            "recommend_forward_conditions",
        )
        for name in direct_tools:
            payload = tool_results.get(name)
            if isinstance(payload, dict) and payload.get("success") and self._payload_has_condition_recommendations(payload):
                return True

        analyze = tool_results.get("analyze_reaction")
        nested = (analyze or {}).get("conditions") if isinstance(analyze, dict) else None
        if isinstance(nested, dict) and nested.get("success") and self._payload_has_condition_recommendations(nested):
            return True

        for _, payload in self._iter_tool_payloads(tool_results):
            if not payload.get("success", False):
                continue
            if self._payload_has_condition_recommendations(payload):
                return True
            for nested_key in ("conditions", "conditions_for_top_disconnection", "conditions_for_top_product"):
                nested_payload = payload.get(nested_key)
                if isinstance(nested_payload, dict) and nested_payload.get("success") and self._payload_has_condition_recommendations(nested_payload):
                    return True
        return False

    def _has_taxonomy_alignment_evidence(self, tool_results: Dict[str, Any]) -> bool:
        """Return True when successful tool outputs include a concrete taxonomy reaction type."""
        def _is_known_reaction_type(value: Any) -> bool:
            text = str(value or "").strip()
            return bool(text) and text.lower() not in {"unknown", "none", "null"}

        for _, payload in self._iter_tool_payloads(tool_results):
            if not payload.get("success", False):
                continue
            if _is_known_reaction_type(payload.get("reaction_type_id")):
                return True
            if _is_known_reaction_type(payload.get("reaction_type")):
                return True
            if _is_known_reaction_type(payload.get("detected_reaction_type")):
                return True
            metadata = payload.get("reaction_type_metadata")
            if isinstance(metadata, dict) and _is_known_reaction_type(metadata.get("id")):
                return True
        return False

    def _answer_claims_condition_recommendations(self, answer: str) -> bool:
        """Detect whether the final answer makes concrete condition recommendations."""
        text = str(answer or "").lower()
        if not text:
            return False
        if ("condition" not in text) and not any(k in text for k in ("catalyst", "ligand", "base", "solvent", "temperature")):
            return False

        negative_markers = (
            "cannot recommend conditions",
            "unable to recommend conditions",
            "no condition recommendation",
            "no condition data",
            "conditions not available",
            "did not return conditions",
            "no reliable condition evidence",
        )
        if any(marker in text for marker in negative_markers):
            return False

        direct_claim_markers = (
            "recommended conditions",
            "condition recommendations",
            "recommend conditions",
            "condition recommendation",
            "suggested conditions",
            "optimal conditions",
        )
        if any(marker in text for marker in direct_claim_markers):
            return True

        slot_hits = sum(1 for k in ("catalyst", "ligand", "base", "solvent", "temperature") if k in text)
        return ("condition" in text) and (slot_hits >= 2)

    def _has_explicit_route_evidence(self, tool_results: Dict[str, Any]) -> bool:
        """Return True when tools produced explicit precursor/product structures for a route."""
        def _has_precursor_pair(row: Dict[str, Any]) -> bool:
            p1 = str(row.get("precursor_1") or "").strip()
            p2 = str(row.get("precursor_2") or "").strip()
            return bool(p1 and p2)

        for _, payload in self._iter_tool_payloads(tool_results):
            if not payload.get("success", False):
                continue

            route = payload.get("route")
            if isinstance(route, list):
                for step in route:
                    if isinstance(step, dict) and _has_precursor_pair(step):
                        return True

            top_disconnection = payload.get("top_disconnection")
            if isinstance(top_disconnection, dict) and _has_precursor_pair(top_disconnection):
                return True

            disconnections = payload.get("disconnections")
            if isinstance(disconnections, dict):
                disconnections = disconnections.get("disconnections")
            if isinstance(disconnections, list):
                for row in disconnections:
                    if isinstance(row, dict) and _has_precursor_pair(row):
                        return True

            products = payload.get("products")
            if isinstance(products, dict):
                products = products.get("products")
            if isinstance(products, list):
                for row in products:
                    if not isinstance(row, dict):
                        continue
                    if str(row.get("product_smiles") or "").strip():
                        return True

            top_product = payload.get("top_product")
            if isinstance(top_product, dict) and str(top_product.get("product_smiles") or "").strip():
                return True

        return False

    def _apply_output_verification_gate(
        self,
        *,
        answer: str,
        tool_results: Dict[str, Any],
        task_type: str,
        active_skill_records: Optional[List[SkillRecord]] = None,
    ) -> Tuple[str, List[str], float, Dict[str, Any]]:
        """Run post-answer validation gate and append explicit verification notes."""
        gate_warnings: List[str] = []
        gate_lines: List[str] = []
        confidence_penalty = 0.0
        safe_answer = str(answer or "")
        active_skill_records = list(active_skill_records or [])
        condition_skill_active = any(
            str(record.manifest.id) == "condition_recommendation"
            for record in active_skill_records
        )

        reaction_candidates = self._extract_reaction_smiles_candidates_from_text(safe_answer)
        invalid_reactions: List[Tuple[str, str]] = []
        if reaction_candidates:
            from .tools._helpers import _validate_reaction_smiles
            from .tools.composite import _evaluate_synthesis_proposal

            for rxn in reaction_candidates:
                cleaned, basic_err = _validate_reaction_smiles(rxn, require_product=True)
                if basic_err:
                    invalid_reactions.append((rxn, basic_err))
                    continue

                try:
                    eval_result = _evaluate_synthesis_proposal(
                        mode="reaction",
                        reaction_smiles=cleaned,
                        include_consistency_checks=False,
                    )
                except Exception as exc:
                    invalid_reactions.append((rxn, f"evaluation failed: {exc}"))
                    continue

                if not isinstance(eval_result, dict) or not eval_result.get("success", False):
                    err = ""
                    if isinstance(eval_result, dict):
                        err = str(eval_result.get("error") or "").strip()
                    invalid_reactions.append((rxn, err or "evaluation returned unsuccessful result"))
                    continue

                verdict = str(eval_result.get("verdict") or "").upper()
                if verdict == "FAIL":
                    critical = list(eval_result.get("critical_failures", []) or [])
                    warns = list(eval_result.get("warnings", []) or [])
                    reason = str(critical[0] if critical else (warns[0] if warns else "evaluator verdict=FAIL"))
                    invalid_reactions.append((rxn, reason))

        if invalid_reactions:
            invalid_smiles = [rxn for rxn, _ in invalid_reactions]
            safe_answer = self._redact_invalid_reaction_smiles(safe_answer, invalid_smiles)
            msg = (
                "Output verification gate removed invalid reaction SMILES from the final answer. "
                f"Count={len(invalid_reactions)}."
            )
            gate_warnings.append(msg)
            confidence_penalty += min(0.30, 0.06 * float(len(invalid_reactions)))
            gate_lines.append(f"- Removed `{len(invalid_reactions)}` invalid reaction SMILES (redacted above).")
            for idx, (rxn, reason) in enumerate(invalid_reactions[:3], 1):
                short_rxn = rxn if len(rxn) <= 80 else f"{rxn[:77]}..."
                gate_lines.append(f"- Invalid reaction {idx}: `{short_rxn}` ({reason}).")
            if len(invalid_reactions) > 3:
                gate_lines.append(f"- Additional invalid reactions omitted: {len(invalid_reactions) - 3}.")

        if self._answer_claims_condition_recommendations(safe_answer):
            if not self._has_condition_recommendation_evidence(tool_results):
                if condition_skill_active:
                    gate_warnings.append(
                        "Output verification: active condition skill requires successful condition-tool evidence before recommending conditions."
                    )
                    confidence_penalty += 0.12
                else:
                    gate_warnings.append(
                        "Output verification: final answer contains condition recommendations without successful condition-tool evidence."
                    )
                    confidence_penalty += 0.08
                gate_lines.append(
                    "- Condition recommendations in this answer are not backed by successful condition-tool outputs."
                )
            if condition_skill_active and not self._has_taxonomy_alignment_evidence(tool_results):
                gate_warnings.append(
                    "Output verification: active condition skill requires taxonomy-aligned reaction identity before recommending conditions."
                )
                confidence_penalty += 0.08
                gate_lines.append(
                    "- Condition recommendations are missing taxonomy-backed reaction identity evidence."
                )

        task_type_norm = str(task_type or "").strip().lower()
        if (
            reaction_candidates
            and len(reaction_candidates) >= 2
            and (("retro" in task_type_norm) or ("route" in task_type_norm))
            and (not self._has_explicit_route_evidence(tool_results))
        ):
            gate_warnings.append(
                "Output verification: answer presents multiple explicit reaction steps, "
                "but tools did not return explicit precursor/product route evidence."
            )
            confidence_penalty += 0.10
            gate_lines.append(
                "- Route details include unverified explicit reaction steps; use tool outputs as authoritative."
            )

        if gate_lines:
            safe_answer += "\n\n---\n⚠ **Output verification gate**\n" + "\n".join(gate_lines)

        report = {
            "invalid_reactions": [
                {"reaction_smiles": rxn, "reason": reason}
                for rxn, reason in invalid_reactions
            ],
            "condition_claim_without_evidence": any(
                ("condition-tool evidence" in str(w).lower()) and ("condition" in str(w).lower())
                for w in gate_warnings
            ),
            "condition_skill_missing_taxonomy_alignment": any(
                "taxonomy-aligned reaction identity" in str(w).lower()
                for w in gate_warnings
            ),
            "unverified_route_steps": any(
                "multiple explicit reaction steps" in str(w).lower()
                for w in gate_warnings
            ),
        }
        return safe_answer, gate_warnings, confidence_penalty, report

    def _build_repair_tool_evidence_summary(
        self,
        tool_results: Dict[str, Any],
        max_chars: int = 3000,
    ) -> str:
        """Build a compact tool-evidence summary for repair prompting."""
        if not isinstance(tool_results, dict) or not tool_results:
            return "(no tool results)"
        parts: List[str] = []
        used = 0
        for tool_name, result in tool_results.items():
            if tool_name == _TOOL_CALL_RESULTS_META_KEY:
                continue
            if not isinstance(result, dict):
                continue
            display = {
                k: v
                for k, v in result.items()
                if k not in ("_warnings",) and v is not None
            }
            try:
                snippet = json.dumps(display, default=str)[:700]
            except Exception:
                snippet = str(display)[:700]
            block = f"[{tool_name}]\n{snippet}"
            if used + len(block) > max_chars:
                break
            parts.append(block)
            used += len(block)
        return "\n\n".join(parts) if parts else "(no usable tool results)"

    def _run_output_repair_pass(
        self,
        *,
        query: str,
        raw_answer: str,
        invalid_reactions: List[Dict[str, str]],
        tool_results: Dict[str, Any],
    ) -> Tuple[str, int, Dict[str, Any], str]:
        """Run a single LLM repair pass and return revised answer + token usage."""
        from langchain_core.messages import HumanMessage, SystemMessage

        token_usage = self._new_token_section("repair", "Repair")
        if not invalid_reactions:
            return raw_answer, 0, token_usage, ""

        invalid_lines = []
        for idx, row in enumerate(invalid_reactions[:6], 1):
            rxn = str(row.get("reaction_smiles") or "")
            reason = str(row.get("reason") or "")
            invalid_lines.append(f"{idx}. {rxn}  | reason: {reason}")
        invalid_block = "\n".join(invalid_lines) if invalid_lines else "(none)"
        tool_evidence = self._build_repair_tool_evidence_summary(tool_results)

        prompt = (
            "The previous chemistry answer contains invalid reaction SMILES and must be repaired.\n\n"
            "Constraints:\n"
            "1) Use ONLY the tool evidence below.\n"
            "2) If a reaction SMILES cannot be made valid from evidence, remove it.\n"
            "3) Do NOT invent precursors, products, or conditions.\n"
            "4) Keep the answer concise and preserve the original intent.\n\n"
            f"Original user query:\n{query}\n\n"
            f"Original answer:\n{raw_answer[:4500]}\n\n"
            f"Invalid reaction SMILES detected:\n{invalid_block}\n\n"
            f"Tool evidence:\n{tool_evidence}\n\n"
            "Return only the revised final answer text."
        )

        try:
            response = self.llm.invoke(
                [
                    SystemMessage(
                        content=(
                            "You are a chemistry assistant repairing an answer after strict structure validation failures."
                        )
                    ),
                    HumanMessage(content=prompt),
                ]
            )
            self._accumulate_token_usage(token_usage, response, llm_calls=1)
            revised = str(self._get_text(response) or "").strip()
            if not revised:
                return raw_answer, 1, token_usage, "Auto-repair returned an empty answer."
            return revised, 1, token_usage, ""
        except Exception as exc:
            return raw_answer, 0, token_usage, f"Auto-repair pass failed: {exc}"

    def _attempt_auto_repair_after_verification(
        self,
        *,
        query: str,
        raw_answer: str,
        verified_answer: str,
        verification_warnings: List[str],
        verification_penalty: float,
        verification_report: Dict[str, Any],
        tool_results: Dict[str, Any],
        task_type: str,
        active_skill_records: Optional[List[SkillRecord]] = None,
    ) -> Tuple[str, List[str], float, List[str], int, Dict[str, Any]]:
        """Try one bounded repair pass after eval failure; keep only improved output."""
        repair_warnings: List[str] = []
        token_usage = self._new_token_section("repair", "Repair")
        invalid_reactions = list((verification_report or {}).get("invalid_reactions") or [])
        if not invalid_reactions:
            return (
                verified_answer,
                list(verification_warnings),
                float(verification_penalty),
                repair_warnings,
                0,
                token_usage,
            )

        revised_answer, repair_calls, repair_token_usage, repair_error = self._run_output_repair_pass(
            query=query,
            raw_answer=raw_answer,
            invalid_reactions=invalid_reactions,
            tool_results=tool_results,
        )
        token_usage = repair_token_usage
        if repair_error:
            repair_warnings.append(repair_error)
        if repair_calls <= 0:
            return (
                verified_answer,
                list(verification_warnings),
                float(verification_penalty),
                repair_warnings,
                repair_calls,
                token_usage,
            )

        (
            repaired_checked_answer,
            repaired_warnings,
            repaired_penalty,
            repaired_report,
        ) = self._apply_output_verification_gate(
            answer=revised_answer,
            tool_results=tool_results,
            task_type=task_type,
            active_skill_records=active_skill_records,
        )

        old_invalid = len(invalid_reactions)
        new_invalid = len(list((repaired_report or {}).get("invalid_reactions") or []))
        improved = (
            (new_invalid < old_invalid and repaired_penalty <= float(verification_penalty))
            or (repaired_penalty + 1e-9 < float(verification_penalty))
        )
        if improved:
            repair_warnings.append(
                f"Auto-repair applied: validation improved (invalid reactions {old_invalid} -> {new_invalid})."
            )
            return (
                repaired_checked_answer,
                repaired_warnings,
                float(repaired_penalty),
                repair_warnings,
                repair_calls,
                token_usage,
            )

        repair_warnings.append(
            "Auto-repair attempted but did not improve validation; kept conservative verified output."
        )
        return (
            verified_answer,
            list(verification_warnings),
            float(verification_penalty),
            repair_warnings,
            repair_calls,
            token_usage,
        )

    def _build_performance_summary(self, tool_results: Dict[str, Any], task_type: str = "") -> str:
        """
        Build a concise HTE-heavy tool timing summary for the final answer.

        Prefers per-call tool results when available so duplicate same-tool calls
        (e.g., multiple route steps) are aggregated correctly.
        """
        raw_calls = tool_results.get(_TOOL_CALL_RESULTS_META_KEY)
        if isinstance(raw_calls, dict) and raw_calls:
            call_results = [
                r for r in raw_calls.values()
                if isinstance(r, dict)
            ]
        else:
            call_results = [
                r for name, r in tool_results.items()
                if name != _TOOL_CALL_RESULTS_META_KEY and isinstance(r, dict)
            ]

        rec_calls = [r for r in call_results if r.get("success") and ("hte_timing_ms" in r or "hte_processing_time_ms" in r)]
        precedent_calls = [r for r in call_results if r.get("success") and "hte_search_timing_ms" in r]
        if not rec_calls and not precedent_calls:
            return ""

        def _f(x: Any) -> float:
            try:
                return float(x)
            except Exception:
                return 0.0

        def _sum_timing(items: List[Dict[str, Any]], field: str, key: str) -> float:
            return round(sum(_f((r.get(key) or {}).get(field)) for r in items), 2)

        def _sum_direct(items: List[Dict[str, Any]], field: str) -> float:
            return round(sum(_f(r.get(field)) for r in items), 2)

        def _max_direct(items: List[Dict[str, Any]], field: str) -> float:
            vals = [_f(r.get(field)) for r in items if r.get(field) is not None]
            return round(max(vals), 2) if vals else 0.0

        def _sum_stage_maps(items: List[Dict[str, Any]], key: str) -> Dict[str, float]:
            agg: Dict[str, float] = {}
            for r in items:
                stage_map = r.get(key) or {}
                if not isinstance(stage_map, dict):
                    continue
                for stage_name, stage_val in stage_map.items():
                    try:
                        val = float(stage_val)
                    except Exception:
                        continue
                    agg[stage_name] = round(float(agg.get(stage_name, 0.0)) + val, 2)
            return agg

        is_route_like = ("retro" in (task_type or "").lower()) or ("route" in (task_type or "").lower())
        title = "Route performance summary" if is_route_like else "Performance summary"
        lines = [f"⚙ **{title}** (HTE-heavy tools)"]

        rec_total_ms = 0.0
        if rec_calls:
            rec_total_ms = _sum_direct(rec_calls, "hte_processing_time_ms") or _sum_timing(rec_calls, "total_ms", "hte_timing_ms")
            rec_max_ms = _max_direct(rec_calls, "hte_processing_time_ms")
            if rec_max_ms <= 0:
                rec_max_ms = max((_f((r.get("hte_timing_ms") or {}).get("total_ms")) for r in rec_calls), default=0.0)
            rec_compute_ms = _sum_timing(rec_calls, "recommend_compute_ms", "hte_timing_ms")
            rec_get_ms = _sum_timing(rec_calls, "recommender_get_ms", "hte_timing_ms")
            lines.append(
                f"- `recommend_conditions`: {len(rec_calls)} call(s), sum {rec_total_ms/1000:.1f}s, "
                f"max {rec_max_ms/1000:.1f}s, compute {rec_compute_ms/1000:.1f}s, init/get {rec_get_ms/1000:.1f}s"
            )
            rec_stage_ms = _sum_stage_maps(rec_calls, "hte_recommender_stage_timing_ms")
            if rec_stage_ms:
                rec_stage_ms.pop("total_ms", None)
                top_stages = sorted(rec_stage_ms.items(), key=lambda kv: kv[1], reverse=True)[:3]
                if top_stages:
                    stage_text = ", ".join(
                        f"{name.replace('_ms', '')} {val/1000:.1f}s" for name, val in top_stages
                    )
                    lines.append(f"- `recommend_conditions` internal stages (sum): {stage_text}")

        precedent_total_ms = 0.0
        if precedent_calls:
            precedent_total_ms = _sum_timing(precedent_calls, "total_ms", "hte_search_timing_ms")
            load_ms = _sum_timing(precedent_calls, "load_family_ms", "hte_search_timing_ms")
            drfp_ms = _sum_timing(precedent_calls, "drfp_rerank_ms", "hte_search_timing_ms")
            lines.append(
                f"- `search_hte_precedent`: {len(precedent_calls)} call(s), sum {precedent_total_ms/1000:.1f}s, "
                f"family-load {load_ms/1000:.1f}s, DRFP {drfp_ms/1000:.1f}s"
            )

        total_hte_ms = round(rec_total_ms + precedent_total_ms, 2)
        if total_hte_ms > 0:
            lines.append(f"- Aggregate HTE tool self-time (sum across calls): {total_hte_ms/1000:.1f}s")
            lines.append("- Note: summed per-call tool times; wall time differs with batching/overlap.")

        return "\n".join(lines)

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
        if not (isinstance(rtype, dict) and rtype.get("success")):
            ar = tool_results.get("analyze_reaction")
            nested = (ar or {}).get("reaction_type") if isinstance(ar, dict) else None
            if isinstance(nested, dict) and nested.get("success"):
                rtype = nested
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
        if not (isinstance(bchg, dict) and bchg.get("success")):
            ar = tool_results.get("analyze_reaction")
            nested = (ar or {}).get("bond_changes") if isinstance(ar, dict) else None
            if isinstance(nested, dict) and nested.get("success"):
                bchg = nested
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
        if not (isinstance(cond, dict) and cond.get("success")):
            ar = tool_results.get("analyze_reaction")
            nested = (ar or {}).get("conditions") if isinstance(ar, dict) else None
            if isinstance(nested, dict) and nested.get("success"):
                cond = nested
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
    ) -> tuple:  # (List[Finding], str, int, Dict[str, Any])
        """
        Phase 6 — Run the adversarial critic pass for retrosynthesis workflows.

        Calls CriticAgent.review() with the main loop's outputs and returns
        structured findings plus a one-sentence verdict.

        Returns:
            (findings, verdict, llm_call_count, token_usage)
            findings is a List[Finding]; verdict is a str; llm_call_count is int.
            Returns ([], "", 0, token_usage) on any failure so the main answer is unaffected.
        """
        from .critic import CriticAgent, Severity

        token_usage = self._new_token_section("critic", "Critic")
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
            self._accumulate_token_usage(
                token_usage,
                getattr(critic, "last_token_usage", None),
                llm_calls=1,
            )
            if self.verbose:
                logger.info(
                    f"[ChemCoworker] Critic: {len(findings)} finding(s); verdict={verdict[:80]!r}"
                )
            return findings, verdict, 1, token_usage
        except Exception as exc:
            logger.warning(f"[ChemCoworker] Critic pass failed: {exc}")
            return [], "", 0, token_usage
        finally:
            self.event_bus.emit(ChemEvent.PHASE_END, phase="critic")


    def _run_revision_pass(
        self,
        query: str,
        original_answer: str,
        findings: List[Any],  # List[Finding]
        critic_verdict: str,
    ) -> tuple:  # (revised_answer: str, llm_call_count: int, Dict[str, Any])
        """
        Phase 6b — Follow-up LLM call that revises the main answer to address
        the critic's findings.

        Returns:
            (revised_answer, llm_call_count, token_usage)
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

        token_usage = self._new_token_section("revision", "Revision")
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
            self._accumulate_token_usage(token_usage, response, llm_calls=1)
            content = getattr(response, "content", response)
            revised = content if isinstance(content, str) else str(content)
            if self.verbose:
                logger.info(
                    f"[ChemCoworker] Revision pass complete "
                    f"({len(revised)} chars → was {len(original_answer)} chars)"
                )
            return revised, 1, token_usage
        except Exception as exc:
            logger.warning(f"[ChemCoworker] Revision pass failed: {exc}")
            return original_answer, 0, token_usage
        finally:
            self.event_bus.emit(ChemEvent.PHASE_END, phase="revision")


def create_coworker(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    verbose: bool = False,
) -> ChemCoworker:
    """Factory function — convenience alternative to ChemCoworker(...)."""
    return ChemCoworker(provider=provider, model=model, verbose=verbose)


