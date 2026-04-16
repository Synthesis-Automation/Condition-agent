"""
ToolExecutor — parallel tool execution engine for ChemCoworker.

Runs ExecutionPlan groups sequentially; within each group, independent
tool calls execute concurrently via ThreadPoolExecutor. This reduces
wall-clock time significantly when multiple tools can run in parallel.

Key design:
  - Groups are sequential (G0 runs to completion before G1 starts)
  - Tools within a group are independent and run in parallel
  - resolve_name_to_smiles is special: its resolved SMILES is automatically
    substituted into placeholder args of all subsequent groups
  - Failures in one tool don't abort other tools in the same group
"""
from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Callable, Dict, List, Optional, Tuple

from .plan import ExecutionPlan, ToolCall
from .event_bus import ChemEvent as _ChemEvent

logger = logging.getLogger(__name__)

# Maximum parallel threads per group
_MAX_WORKERS = 4
_HTE_HEAVY_TOOLS = {
    "recommend_conditions",
    "get_literature_condition_evidence",
    "get_motif_condition_evidence",
    "get_rule_condition_evidence",
    "get_similarity_condition_evidence",
    "search_hte_precedent",
}

# ---------------------------------------------------------------------------
# Import pre-warmer
# ---------------------------------------------------------------------------

_prewarmed = False


def _prewarm_tool_imports() -> None:
    """Import all chemtools modules used by tools in the main thread.

    Prevents _ModuleLock deadlock: when ThreadPoolExecutor runs tools
    concurrently, multiple threads race to lazily import the same chemtools
    submodules. Python's per-module import lock causes a deadlock if two
    threads each hold a lock the other needs.

    Calling this once in the main thread before any ThreadPoolExecutor
    ensures all modules are in sys.modules, so thread-side imports are
    instant no-ops with no lock acquisition.
    """
    global _prewarmed
    if _prewarmed:
        return
    _prewarmed = True
    _mods = [
        "chemtools.featurizers.analysis.smiles",  # normalize_reaction
        "chemtools.featurizers.unified",           # detect_reaction_type, featurize_reaction
        "chemtools._atom_mapping",                 # analyze_bond_changes
        "chemtools.util.functional_groups",        # inspect_functional_groups
        "chemtools.recommend.hte_adapter",         # recommend_conditions
        "chemtools.reagent.lookup",                # lookup_reagent, list_reagents_by_role
    ]
    for mod in _mods:
        try:
            __import__(mod)
        except Exception:
            pass  # Missing optional dependency — let the tool handle it


def _extract_provides_smiles(
    group_results: Dict[str, Any],
    plugins: Dict[str, Any],  # {name: ToolPlugin}
) -> Optional[str]:
    """Data-contract SMILES resolution — finds SMILES from executed tools that
    declare provides=["resolved_smiles"] (or "smiles" / "canonical_smiles").
    """
    for tool_name, result in group_results.items():
        if not isinstance(result, dict) or not result.get("success", False):
            continue
        plugin = plugins.get(tool_name)
        if plugin is None:
            continue
        for provided_key in plugin.provides:
            if provided_key in ("resolved_smiles", "smiles", "canonical_smiles"):
                val = result.get(provided_key)
                if val and isinstance(val, str) and val.lower() not in ("(none)", "none", ""):
                    return val
    return None


class ToolExecutor:
    """
    Executes an ExecutionPlan by running tool groups in sequence,
    with tools within each group running in parallel.
    """

    def __init__(
        self,
        max_workers: int = _MAX_WORKERS,
        verbose: bool = False,
        event_bus: Optional[Any] = None,   # EventBus — typed as Any to avoid circular import
        hooks: Optional[Any] = None,   # HookRegistry — typed as Any to avoid circular import
        registry: Optional[Any] = None,  # ToolRegistry — for data-contract resolution (Phase 1)
        runtime_context: Optional[Any] = None,  # Phase 11: per-run tool runtime context
    ):
        self.max_workers = max_workers
        self.verbose = verbose
        self.event_bus = event_bus
        self.hooks = hooks
        self.registry = registry
        self.runtime_context = runtime_context

    def _run_parallel(
        self,
        calls: List[ToolCall],
        callables: Dict[str, Callable[..., Any]],
        runtime_context: Optional[Any] = None,
        return_call_results: bool = False,
    ) -> Any:
        """Execute a group of ToolCalls concurrently, firing progress events."""
        if not calls:
            return ({}, {}) if return_call_results else {}

        # Fire TOOL_START for every tool in the group before submitting
        if self.event_bus:
            for call in calls:
                self.event_bus.emit(_ChemEvent.TOOL_START, tool_name=call.name)

        # Single call — skip thread overhead
        if len(calls) == 1:
            call = calls[0]
            ok, payload, elapsed = self._execute_one_timed(
                call, callables, runtime_context=runtime_context
            )
            if ok:
                result = payload
            else:
                exc = payload
                logger.warning(f"[Executor] Tool '{call.name}' raised: {exc}")
                result = {"success": False, "error": str(exc)}
            if self.event_bus:
                has_error = isinstance(result, dict) and not result.get("success", True)
                if has_error:
                    self.event_bus.emit(
                        _ChemEvent.TOOL_ERROR,
                        tool_name=calls[0].name,
                        elapsed_s=elapsed,
                        error=str(result.get("error", "")),
                    )
                else:
                    self.event_bus.emit(
                        _ChemEvent.TOOL_DONE,
                        tool_name=calls[0].name,
                        elapsed_s=elapsed,
                    )
            by_name = {call.name: result}
            by_call_id = {(call.call_id or call.name): result}
            return (by_name, by_call_id) if return_call_results else by_name

        results: Dict[str, Any] = {}
        results_by_call_id: Dict[str, Any] = {}
        workers = self._compute_parallel_workers(calls)
        future_to_call: Dict[Any, ToolCall] = {}

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for call in calls:
                fut = pool.submit(self._execute_one_timed, call, callables, runtime_context)
                future_to_call[fut] = call
            for future in as_completed(future_to_call):
                call = future_to_call[future]
                name = call.name
                try:
                    ok, payload, elapsed = future.result()
                    if ok:
                        result_obj = payload
                    else:
                        raise payload
                    results[name] = result_obj
                    results_by_call_id[call.call_id or name] = result_obj
                    has_error = isinstance(result_obj, dict) and not result_obj.get("success", True)
                    if self.event_bus:
                        if has_error:
                            self.event_bus.emit(
                                _ChemEvent.TOOL_ERROR,
                                tool_name=name,
                                elapsed_s=elapsed,
                                error=str(result_obj.get("error", "")),
                            )
                        else:
                            self.event_bus.emit(
                                _ChemEvent.TOOL_DONE, tool_name=name, elapsed_s=elapsed
                            )
                except Exception as exc:
                    logger.warning(f"[Executor] Tool '{name}' raised: {exc}")
                    err = {"success": False, "error": str(exc)}
                    results[name] = err
                    results_by_call_id[call.call_id or name] = err
                    if self.event_bus:
                        self.event_bus.emit(
                            _ChemEvent.TOOL_ERROR,
                            tool_name=name,
                            elapsed_s=elapsed,
                            error=str(exc),
                        )

        return (results, results_by_call_id) if return_call_results else results

    def _execute_one_timed(
        self,
        call: ToolCall,
        callables: Dict[str, Callable[..., Any]],
        runtime_context: Optional[Any] = None,
    ) -> Tuple[bool, Any, float]:
        """Execute one tool and return (ok, result_or_exc, elapsed_s)."""
        start = time.monotonic()
        try:
            result = self._execute_one(call, callables, runtime_context=runtime_context)
            return True, result, (time.monotonic() - start)
        except Exception as exc:  # pragma: no cover - defensive; _execute_one usually isolates
            return False, exc, (time.monotonic() - start)

    def _compute_parallel_workers(self, calls: List[ToolCall]) -> int:
        """
        Compute effective worker count for a batch.

        Safeguard: if multiple HTE-heavy tools are present in the same batch,
        serialize the batch to avoid duplicated dataset initialization and heavy
        contention on the HTE recommender path.
        """
        if not calls:
            return 1
        base = max(1, min(len(calls), self.max_workers))
        heavy_count = sum(1 for c in calls if c.name in _HTE_HEAVY_TOOLS)
        if heavy_count >= 2:
            if self.verbose:
                logger.info(
                    "[Executor] Serializing batch with %d HTE-heavy tool calls: %s",
                    heavy_count,
                    [c.name for c in calls],
                )
            return 1
        return base

    def _execute_one(
        self,
        call: ToolCall,
        callables: Dict[str, Callable[..., Any]],
        runtime_context: Optional[Any] = None,
    ) -> Any:
        """Execute a single ToolCall with error isolation and optional hooks (A1)."""
        fn = callables.get(call.name)
        if fn is None:
            return {"success": False, "error": f"Tool '{call.name}' not registered"}

        from .hooks import HookAbort, HookContext
        ctx = HookContext(tool_name=call.name, args=dict(call.args))

        try:
            # Pre-hooks — may raise HookAbort to cancel the call
            if self.hooks:
                self.hooks.fire_pre(ctx)

            active_runtime_context = runtime_context if runtime_context is not None else self.runtime_context
            runtime_token = None
            if active_runtime_context is not None:
                try:
                    from .tool_runtime import set_current_tool_runtime_context
                    runtime_token = set_current_tool_runtime_context(active_runtime_context)
                except Exception:
                    runtime_token = None

            start = time.monotonic()
            try:
                result = fn(**call.args)
            finally:
                if runtime_token is not None:
                    try:
                        from .tool_runtime import reset_current_tool_runtime_context
                        reset_current_tool_runtime_context(runtime_token)
                    except Exception:
                        pass
            ctx.result = result

            if self.verbose:
                elapsed = time.monotonic() - start
                logger.debug(f"[Executor] {call.name} completed in {elapsed:.3f}s")

            # Post-hooks — may return a replacement result
            if self.hooks:
                override = self.hooks.fire_post(ctx)
                if override is not None:
                    result = override

            # Phase 1: run tool validators registered on ToolPlugin
            if self.registry and call.name in self.registry._plugins:
                plugin = self.registry._plugins[call.name]
                for validator_fn in plugin.validators:
                    try:
                        warning = validator_fn(result)
                        if warning and isinstance(result, dict):
                            result.setdefault("_warnings", []).append(warning)
                    except Exception as vexc:
                        logger.debug(f"[Executor] Validator for {call.name} raised: {vexc}")

            return result

        except HookAbort as e:
            return {"success": False, "error": f"[hook aborted] {e}"}
        except Exception as exc:
            logger.warning(f"[Executor] {call.name} failed: {exc}")
            return {"success": False, "error": str(exc), "tool": call.name}
