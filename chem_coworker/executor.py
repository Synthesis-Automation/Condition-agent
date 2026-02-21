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

from .plan import ExecutionPlan, ToolCall, _SMILES_PLACEHOLDER_RE

logger = logging.getLogger(__name__)

# Maximum parallel threads per group
_MAX_WORKERS = 4

# Callback types (A3 — streaming, A1 — hooks)
ProgressCallback = Callable[[str, str, float], None]
# signature: (event: "start"|"done"|"error", tool_name: str, elapsed_s: float) -> None

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


# ---------------------------------------------------------------------------
# Inter-group SMILES propagation helpers
# ---------------------------------------------------------------------------

def _extract_resolved_smiles(group_results: Dict[str, Any]) -> Optional[str]:
    """Return the SMILES resolved by a name-resolver tool, or None.

    Called after each group executes. When G0 contains resolve_to_smiles or
    resolve_name_to_smiles (legacy) and it succeeds, its 'smiles' field is
    propagated to patch G1+ tool args that still hold a placeholder string.
    """
    # Check both new consolidated name and old name for backward compat
    r = group_results.get("resolve_to_smiles") or group_results.get("resolve_name_to_smiles")
    if not isinstance(r, dict) or not r.get("success", False):
        return None
    for key in ("smiles", "resolved_smiles", "canonical_smiles"):
        val = r.get(key)
        if val and isinstance(val, str) and val.lower() not in ("(none)", "none", ""):
            return val
    return None


def _substitute_in_remaining_groups(
    plan: ExecutionPlan,
    completed_group_idx: int,
    resolved_smiles: str,
) -> None:
    """Replace SMILES placeholder strings in all groups after completed_group_idx.

    Modifies plan.groups in place. Only replaces arg values that match
    _SMILES_PLACEHOLDER_RE — real SMILES already in args are left untouched.
    """
    for group in plan.groups[completed_group_idx + 1 :]:
        for call in group:
            call.args = {
                k: resolved_smiles if (
                    isinstance(v, str) and _SMILES_PLACEHOLDER_RE.match(v)
                ) else v
                for k, v in call.args.items()
            }


class ToolExecutor:
    """
    Executes an ExecutionPlan by running tool groups in sequence,
    with tools within each group running in parallel.
    """

    def __init__(
        self,
        max_workers: int = _MAX_WORKERS,
        verbose: bool = False,
        progress_cb: Optional[ProgressCallback] = None,
        hooks: Optional[Any] = None,   # HookRegistry — typed as Any to avoid circular import
    ):
        self.max_workers = max_workers
        self.verbose = verbose
        self.progress_cb = progress_cb
        self.hooks = hooks

    def run_plan(
        self,
        plan: ExecutionPlan,
        callables: Dict[str, Callable[..., Any]],
    ) -> Dict[str, Any]:
        """
        Execute an ExecutionPlan and return all tool results.

        Args:
            plan: ExecutionPlan from PlanParser.
            callables: {tool_name: fn} dict from ToolRegistry.get_callables().

        Returns:
            {tool_name: result_payload} for all tools that were executed.
        """
        results: Dict[str, Any] = {}
        total_start = time.monotonic()

        # Pre-warm all tool-side chemtools imports in the main thread so that
        # concurrent threads inside ThreadPoolExecutor never race on _ModuleLock.
        _prewarm_tool_imports()

        for group_idx, group in enumerate(plan.groups):
            if self.verbose:
                names = [tc.name for tc in group]
                logger.info(f"[Executor] Group {group_idx}: {names}")

            group_start = time.monotonic()
            group_results = self._run_parallel(group, callables)
            results.update(group_results)

            # Propagate resolved SMILES: if resolve_name_to_smiles just ran,
            # patch placeholder args in all remaining groups so they receive
            # the real SMILES instead of the literal placeholder string.
            resolved = _extract_resolved_smiles(group_results)
            if resolved:
                _substitute_in_remaining_groups(plan, group_idx, resolved)
                if self.verbose:
                    logger.info(
                        f"[Executor] Propagated resolved SMILES '{resolved}' "
                        f"to groups {group_idx + 1}+"
                    )

            if self.verbose:
                elapsed = time.monotonic() - group_start
                logger.info(
                    f"[Executor] Group {group_idx} done in {elapsed:.2f}s: "
                    f"{list(group_results.keys())}"
                )

        if self.verbose:
            total = time.monotonic() - total_start
            logger.info(
                f"[Executor] All {len(plan.all_tool_calls)} tools done in {total:.2f}s. "
                f"Groups: {len(plan.groups)}"
            )

        return results

    def _run_parallel(
        self,
        calls: List[ToolCall],
        callables: Dict[str, Callable[..., Any]],
    ) -> Dict[str, Any]:
        """Execute a group of ToolCalls concurrently, firing progress events."""
        if not calls:
            return {}

        t0 = time.monotonic()

        # Fire "start" for every tool in the group before submitting (A3)
        if self.progress_cb:
            for call in calls:
                self.progress_cb("start", call.name, 0.0)

        # Single call — skip thread overhead
        if len(calls) == 1:
            result = self._execute_one(calls[0], callables)
            if self.progress_cb:
                elapsed = time.monotonic() - t0
                has_error = isinstance(result, dict) and not result.get("success", True)
                self.progress_cb("error" if has_error else "done", calls[0].name, elapsed)
            return {calls[0].name: result}

        results: Dict[str, Any] = {}
        workers = min(len(calls), self.max_workers)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_name = {
                pool.submit(self._execute_one, call, callables): call.name
                for call in calls
            }
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                elapsed = time.monotonic() - t0
                try:
                    results[name] = future.result()
                    has_error = isinstance(results[name], dict) and not results[name].get("success", True)
                    event = "error" if has_error else "done"
                except Exception as exc:
                    logger.warning(f"[Executor] Tool '{name}' raised: {exc}")
                    results[name] = {"success": False, "error": str(exc)}
                    event = "error"
                if self.progress_cb:
                    self.progress_cb(event, name, elapsed)

        return results

    def _execute_one(
        self,
        call: ToolCall,
        callables: Dict[str, Callable[..., Any]],
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

            start = time.monotonic()
            result = fn(**call.args)
            ctx.result = result

            if self.verbose:
                elapsed = time.monotonic() - start
                logger.debug(f"[Executor] {call.name} completed in {elapsed:.3f}s")

            # Post-hooks — may return a replacement result
            if self.hooks:
                override = self.hooks.fire_post(ctx)
                if override is not None:
                    result = override

            return result

        except HookAbort as e:
            return {"success": False, "error": f"[hook aborted] {e}"}
        except Exception as exc:
            logger.warning(f"[Executor] {call.name} failed: {exc}")
            return {"success": False, "error": str(exc), "tool": call.name}
