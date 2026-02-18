"""
ToolExecutor — parallel tool execution engine for ChemCoworker.

Runs ExecutionPlan groups sequentially; within each group, independent
tool calls execute concurrently via ThreadPoolExecutor. This reduces
wall-clock time significantly when multiple tools can run in parallel.

Key design:
  - Groups are sequential (G0 runs to completion before G1 starts)
  - Tools within a group are independent and run in parallel
  - Results from earlier groups are NOT passed to later groups automatically
    (the LLM-produced plan should include all args explicitly)
  - Failures in one tool don't abort other tools in the same group
"""
from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Callable, Dict, List, Optional, Tuple

from .plan import ExecutionPlan, ToolCall

logger = logging.getLogger(__name__)

# Maximum parallel threads per group
_MAX_WORKERS = 4


class ToolExecutor:
    """
    Executes an ExecutionPlan by running tool groups in sequence,
    with tools within each group running in parallel.
    """

    def __init__(self, max_workers: int = _MAX_WORKERS, verbose: bool = False):
        self.max_workers = max_workers
        self.verbose = verbose

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

        for group_idx, group in enumerate(plan.groups):
            if self.verbose:
                names = [tc.name for tc in group]
                logger.info(f"[Executor] Group {group_idx}: {names}")

            group_start = time.monotonic()
            group_results = self._run_parallel(group, callables)
            results.update(group_results)

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
        """Execute a group of ToolCalls concurrently."""
        if not calls:
            return {}

        # Single call — skip overhead
        if len(calls) == 1:
            return {calls[0].name: self._execute_one(calls[0], callables)}

        results: Dict[str, Any] = {}
        workers = min(len(calls), self.max_workers)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            future_to_name = {
                pool.submit(self._execute_one, call, callables): call.name
                for call in calls
            }
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    results[name] = future.result()
                except Exception as exc:
                    logger.warning(f"[Executor] Tool '{name}' raised: {exc}")
                    results[name] = {"success": False, "error": str(exc)}

        return results

    def _execute_one(
        self,
        call: ToolCall,
        callables: Dict[str, Callable[..., Any]],
    ) -> Any:
        """Execute a single ToolCall, with error isolation."""
        fn = callables.get(call.name)
        if fn is None:
            return {"success": False, "error": f"Tool '{call.name}' not registered"}
        try:
            start = time.monotonic()
            result = fn(**call.args)
            if self.verbose:
                elapsed = time.monotonic() - start
                logger.debug(f"[Executor] {call.name} completed in {elapsed:.3f}s")
            return result
        except Exception as exc:
            logger.warning(f"[Executor] {call.name} failed: {exc}")
            return {"success": False, "error": str(exc), "tool": call.name}
