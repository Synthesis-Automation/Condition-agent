"""
Tool hook system for ChemCoworker — inspired by Claude Code's hook architecture.

Hooks are lightweight callbacks that fire before (pre) or after (post) every
tool execution inside ToolExecutor. They enable:

  Pre-hooks:
    - SMILES validation before a tool receives bad input
    - Rate limiting / circuit breakers
    - Argument sanitization
    - Block a tool entirely (raise HookAbort)

  Post-hooks:
    - Audit logging of every tool call and its outcome
    - Caching successful results for repeated calls
    - Auto-appending successful conditions to notes
    - Modify or enrich the result dict returned to the LLM

Usage:
    from chem_coworker.hooks import HookRegistry, HookContext, HookAbort

    hooks = HookRegistry()

    # Log every tool call
    def audit(ctx: HookContext) -> None:
        print(f"[HOOK] {ctx.tool_name} -> {'ok' if ctx.result.get('success') else 'FAIL'}")
    hooks.add_post(audit)

    # Block a specific tool
    def block_expensive(ctx: HookContext) -> None:
        if ctx.tool_name == "get_molecular_descriptors":
            raise HookAbort("Blocked by policy.")
    hooks.add_pre(block_expensive)

    coworker = ChemCoworker(..., hooks=hooks)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional


class HookAbort(Exception):
    """
    Raise from a pre-hook to cancel a tool call.
    The tool will not run; the result will be {"success": False, "error": "[hook aborted] <msg>"}.
    """


@dataclass
class HookContext:
    """
    Shared context object passed to every hook.

    In pre-hooks:  result is None (tool hasn't run yet).
    In post-hooks: result holds the tool's return value.
    """
    tool_name: str
    args: Dict[str, Any]
    result: Optional[Dict[str, Any]] = None   # None in pre-hooks
    error: Optional[str] = None               # set if tool raised an exception


# Type aliases
PreHook  = Callable[[HookContext], None]
PostHook = Callable[[HookContext], Optional[Dict[str, Any]]]


@dataclass
class HookRegistry:
    """
    Registry of pre- and post-tool hooks.

    Pre-hooks run before the tool function is called.
    Post-hooks run after the tool returns and may optionally return a replacement result.
    Methods return self to enable chaining: registry.add_pre(f).add_post(g)
    """
    pre_hooks:  List[PreHook]  = field(default_factory=list)
    post_hooks: List[PostHook] = field(default_factory=list)

    def add_pre(self, fn: PreHook) -> "HookRegistry":
        """Register a pre-hook. Returns self for chaining."""
        self.pre_hooks.append(fn)
        return self

    def add_post(self, fn: PostHook) -> "HookRegistry":
        """Register a post-hook. Returns self for chaining."""
        self.post_hooks.append(fn)
        return self

    def fire_pre(self, ctx: HookContext) -> None:
        """Run all pre-hooks. Any hook may raise HookAbort to cancel the tool call."""
        for fn in self.pre_hooks:
            fn(ctx)

    def fire_post(self, ctx: HookContext) -> Optional[Dict[str, Any]]:
        """
        Run all post-hooks. If any hook returns a non-None dict, it replaces
        the tool result. The last non-None return wins.
        """
        override: Optional[Dict[str, Any]] = None
        for fn in self.post_hooks:
            result = fn(ctx)
            if result is not None:
                override = result
        return override
