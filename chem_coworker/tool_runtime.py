"""
Per-tool-call runtime context accessors.

This module provides a contextvar-backed runtime context object that the
executor sets for each tool invocation. Tool implementations can read the
current context without changing their public function signatures (important
for LangChain tool schema generation).
"""
from __future__ import annotations

from contextvars import ContextVar, Token
from typing import Any, Optional


_CURRENT_TOOL_RUNTIME_CONTEXT: ContextVar[Any] = ContextVar(
    "chem_coworker_current_tool_runtime_context",
    default=None,
)


def set_current_tool_runtime_context(ctx: Any) -> Token:
    """Set current tool runtime context for this execution/thread."""
    return _CURRENT_TOOL_RUNTIME_CONTEXT.set(ctx)


def reset_current_tool_runtime_context(token: Token) -> None:
    """Reset current tool runtime context using the token returned by set()."""
    _CURRENT_TOOL_RUNTIME_CONTEXT.reset(token)


def get_current_tool_runtime_context() -> Optional[Any]:
    """Return the current tool runtime context, if one has been set."""
    return _CURRENT_TOOL_RUNTIME_CONTEXT.get()


__all__ = [
    "get_current_tool_runtime_context",
    "reset_current_tool_runtime_context",
    "set_current_tool_runtime_context",
]
