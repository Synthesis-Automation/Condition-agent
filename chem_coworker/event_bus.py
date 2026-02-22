"""
Phase 2: EventBus — replaces 5 individual callback params with a single observer.

Usage:
    bus = EventBus()

    @bus.on(ChemEvent.TOOL_DONE)
    def handle_done(tool_name, elapsed_s, **_):
        print(f"{tool_name} finished in {elapsed_s:.2f}s")

    # or via subscribe (chainable):
    bus.subscribe(ChemEvent.STREAM_TOKEN, lambda token, **_: print(token, end=""))

    bus.emit(ChemEvent.TOOL_DONE, tool_name="detect_reaction_type", elapsed_s=0.12)
"""

from __future__ import annotations

import logging
from collections import defaultdict
from enum import Enum, auto
from typing import Any, Callable, Dict, List

logger = logging.getLogger(__name__)


class ChemEvent(Enum):
    TOOL_START = auto()    # kwargs: tool_name: str
    TOOL_DONE = auto()     # kwargs: tool_name: str, elapsed_s: float
    TOOL_ERROR = auto()    # kwargs: tool_name: str, elapsed_s: float, error: str
    PHASE_START = auto()   # kwargs: phase: str  e.g. "reason" | "synth" | "observe" | "compact"
    PHASE_END = auto()     # kwargs: phase: str
    PRE_SYNTH = auto()     # kwargs: hypothesis, confidence, rationale, tools_called, plan_revised
    STREAM_TOKEN = auto()  # kwargs: token: str
    COMPACT_START = auto() # kwargs: (none)
    COMPACT_END = auto()   # kwargs: (none)


Handler = Callable[..., None]


class EventBus:
    """
    Lightweight synchronous event bus.

    - Handlers are called in subscription order.
    - Exceptions inside handlers are logged and swallowed so that one broken
      handler never aborts the chemistry pipeline.
    - `bool(bus)` is True only when at least one handler is registered (lets
      callers skip the emit() call entirely as an optimization).
    """

    def __init__(self) -> None:
        self._handlers: Dict[ChemEvent, List[Handler]] = defaultdict(list)

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def subscribe(self, event: ChemEvent, handler: Handler) -> "EventBus":
        """Register *handler* for *event*. Returns self for chaining."""
        self._handlers[event].append(handler)
        return self

    def on(self, event: ChemEvent) -> Callable[[Handler], Handler]:
        """Decorator: @bus.on(ChemEvent.TOOL_DONE)"""
        def decorator(fn: Handler) -> Handler:
            self.subscribe(event, fn)
            return fn
        return decorator

    # ------------------------------------------------------------------
    # Emission
    # ------------------------------------------------------------------

    def emit(self, event: ChemEvent, **kwargs: Any) -> None:
        """
        Call every handler registered for *event* with **kwargs.

        Errors inside handlers are caught, logged, and ignored so a UI bug
        can never crash the chemistry pipeline.
        """
        for handler in self._handlers[event]:
            try:
                handler(**kwargs)
            except Exception as exc:  # noqa: BLE001
                logger.warning(
                    "[EventBus] Handler %s for %s raised: %s",
                    getattr(handler, "__name__", repr(handler)),
                    event.name,
                    exc,
                )

    # ------------------------------------------------------------------
    # Truthiness
    # ------------------------------------------------------------------

    def __bool__(self) -> bool:
        """True when at least one handler has been registered on any event."""
        return any(self._handlers.values())
