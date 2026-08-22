"""Concise stderr progress reporting for long-running CLI operations."""

from __future__ import annotations

import sys
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Iterator, TextIO


def _format_elapsed(seconds: float) -> str:
    """Return a compact human-readable elapsed duration."""

    seconds = max(0.0, seconds)
    if seconds < 10:
        return f"{seconds:.1f}s"
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, remaining = divmod(int(seconds), 60)
    return f"{minutes}m {remaining:02d}s"


@dataclass
class ProgressTask:
    """Render one live TTY line or two concise redirected-output lines."""

    message: str
    enabled: bool = True
    stream: TextIO | None = None
    refresh_interval: float = 0.5
    _started_at: float = field(init=False, default=0.0)
    _stop: threading.Event = field(init=False, default_factory=threading.Event)
    _thread: threading.Thread | None = field(init=False, default=None)
    _is_tty: bool = field(init=False, default=False)
    _rendered_width: int = field(init=False, default=0)
    _result_detail: str = field(init=False, default="")

    def __enter__(self) -> "ProgressTask":
        """Start reporting without blocking the operation being observed."""

        if not self.enabled:
            return self
        if self.refresh_interval <= 0:
            raise ValueError("progress refresh interval must be positive")
        self.stream = self.stream or sys.stderr
        self._started_at = time.monotonic()
        try:
            self._is_tty = bool(self.stream.isatty())
        except (AttributeError, OSError):
            self._is_tty = False
        if self._is_tty:
            self._render_live()
            self._thread = threading.Thread(
                target=self._refresh,
                name="chem-coworker-progress",
                daemon=True,
            )
            self._thread.start()
        else:
            self.stream.write(f"{self.message}...\n")
            self.stream.flush()
        return self

    def set_result(self, detail: str) -> None:
        """Add a short measured result to the completion line."""

        self._result_detail = " ".join(detail.split())

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        """Stop reporting and emit one completion, stop, or failure line."""

        del exc_value, traceback
        if not self.enabled:
            return False
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.refresh_interval * 2))
        elapsed = _format_elapsed(time.monotonic() - self._started_at)
        if exc_type is None:
            status = f"done in {elapsed}"
            if self._result_detail:
                status = f"{status}; {self._result_detail}"
        elif issubclass(exc_type, KeyboardInterrupt):
            status = f"stopped after {elapsed}"
        else:
            status = f"failed after {elapsed}"
        line = f"{self.message} - {status}"
        if self._is_tty:
            self._write_live(line, final=True)
        else:
            assert self.stream is not None
            self.stream.write(f"{line}\n")
            self.stream.flush()
        return False

    def _refresh(self) -> None:
        while not self._stop.wait(self.refresh_interval):
            self._render_live()

    def _render_live(self) -> None:
        elapsed = _format_elapsed(time.monotonic() - self._started_at)
        self._write_live(f"{self.message} - {elapsed}")

    def _write_live(self, line: str, *, final: bool = False) -> None:
        assert self.stream is not None
        padding = " " * max(0, self._rendered_width - len(line))
        self.stream.write(f"\r{line}{padding}{'\n' if final else ''}")
        self.stream.flush()
        self._rendered_width = len(line)


@contextmanager
def concise_progress(
    message: str,
    *,
    enabled: bool = True,
    stream: TextIO | None = None,
) -> Iterator[ProgressTask]:
    """Report live elapsed time while preserving stdout for result data."""

    with ProgressTask(message=message, enabled=enabled, stream=stream) as task:
        yield task


__all__ = ["ProgressTask", "concise_progress"]
