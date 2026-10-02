"""Shared process control for recorded scientific execution and agent runtimes."""

from __future__ import annotations

import os
import signal
import subprocess
from typing import Any


def hidden_process_options() -> dict[str, Any]:
    """Keep Windows child processes hidden; return no extra POSIX flags."""
    return {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}

def stop_process_tree(process: subprocess.Popen[Any]) -> None:
    """Stop only the supplied child process tree, preserving existing platform behavior."""
    if process.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True, timeout=15, **hidden_process_options(),
        )
    else:
        os.killpg(process.pid, signal.SIGKILL)
    if process.poll() is None:
        process.kill()
    process.wait(timeout=15)
