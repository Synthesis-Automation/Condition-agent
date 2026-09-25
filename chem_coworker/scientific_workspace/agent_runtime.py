"""Optional agent harness adapters; deterministic scientific packages stay model-free."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
from threading import Event
from time import monotonic
from typing import Any, Callable, Protocol


from .answer_contracts import ANSWER_SCHEMA


@dataclass(frozen=True)
class AgentResult:
    """An unreviewed answer from a runtime, never a chemistry validation result."""

    answer: dict[str, Any]
    thread_id: str
    usage: dict[str, Any]


class AgentRuntime(Protocol):
    """Replaceable conversation harness boundary, independent of scientific operations."""

    def describe(self) -> dict[str, Any]:
        """Report configuration and local availability, without exposing credentials."""

    def run(
        self, *, prompt: str, workspace: Path, turn_directory: Path,
        thread_id: str | None, cancel: Event,
        on_event: Callable[[dict[str, Any]], None],
    ) -> AgentResult:
        """Run one complete agent turn with iterative tools and optional conversation memory."""


class AgentStopped(RuntimeError):
    """Cancellation or elapsed wall-clock budget; partial evidence remains on disk."""

    def __init__(self, status: str) -> None:
        super().__init__(status)
        self.status = status


def _hidden_process_options() -> dict[str, Any]:
    return {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}


def find_codex(executable: str | None = None) -> tuple[str, str]:
    """Find a modern CLI, including an IDE-bundled binary when PATH is obsolete."""
    selected = executable or os.environ.get("SCIENTIFIC_CODEX_PATH")
    candidates = [selected] if selected else [shutil.which("codex.exe"), shutil.which("codex")]
    if not selected:
        candidates.extend(str(path) for path in sorted(
            (Path.home() / ".vscode" / "extensions").glob("openai.chatgpt-*/bin/*/codex.exe"),
            key=lambda path: path.stat().st_mtime, reverse=True,
        ))
    for candidate in dict.fromkeys(item for item in candidates if item):
        # Never launch a PowerShell/batch shim through a shell. Native binaries or
        # POSIX executable scripts are supported; prompts only travel over stdin.
        if Path(candidate).suffix.lower() in {".cmd", ".bat", ".ps1"}:
            continue
        try:
            help_result = subprocess.run(
                [candidate, "exec", "--help"], capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=10, **_hidden_process_options(),
            )
            if help_result.returncode or not all(
                flag in help_result.stdout for flag in ("--output-schema", "resume", "--json")
            ):
                continue
            version = subprocess.run(
                [candidate, "--version"], capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=10, check=True,
                **_hidden_process_options(),
            ).stdout.strip()
            return str(Path(candidate).resolve()), version
        except (OSError, subprocess.SubprocessError):
            continue
    raise RuntimeError(
        "A current native Codex CLI with exec, resume and --output-schema is required. "
        "Install/update Codex, run codex login, or set SCIENTIFIC_CODEX_PATH to its executable."
    )


def _stop_process(process: subprocess.Popen[Any]) -> None:
    """Stop only the runtime process tree started for this turn."""
    if process.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True, timeout=15, **_hidden_process_options(),
        )
    else:
        os.killpg(process.pid, signal.SIGKILL)
    if process.poll() is None:
        process.kill()
    process.wait(timeout=15)


class CodexRuntime:
    """Use Codex's tool loop and saved thread, rather than implement an LLM planner."""

    def __init__(
        self, *, executable: str | None = None, model: str | None = None,
        timeout_seconds: float = 900,
    ) -> None:
        if not 1 <= timeout_seconds <= 7200:
            raise ValueError("timeout_seconds must be between 1 and 7200")
        self.executable, self.version = find_codex(executable)
        self.model = model
        self.timeout_seconds = timeout_seconds

    def describe(self) -> dict[str, Any]:
        """Local discovery does not assert provider authentication or model availability."""
        return {
            "runtime": "codex_exec", "version": self.version,
            "model": self.model or "Codex configuration default",
            "sandbox": "workspace-write", "timeout_seconds": self.timeout_seconds,
            "authentication": "existing Codex login; checked when a turn runs",
        }

    def command(self, workspace: Path, turn_directory: Path, thread_id: str | None) -> list[str]:
        """Build explicit arguments; never interpolate a user's question into a shell."""
        command = [
            self.executable, "exec", "--sandbox", "workspace-write",
            "-c", 'approval_policy="never"', "--cd", str(workspace),
        ]
        if thread_id:
            command.extend(["resume", thread_id])
        command.extend([
            "--json", "--skip-git-repo-check", "--output-schema",
            str(turn_directory / "answer-schema.json"), "--output-last-message",
            str(turn_directory / "agent-final.json"),
        ])
        if self.model:
            command.extend(["--model", self.model])
        command.append("-")
        return command

    def run(
        self, *, prompt: str, workspace: Path, turn_directory: Path,
        thread_id: str | None, cancel: Event,
        on_event: Callable[[dict[str, Any]], None],
    ) -> AgentResult:
        """Capture JSONL events, enforce a deadline, and preserve complete local logs."""
        (turn_directory / "answer-schema.json").write_text(json.dumps(ANSWER_SCHEMA), "utf-8")
        prompt_path = turn_directory / "prompt.txt"
        prompt_path.write_text(prompt, "utf-8")
        output_path = turn_directory / "runtime.jsonl"
        stderr_path = turn_directory / "runtime.stderr.txt"
        environment = os.environ.copy()
        repository = str(Path(__file__).resolve().parents[2])
        environment["PYTHONPATH"] = repository + os.pathsep + environment.get("PYTHONPATH", "")
        environment["PYTHONIOENCODING"] = "utf-8"
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        started = monotonic()
        completed = False
        failed = False
        usage: dict[str, Any] = {}
        runtime_thread = thread_id
        pending = ""
        process_options = _hidden_process_options()
        if os.name != "nt":
            process_options["start_new_session"] = True
        with prompt_path.open("rb") as source, output_path.open("wb") as output, \
                stderr_path.open("wb") as errors, output_path.open("r", encoding="utf-8", errors="replace") as reader:
            process = subprocess.Popen(
                self.command(workspace, turn_directory, thread_id), cwd=workspace,
                stdin=source, stdout=output, stderr=errors, env=environment,
                **process_options,
            )
            try:
                while True:
                    exit_code = process.poll()
                    if cancel.is_set():
                        raise AgentStopped("cancelled")
                    if monotonic() - started >= self.timeout_seconds:
                        raise AgentStopped("timed_out")
                    pending += reader.read()
                    lines = pending.split("\n")
                    pending = lines.pop()
                    for line in lines:
                        try:
                            event = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        if not isinstance(event, dict):
                            continue
                        on_event(event)
                        if event.get("type") == "thread.started":
                            runtime_thread = event.get("thread_id")
                        elif event.get("type") == "turn.completed":
                            completed = True
                            usage = event.get("usage", {})
                        elif event.get("type") == "turn.failed":
                            failed = True
                    if exit_code is not None:
                        break
                    cancel.wait(0.1)
            finally:
                _stop_process(process)
        if process.returncode or failed or not completed or not runtime_thread:
            detail = stderr_path.read_text("utf-8", errors="replace")[-2000:]
            raise RuntimeError(f"Codex turn did not complete (exit {process.returncode}). {detail}")
        answer_path = turn_directory / "agent-final.json"
        if not answer_path.is_file():
            raise RuntimeError("Codex completed without a structured final answer")
        return AgentResult(json.loads(answer_path.read_text("utf-8")), runtime_thread, usage)
