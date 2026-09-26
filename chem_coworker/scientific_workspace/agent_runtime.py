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
from .research_profiles import resolve_research_profile


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
        profile: str = "inherit", reasoning_effort: str | None = None,
        web_search: str | None = None, timeout_seconds: float | None = None,
    ) -> None:
        self.settings = resolve_research_profile(
            profile, reasoning_effort=reasoning_effort,
            web_search=web_search, timeout_seconds=timeout_seconds,
        )
        if model is not None and (not isinstance(model, str) or not model.strip()):
            raise ValueError("model must be a nonempty identifier or None to inherit")
        self.executable, self.version = find_codex(executable)
        self.model = model
        self.timeout_seconds = self.settings.timeout_seconds

    def describe(self) -> dict[str, Any]:
        """Local discovery does not assert provider authentication or model availability."""
        model = getattr(self, "model", None)
        settings = getattr(self, "settings", None)
        if settings is None:
            settings = resolve_research_profile(timeout_seconds=self.timeout_seconds)
        return {
            "runtime": "codex_exec", "version": getattr(self, "version", "unreported"),
            "model": model or "Codex configuration default",
            "sandbox": "workspace-write", "timeout_seconds": self.timeout_seconds,
            "authentication": "existing Codex login; checked when a turn runs",
            "research_settings": settings.to_dict(),
            "model_request": {"value": model, "source": "override" if model else "inherited"},
            "effective_configuration_status": "unconfirmed",
            "capabilities": {"web_search": "not_checked", "code_execution": "not_checked"},
            "limitations": [
                "Requested settings do not establish provider support or successful tool access.",
                "The effective model and reasoning effort are not resolved by this adapter.",
                "Disabling the web_search tool is not a network-isolation policy for shell or MCP tools.",
            ],
        }

    def command(self, workspace: Path, turn_directory: Path, thread_id: str | None) -> list[str]:
        """Build explicit arguments; never interpolate a user's question into a shell."""
        command = [
            self.executable, "exec", "--sandbox", "workspace-write",
            "-c", 'approval_policy="never"', "--cd", str(workspace),
        ]
        if thread_id:
            command.extend(["resume", thread_id])
        settings = getattr(self, "settings", None)
        if settings is not None:
            for key, value in (("model_reasoning_effort", settings.reasoning_effort),
                               ("web_search", settings.web_search)):
                if value is not None:
                    # This is a direct argv item, not shell-escaped command text.
                    command.extend(["-c", f"{key}={json.dumps(value)}"])
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
        (turn_directory / "runtime-request.json").write_text(json.dumps({
            "schema_version": "scientific_runtime_request.v1",
            "configuration": self.describe(), "resumed_thread_id": thread_id,
        }, indent=2), "utf-8")
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
        observed_items: dict[str, dict[str, Any]] = {}
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
                    if exit_code is not None and pending:
                        lines.append(pending)
                        pending = ""
                    for line in lines:
                        try:
                            event = json.loads(line)
                        except json.JSONDecodeError:
                            continue
                        if not isinstance(event, dict):
                            continue
                        item = event.get("item")
                        if event.get("type") in {"item.started", "item.updated", "item.completed"} and isinstance(item, dict):
                            identity = item.get("id")
                            if isinstance(identity, str):
                                observed_items[identity] = {
                                    **observed_items.get(identity, {}), **item,
                                    "last_event_type": event["type"],
                                }
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
                counts = {}
                for kind in ("web_search", "command_execution", "mcp_tool_call"):
                    items = [item for item in observed_items.values() if item.get("type") == kind]
                    counts[kind] = {
                        "observed": len(items),
                        "completed_events": sum(item["last_event_type"] == "item.completed" for item in items),
                        "reported_failed": sum(item.get("status") == "failed" for item in items),
                    }
                (turn_directory / "runtime-observations.json").write_text(json.dumps({
                    "schema_version": "scientific_runtime_observations.v1",
                    "thread_id": runtime_thread,
                    "elapsed_seconds": round(monotonic() - started, 3),
                    "process_exit_code": process.returncode,
                    "turn_completed_event": completed, "turn_failed_event": failed,
                    "tool_events": counts,
                    "limitations": [
                        "Counts reflect unique tool item IDs in observed JSONL events, not successful scientific checks.",
                        "A completed event does not establish successful access or the correctness of its result.",
                        "Missing events do not establish that a capability was unavailable.",
                        "Effective model, reasoning effort, and search mode are not inferred from these counts.",
                    ],
                }, indent=2), "utf-8")
        if process.returncode or failed or not completed or not runtime_thread:
            detail = stderr_path.read_text("utf-8", errors="replace")[-2000:]
            raise RuntimeError(f"Codex turn did not complete (exit {process.returncode}). {detail}")
        answer_path = turn_directory / "agent-final.json"
        if not answer_path.is_file():
            raise RuntimeError("Codex completed without a structured final answer")
        return AgentResult(json.loads(answer_path.read_text("utf-8")), runtime_thread, usage)
