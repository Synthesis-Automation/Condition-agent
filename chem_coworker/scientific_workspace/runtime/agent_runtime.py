"""Optional agent harness adapters; deterministic scientific packages stay model-free."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from threading import Event
from tempfile import TemporaryDirectory
from time import monotonic
from typing import Any, Callable, Protocol

from ..answers.answer_contracts import ANSWER_SCHEMA
from ..answers.answer_handoff import (
    ANSWER_HANDOFF_SCHEMA,
    ANSWER_HANDOFF_VERSION,
    AnswerSubmissionError,
    answer_handoff_prompt,
    load_answer_handoff,
    load_text_answer,
    prepare_answer_handoff,
)
from ..core.process_utils import hidden_process_options, stop_process_tree
from ..paths import REPOSITORY_ROOT
from .research_profiles import resolve_research_profile
from .runtime_environment import discover_ripgrep, runtime_environment
from .workspace_modes import WorkspaceMode, mode_policy, runtime_configuration


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
        mode: str | WorkspaceMode = WorkspaceMode.NORMAL,
    ) -> AgentResult:
        """Run one complete agent turn with iterative tools and optional conversation memory."""


class AgentStopped(RuntimeError):
    """Cancellation or elapsed wall-clock budget; partial evidence remains on disk."""

    def __init__(self, status: str) -> None:
        super().__init__(status)
        self.status = status


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
                encoding="utf-8", errors="replace", timeout=10, **hidden_process_options(),
            )
            if help_result.returncode or not all(
                flag in help_result.stdout for flag in ("--output-schema", "resume", "--json")
            ):
                continue
            version = subprocess.run(
                [candidate, "--version"], capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=10, check=True,
                **hidden_process_options(),
            ).stdout.strip()
            return str(Path(candidate).resolve()), version
        except (OSError, subprocess.SubprocessError):
            continue
    raise RuntimeError(
        "A current native Codex CLI with exec, resume and --output-schema is required. "
        "Install/update Codex, run codex login, or set SCIENTIFIC_CODEX_PATH to its executable."
    )


def _runtime_error_message(value: Any) -> str:
    """Extract a bounded message from CLI errors, including encoded provider JSON."""
    for _ in range(6):
        if isinstance(value, dict):
            value = value.get("error") or value.get("message")
        elif isinstance(value, str):
            try:
                decoded = json.loads(value)
            except json.JSONDecodeError:
                return value.strip()[:2000]
            if not isinstance(decoded, dict):
                return value.strip()[:2000]
            value = decoded
        else:
            return ""
    return value.strip()[:2000] if isinstance(value, str) else ""


def _execution_setup_error(item: dict[str, Any]) -> str:
    """Recognize an execution setup failure from tool output, never agent prose."""
    if item.get("status") != "failed":
        return ""
    if item.get("type") == "command_execution":
        detail = item.get("aggregated_output", "")
    elif item.get("type") == "mcp_tool_call":
        detail = _runtime_error_message(item.get("error"))
    else:
        return ""
    if isinstance(detail, str) and "helper_unknown_error: setup refresh had errors" in detail:
        return detail.strip()[:2000]
    return ""


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
        self.search_tool = discover_ripgrep(runtime_executable=self.executable)
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
            "executable": getattr(self, "executable", None),
            "model": model or "Codex configuration default",
            "sandbox": "workspace-write", "timeout_seconds": self.timeout_seconds,
            "authentication": "existing Codex login; checked when a turn runs",
            "research_settings": settings.to_dict(),
            "model_request": {"value": model, "source": "override" if model else "inherited"},
            "effective_configuration_status": "unconfirmed",
            "capabilities": {"web_search": "not_checked", "code_execution": "not_checked"},
            "local_tools": {"rg": getattr(self, "search_tool", {"status": "not_checked"})},
            "answer_transport": {"schema_version": ANSWER_HANDOFF_VERSION, "answer_file": "answer-draft.json"},
            "limitations": [
                "Requested settings do not establish provider support or successful tool access.",
                "The effective model and reasoning effort are not resolved by this adapter.",
                "Disabling the web_search tool is not a network-isolation policy for shell or MCP tools.",
            ],
        }

    def command(
        self, workspace: Path, turn_directory: Path, thread_id: str | None,
        mode: str | WorkspaceMode = WorkspaceMode.NORMAL,
    ) -> list[str]:
        """Build explicit arguments; never interpolate a user's question into a shell."""
        mode = WorkspaceMode(mode)
        command = [
            self.executable, "exec", "--sandbox", "read-only" if mode == WorkspaceMode.PURE else "workspace-write",
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
        if not mode.task_guidance:
            # Prevent implicit repository AGENTS.md or personal memory from
            # reintroducing removed application instructions.
            for setting in ('project_doc_max_bytes=0', 'developer_instructions=""',
                            'features.memories=false', 'features.skip_host_skill_discovery=true',
                            'features.plugins=false'):
                command.extend(["-c", setting])
        if mode == WorkspaceMode.PURE:
            command.extend(self._pure_tool_overrides(workspace))
        command.extend(["--json", "--skip-git-repo-check"])
        if mode.structured_answer:
            command.extend(["--output-schema", str(turn_directory / "answer-handoff-schema.json")])
        command.extend(["--output-last-message", str(turn_directory / (
            "agent-final.json" if mode.structured_answer else "agent-final.txt"))])
        if self.model:
            command.extend(["--model", self.model])
        command.append("-")
        return command

    def _pure_tool_overrides(self, workspace: Path) -> list[str]:
        """Fail closed unless local tool paths can be disabled on this CLI."""
        disabled = ("shell_tool", "unified_exec", "apps", "plugins", "multi_agent",
                    "browser_use", "computer_use", "view_image",
                    "hooks", "memories", "image_generation", "skill_search",
                    "in_app_browser", "in_app_local_automation", "browser_use_external",
                    "multi_agent_v2", "remote_plugin")
        feature_result = subprocess.run(
            [self.executable, "features", "list"], cwd=workspace,
            capture_output=True, text=True, encoding="utf-8", timeout=10,
            check=True, **hidden_process_options(),
        )
        available = {line.split()[0] for line in feature_result.stdout.splitlines() if line.split()}
        required = {"shell_tool", "apps", "plugins", "skip_host_skill_discovery", "view_image"}
        if not required.issubset(available):
            raise RuntimeError("Pure agent mode requires a current Codex CLI with local-tool and skill isolation controls")
        overrides = []
        # Code mode dispatches native tools (including web search); it is not
        # local shell execution. Disabling its host also breaks allowed tools.
        # Restrict the underlying local-access tools, not their shared router.
        if "code_mode_host" in available:
            overrides.extend(["-c", "features.code_mode_host=true"])
        for feature in disabled:
            if feature in available:
                overrides.extend(["-c", f"features.{feature}=false"])
        # An empty table override can merge with inherited servers. Disable each
        # resolved server explicitly instead; never log its credentials/config.
        servers = subprocess.run(
            [self.executable, *overrides, "mcp", "list", "--json"], cwd=workspace,
            capture_output=True, text=True, encoding="utf-8", timeout=10,
            check=True, **hidden_process_options(),
        )
        records = json.loads(servers.stdout)
        if not isinstance(records, list):
            raise RuntimeError("Cannot establish MCP isolation for pure agent mode")
        for server in records:
            name = server["name"]
            # CLI override keys split on dots; TOML-quoting a component creates
            # a different key on supported CLI releases. Refuse ambiguous names.
            if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", name):
                raise RuntimeError("Cannot isolate an MCP server with an unsupported configuration name")
            overrides.extend(["-c", f"mcp_servers.{name}.enabled=false"])
        verified = subprocess.run(
            [self.executable, *overrides, "mcp", "list", "--json"], cwd=workspace,
            capture_output=True, text=True, encoding="utf-8", timeout=10,
            check=True, **hidden_process_options(),
        )
        remaining = json.loads(verified.stdout)
        if not isinstance(remaining, list) or any(server.get("enabled") is not False for server in remaining):
            raise RuntimeError("Pure agent mode could not disable all configured MCP servers")
        return overrides

    def run(
        self, *, prompt: str, workspace: Path, turn_directory: Path,
        thread_id: str | None, cancel: Event,
        on_event: Callable[[dict[str, Any]], None],
        mode: str | WorkspaceMode = WorkspaceMode.NORMAL,
    ) -> AgentResult:
        """Run one mode without mutating shared runtime configuration."""
        mode = WorkspaceMode(mode)
        arguments = dict(prompt=prompt, workspace=workspace, turn_directory=turn_directory,
                         thread_id=thread_id, cancel=cancel, on_event=on_event, mode=mode)
        if mode == WorkspaceMode.PURE:
            # Outside the checkout: no project configuration or ancestor files.
            with TemporaryDirectory(prefix="scientific-pure-") as scratch:
                return self._run_turn(**arguments, execution_directory=Path(scratch))
        return self._run_turn(**arguments, execution_directory=workspace)

    def _run_turn(
        self, *, prompt: str, workspace: Path, turn_directory: Path,
        thread_id: str | None, cancel: Event,
        on_event: Callable[[dict[str, Any]], None], mode: WorkspaceMode,
        execution_directory: Path,
    ) -> AgentResult:
        """Capture JSONL events, enforce a deadline, and preserve complete local logs."""
        turn_directory = prepare_answer_handoff(workspace, turn_directory)
        if mode.structured_answer:
            (turn_directory / "answer-schema.json").write_text(json.dumps(ANSWER_SCHEMA), "utf-8")
            (turn_directory / "answer-handoff-schema.json").write_text(json.dumps(ANSWER_HANDOFF_SCHEMA), "utf-8")
        (turn_directory / "runtime-request.json").write_text(json.dumps({
            "schema_version": "scientific_runtime_request.v1",
            "configuration": runtime_configuration(self.describe(), mode), "resumed_thread_id": thread_id,
            "mode_policy": mode_policy(mode),
        }, indent=2), "utf-8")
        prompt_path = turn_directory / "prompt.txt"
        prompt_path.write_text(answer_handoff_prompt(prompt, turn_directory)
                               if mode.structured_answer else prompt, "utf-8")
        output_path = turn_directory / "runtime.jsonl"
        stderr_path = turn_directory / "runtime.stderr.txt"
        environment = runtime_environment(
            None if mode == WorkspaceMode.PURE else REPOSITORY_ROOT,
            search_tool={} if mode == WorkspaceMode.PURE else getattr(self, "search_tool", {}),
        )
        started = monotonic()
        last_heartbeat = started
        completed = False
        failed = False
        usage: dict[str, Any] = {}
        runtime_thread = thread_id
        pending = ""
        observed_items: dict[str, dict[str, Any]] = {}
        failure_detail = ""
        event_error = ""
        item_error = ""
        execution_setup_error = ""
        process_options = hidden_process_options()
        if os.name != "nt":
            process_options["start_new_session"] = True
        with prompt_path.open("rb") as source, output_path.open("wb") as output, \
                stderr_path.open("wb") as errors, output_path.open("r", encoding="utf-8", errors="replace") as reader:
            process = subprocess.Popen(
                (self.command(execution_directory, turn_directory, thread_id)
                 if mode == WorkspaceMode.NORMAL else
                 self.command(execution_directory, turn_directory, thread_id, mode)), cwd=execution_directory,
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
                        if event.get("type") == "item.completed" and isinstance(item, dict):
                            setup_error = _execution_setup_error(item)
                            if setup_error:
                                execution_setup_error = setup_error
                            elif (
                                item.get("type") == "command_execution"
                                and item.get("exit_code") == 0
                                and item.get("status") == "completed"
                            ):
                                execution_setup_error = ""
                        if event.get("type") in {"item.started", "item.updated", "item.completed"} and isinstance(item, dict):
                            identity = item.get("id")
                            if isinstance(identity, str):
                                observed_items[identity] = {
                                    **observed_items.get(identity, {}), **item,
                                    "last_event_type": event["type"],
                                }
                            if item.get("type") == "error":
                                item_error = _runtime_error_message(item) or item_error
                        on_event(event)
                        if event.get("type") == "thread.started":
                            runtime_thread = event.get("thread_id")
                        elif event.get("type") == "turn.completed":
                            completed = True
                            usage = event.get("usage", {})
                        elif event.get("type") == "turn.failed":
                            failed = True
                            failure_detail = _runtime_error_message(event) or failure_detail
                        elif event.get("type") == "error":
                            event_error = _runtime_error_message(event) or event_error
                    if exit_code is not None:
                        break
                    if monotonic() - last_heartbeat >= 1:
                        # An internal polling signal, not an invented agent update.
                        # Consumers can observe scientific artifacts created by
                        # a still-running shell command without waiting for exit.
                        on_event({"type": "runtime.heartbeat"})
                        last_heartbeat = monotonic()
                    cancel.wait(0.1)
            finally:
                stop_process_tree(process)
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
                    "runtime_error": failure_detail or event_error or (
                        item_error if process.returncode or not completed else None
                    ),
                    "execution_setup_error": execution_setup_error or None,
                    "tool_events": counts,
                    "limitations": [
                        "Counts reflect unique tool item IDs in observed JSONL events, not successful scientific checks.",
                        "A completed event does not establish successful access or the correctness of its result.",
                        "Missing events do not establish that a capability was unavailable.",
                        "Effective model, reasoning effort, and search mode are not inferred from these counts.",
                    ],
                }, indent=2), "utf-8")
        if process.returncode or failed or not completed or not runtime_thread:
            detail = failure_detail or event_error or item_error or stderr_path.read_text(
                "utf-8", errors="replace",
            )[-2000:].strip()
            message = f"Codex turn did not complete (exit {process.returncode})."
            if detail:
                message += f" {detail}"
            if "model" in detail.lower() and "not supported" in detail.lower():
                message += (
                    " Restart the scientific-chat server with --agent-model MODEL_ID"
                    " available to this Codex login, or use --codex PATH to select"
                    " a CLI that supports the requested model."
                )
            raise RuntimeError(message)
        if mode.structured_answer:
            try:
                answer = load_answer_handoff(workspace, turn_directory, thread_id=runtime_thread, usage=usage)
            except AnswerSubmissionError as exc:
                if (execution_setup_error and exc.payload.get("issue") == "FileNotFoundError"
                        and not (turn_directory / "answer-draft.json").exists()):
                    raise RuntimeError(
                        "Local execution failed during Windows sandbox setup; no answer was submitted. "
                        "Answer correction cannot repair this runtime failure. "
                        f"{execution_setup_error} "
                        "Inspect the Codex .sandbox logs for the underlying setup error, "
                        "restore local execution, then retry the investigation."
                    ) from exc
                raise
        else:
            answer = load_text_answer(workspace, turn_directory)
        return AgentResult(answer, runtime_thread, usage)
