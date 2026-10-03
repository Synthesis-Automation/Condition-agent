"""Persistent local conversations that compose an agent runtime and scientific workspace."""

from __future__ import annotations

import json
import os
import re
import traceback
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from threading import Event, Lock
from typing import Any, Mapping
from uuid import uuid4

from ..answers.answer_contracts import ScientificAnswer, validate_answer_evidence
from ..answers.answer_handoff import AnswerSubmissionError
from ..answers.evidence_review import matching_evidence_review
from ..core.baseline import sha256_file, verify_baseline
from ..core.store import _read_json, _write_json
from ..agent_context.context import record_application_context
from ..agent_context.prompts import investigation_prompt
from ..paths import REPOSITORY_ROOT
from ..workspace import ScientificWorkspace
from .activity import ACTIVITY_VERSION, ActivityHistory, recover_activity
from .agent_runtime import AgentResult, AgentRuntime, AgentStopped
from .workspace_modes import WorkspaceMode, available_modes, mode_policy, runtime_configuration


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _identifier(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{32}", value):
        raise ValueError("Invalid conversation or turn identifier")
    return value


class ConversationService:
    """One active local turn, saved conversations, explicit cancellation and follow-up."""

    def __init__(
        self, root: str | Path, *, runtime: AgentRuntime,
        repository: str | Path | None = None,
        artifacts: Mapping[str, str | Path] | None = None,
    ) -> None:
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.repository = Path(repository or REPOSITORY_ROOT).resolve()
        self.artifacts = dict(artifacts or {})
        self.runtime = runtime
        self._activity_cache: dict[Path, tuple[Any, list[dict[str, Any]]]] = {}
        self._mutex = Lock()
        self._active: tuple[str, str, Event] | None = None
        self._pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="scientific-chat")

    def describe(self) -> dict[str, Any]:
        """Report runtime and configured data availability without loading large indexes."""
        return {
            **self.runtime.describe(), "local_only": True,
            "validation_status": "development_snapshot_not_release_validated",
            "artifacts": {name: (self.repository / path).is_file() for name, path in self.artifacts.items()},
            "workspace_modes": available_modes(), "default_workspace_mode": WorkspaceMode.NORMAL.value,
        }

    def _directory(self, conversation_id: str) -> Path:
        directory = (self.root / _identifier(conversation_id)).resolve()
        if directory.parent != self.root:
            raise ValueError("Conversation path escapes configured root")
        return directory

    def submit(
        self, question: str, conversation_id: str | None = None,
        *, mode: str | WorkspaceMode | None = None,
    ) -> dict[str, str]:
        """Queue a user turn; baseline creation happens in the worker with visible progress."""
        question = question.strip()
        if not question or len(question) > 20000:
            raise ValueError("Question must contain 1 to 20000 characters")
        with self._mutex:
            if self._active is not None:
                raise RuntimeError("An investigation is running; wait or cancel it first")
            identity = conversation_id or uuid4().hex
            directory = self._directory(identity)
            if conversation_id:
                if not (directory / "conversation.json").is_file():
                    raise FileNotFoundError("Conversation does not exist")
                saved = _read_json(directory / "conversation.json")
                selected = WorkspaceMode(saved.get("workspace_mode", "normal"))
                if mode is not None and WorkspaceMode(mode) != selected:
                    raise ValueError("Workspace mode is fixed per conversation; start a new chat to change it")
                if selected != WorkspaceMode.PURE and not (directory / "investigation.json").is_file():
                    raise ValueError("Investigation preparation failed; start a new conversation")
            else:
                selected = WorkspaceMode(WorkspaceMode.NORMAL if mode is None else mode)
                directory.mkdir(exist_ok=False)
                _write_json(directory / "conversation.json", {
                    "id": identity, "title": question[:100], "created_at": _now(),
                    "runtime": runtime_configuration(self.runtime.describe(), selected),
                    "workspace_mode": selected.value, "mode_policy": mode_policy(selected),
                })
                (directory / "turns").mkdir()
            lock_path = directory / ".conversation.lock"
            try:
                descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError as exc:
                raise RuntimeError("Conversation is locked; inspect an interrupted worker before removing its lock") from exc
            with os.fdopen(descriptor, "w") as lock:
                lock.write(str(os.getpid()))
            turn_id = uuid4().hex
            turn = directory / "turns" / turn_id
            turn.mkdir()
            _write_json(turn / "turn.json", {
                "id": turn_id, "conversation_id": identity, "question": question,
                "status": "queued", "created_at": _now(), "progress": [],
                "worker_process_id": os.getpid(),
                "runtime_requested": runtime_configuration(self.runtime.describe(), selected),
                "workspace_mode": selected.value, "mode_policy": mode_policy(selected),
            })
            cancel = Event()
            self._active = (identity, turn_id, cancel)
            self._pool.submit(self._run, identity, turn_id, cancel)
            return {"conversation_id": identity, "turn_id": turn_id}

    def _run(self, identity: str, turn_id: str, cancel: Event) -> None:
        from .activity import ScientificActivityCursor

        directory = self._directory(identity)
        turn = directory / "turns" / turn_id
        state = _read_json(turn / "turn.json")
        mode = WorkspaceMode(state.get("workspace_mode", "normal"))
        workspace: ScientificWorkspace | None = None
        history = ActivityHistory()
        attempt = 0
        scientific_cursor: ScientificActivityCursor | None = None
        last_activity_error: str | None = None

        def debug(kind: str, **details: Any) -> None:
            record = {"schema_version": "scientific_progress_log.v1", "at": _now(),
                      "turn_id": turn_id, "attempt": attempt, "kind": kind, **details}
            with (turn / "progress.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(record, ensure_ascii=False) + "\n")

        def save(**updates: Any) -> None:
            previous_status = state["status"]
            state.update(updates, updated_at=_now())
            _write_json(turn / "turn.json", state)
            if state["status"] != previous_status:
                debug("turn_status", status=state["status"])

        def scientific_progress() -> bool:
            nonlocal last_activity_error
            if scientific_cursor is None:
                return False
            try:
                rows = scientific_cursor.drain(history)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                # Observability failure must not interrupt the scientific work.
                message = f"{type(exc).__name__}: {exc}"
                if message != last_activity_error:
                    debug("activity_log_error", error=message)
                    last_activity_error = message
                return False
            last_activity_error = None
            for row in rows:
                debug(row["kind"], activity=row, event_sequence=row["event_sequence"],
                      artifact_ref=row["artifact_ref"])
            if rows:
                state["progress"] = history.rows
                state["activity_version"] = ACTIVITY_VERSION
            return bool(rows)

        def progress(event: dict[str, Any]) -> None:
            changed = scientific_progress()
            # Heartbeats only observe committed events; they never impersonate
            # agent commentary or move the timestamp of the last actual activity.
            if event.get("type") == "runtime.heartbeat":
                if changed:
                    save()
                return
            if event.get("type") == "thread.started":
                state["thread_id"] = event.get("thread_id")
            row = history.observe(event, _now(), scope=str(attempt))
            if row is not None:
                debug(row["kind"], event_type=event["type"], activity=row,
                      runtime_log="runtime.jsonl" if attempt == 0 else "repair-1/runtime.jsonl")
            elif event.get("type") in {"error", "turn.failed"}:
                debug("runtime_error", event_type=event["type"], error=event.get("error") or event.get("message"))
            state["progress"] = history.rows
            state["activity_version"] = ACTIVITY_VERSION
            save()

        try:
            save(status="preparing")
            if mode == WorkspaceMode.PURE:
                if cancel.is_set():
                    raise AgentStopped("cancelled")
                prior = next((item for item in reversed(self.get(identity)["turns"])
                              if item["id"] != turn_id and item.get("thread_id")), None)
                same_runtime = prior is not None and prior.get("runtime_requested") == state["runtime_requested"]
                thread_id = prior["thread_id"] if same_runtime else None
                prompt = state["question"]
                prompt_path = turn / "investigation-prompt.txt"
                prompt_path.write_bytes(prompt.encode("utf-8"))
                save(status="running", prompt_sha256=sha256_file(prompt_path),
                     thread_resume="resumed" if thread_id else "new_thread")
                result = self.runtime.run(
                    prompt=prompt, workspace=directory, turn_directory=turn,
                    thread_id=thread_id, cancel=cancel, on_event=progress, mode=mode,
                )
                if cancel.is_set():
                    raise AgentStopped("cancelled")
                save(status="completed", thread_id=result.thread_id,
                     answer=self._free_answer(result, mode))
                return
            if not (directory / "investigation.json").exists():
                # The conversation envelope already exists. Initialize the scientific
                # store separately, then move only its newly created files into place.
                from ..core.baseline import capture_baseline
                from ..core.store import InvestigationStore

                paths = {name: (self.repository / path).resolve() for name, path in self.artifacts.items()}
                baseline = (capture_baseline(self.repository, paths) if mode.task_guidance
                            else capture_baseline(self.repository, paths, include_guidance=False))
                staging = directory / "prepared"
                InvestigationStore.create(
                    staging, objective=state["question"], baseline=baseline,
                    agent_metadata={**self.runtime.describe(), "workspace_mode": mode.value},
                )
                for path in staging.iterdir():
                    path.rename(directory / path.name)
                staging.rmdir()
            workspace = ScientificWorkspace(directory)
            if cancel.is_set():
                raise AgentStopped("cancelled")
            verify_baseline(workspace.store.manifest["baseline"])
            if workspace.store.summary()["status"] != "active":
                raise ValueError("Investigation is stopped; resume its lifecycle explicitly first")
            prior_turns = [item for item in self.get(identity)["turns"] if item["id"] != turn_id]
            prior = next((item for item in reversed(prior_turns) if item.get("thread_id")), None)
            previous_config = prior.get("runtime_requested") if prior else None
            if prior and previous_config is None:
                previous_config = _read_json(directory / "conversation.json").get("runtime")
            changed_runtime = prior is not None and previous_config != state["runtime_requested"]
            context_event = record_application_context(
                workspace.store, turn_id=turn_id, runtime=state["runtime_requested"],
            )
            context = workspace.store.read_artifact(context_event.artifact_ref)
            previous_context = (
                workspace.store.read_artifact(prior["application_context_ref"])
                if prior and prior.get("application_context_ref") else None
            )
            changed_context = prior is not None and (
                previous_context is None or previous_context["layers"] != context["layers"]
            )
            thread_id = prior["thread_id"] if prior and not changed_runtime and not changed_context else None
            save(thread_resume=("configuration_changed_new_thread" if changed_runtime else
                                "application_context_changed_new_thread" if changed_context else
                                "resumed" if thread_id else "new_thread"),
                 application_context_ref=context_event.artifact_ref,
                 scientific_identity=context["scientific_identity"])
            _write_json(turn / "capabilities.json", workspace.capabilities())
            user_event = workspace.store.append("user_message", {
                "turn_id": turn_id, "text": state["question"], "origin": "user",
            })
            scientific_cursor = ScientificActivityCursor(workspace.store, after_sequence=user_event.sequence)
            save(status="running", user_ref=user_event.artifact_ref)
            prompt = investigation_prompt(workspace, state["question"], mode=mode)
            attempt_usage = []
            for attempt in range(2):
                attempt_directory = turn if attempt == 0 else turn / "repair-1"
                attempt_directory.mkdir(exist_ok=True)
                if cancel.is_set():
                    raise AgentStopped("cancelled")
                verify_baseline(workspace.store.manifest["baseline"])
                prompt_path = attempt_directory / "investigation-prompt.txt"
                prompt_path.write_bytes(prompt.encode("utf-8"))
                save(prompt_sha256=sha256_file(prompt_path))
                submission_error = None
                try:
                    result = self.runtime.run(
                        prompt=prompt, workspace=directory, turn_directory=attempt_directory,
                        thread_id=thread_id, cancel=cancel, on_event=progress,
                        **({"mode": mode} if mode != WorkspaceMode.NORMAL else {}),
                    )
                    submitted_answer = result.answer
                    submitted_thread = result.thread_id
                    submitted_usage = result.usage
                except AnswerSubmissionError as exc:
                    # A completed runtime can submit a missing/malformed answer
                    # file. Preserve it for the same single correction as an
                    # invalid scientific answer; execution failures still escape.
                    submission_error = exc
                    submitted_answer = exc.payload
                    submitted_thread = exc.thread_id
                    submitted_usage = exc.usage
                if scientific_progress():
                    save()
                attempt_usage.append(submitted_usage)
                if cancel.is_set():
                    raise AgentStopped("cancelled")
                verify_baseline(workspace.store.manifest["baseline"])
                if not mode.structured_answer:
                    if submission_error is not None:
                        raise submission_error
                    event = workspace.store.append("agent_answer", {
                        **self._free_answer(result, mode), "turn_id": turn_id,
                        "application_context_ref": context_event.artifact_ref,
                        "scientific_identity": context["scientific_identity"],
                    })
                    save(status="completed", answer_ref=event.artifact_ref,
                         answer=workspace.store.read_artifact(event.artifact_ref), thread_id=result.thread_id)
                    return
                try:
                    if submission_error is not None:
                        raise submission_error
                    answer = ScientificAnswer.model_validate(submitted_answer)
                    cited = validate_answer_evidence(answer, workspace.store)
                    break
                except (ValueError, FileNotFoundError) as exc:
                    debug("answer_validation_error", error={"type": type(exc).__name__, "message": str(exc)})
                    workspace.store.append("agent_answer_rejected", {
                        "turn_id": turn_id, "attempt": attempt + 1, "answer": submitted_answer,
                        "error": str(exc), "usage": submitted_usage,
                    })
                    if attempt:
                        raise
                    rejected_path = turn / "rejected-answer.json"
                    _write_json(rejected_path, submitted_answer)
                    thread_id = submitted_thread
                    save(repair_attempts=1, validation_error=str(exc))
                    prompt = investigation_prompt(workspace, state["question"], mode=mode) + (
                        f"\nYour previous submission was rejected: {str(exc)[:4000]}\n"
                        f"Rejected answer or transport diagnostics: {rejected_path}. "
                        f"The previous attempt's draft, if created, is at {attempt_directory / 'answer-draft.json'}.\n"
                        "This is the only correction attempt. Inspect saved evidence, correct the draft, "
                        "and save/validate it at this attempt's runtime-provided path before submitting "
                        "through the runtime handoff. Diagnostics are not an answer draft. Never fabricate evidence, "
                        "weaken checks, or label a proposal as an observation to make validation pass."
                    )
            answer.evidence_refs = cited
            review_ref = matching_evidence_review(workspace.store, answer, after_sequence=user_event.sequence)
            trace = {
                path.relative_to(turn).as_posix(): sha256_file(path) for path in turn.rglob("*")
                # Operational logs continue through final status persistence.
                if path.is_file() and path.name not in {"turn.json", "progress.jsonl"}
            }
            event = workspace.store.append("agent_answer", {
                **answer.model_dump(), "turn_id": turn_id, "thread_id": result.thread_id,
                "origin": "agent_authored", "review_status": "unreviewed",
                "workspace_mode": mode.value, "mode_policy": mode_policy(mode),
                "evidence_status": "linked_unreviewed" if cited else "no_local_evidence",
                "runtime": runtime_configuration(self.runtime.describe(), mode), "usage": result.usage,
                "attempt_usage": attempt_usage,
                "application_context_ref": context_event.artifact_ref,
                "scientific_identity": context["scientific_identity"],
                "self_review_status": "recorded_for_final_draft" if review_ref else "not_recorded_for_final_draft",
                "self_review_ref": review_ref,
                "trace_files_sha256": trace,
            }, evidence_refs=tuple(answer.evidence_refs))
            save(
                status="completed", answer_ref=event.artifact_ref,
                answer=workspace.store.read_artifact(event.artifact_ref),
                thread_id=result.thread_id,
            )
        except Exception as exc:
            scientific_progress()
            status = exc.status if isinstance(exc, AgentStopped) else "failed"
            error = {"type": type(exc).__name__, "message": str(exc)}
            debug("turn_error", status=status, error=error, traceback=traceback.format_exc())
            if workspace is not None:
                try:
                    workspace.store.append("agent_turn_error", {"turn_id": turn_id, "status": status, "error": error})
                except (OSError, ValueError, RuntimeError):
                    pass
            save(status=status, error=error)
        finally:
            if scientific_progress():
                save()
            if mode.task_guidance and workspace is not None and workspace.store.manifest["baseline"].get("learning_context"):
                try:
                    published = workspace.publish_lessons()
                    if published["published"]:
                        debug("lessons_published", **published)
                except Exception as exc:
                    # Advice publication must not change the scientific answer or
                    # conceal a completed/failed turn; recorded lessons remain retryable.
                    debug("lesson_publication_error", error={"type": type(exc).__name__, "message": str(exc)})
            (directory / ".conversation.lock").unlink(missing_ok=True)
            with self._mutex:
                self._active = None

    def _free_answer(self, result: AgentResult, mode: WorkspaceMode) -> dict[str, Any]:
        """Wrap verbatim text for storage, without scientific validation or repair."""
        text = result.answer.get("answer_markdown")
        if not isinstance(text, str) or not text.strip():
            raise ValueError("Agent returned no final text")
        return {"schema_version": "agent_text.v1", "answer_markdown": text,
                "workspace_mode": mode.value, "mode_policy": mode_policy(mode),
                "origin": "agent_authored", "review_status": "unreviewed",
                "evidence_status": "not_validated", "thread_id": result.thread_id,
                "runtime": runtime_configuration(self.runtime.describe(), mode), "usage": result.usage}

    def get(self, conversation_id: str) -> dict[str, Any]:
        """Read persisted conversation and progress; safe to reopen after a normal restart."""
        directory = self._directory(conversation_id)
        metadata = json.loads((directory / "conversation.json").read_text("utf-8"))
        turns = [_read_json(path) for path in (directory / "turns").glob("*/turn.json")]
        turns.sort(key=lambda item: (item["created_at"], item["id"]))
        for turn in turns:
            turn_directory = directory / "turns" / _identifier(turn["id"])
            turn["debug_log_available"] = (turn_directory / "progress.jsonl").is_file()
            if turn.get("activity_version") != ACTIVITY_VERSION:
                paths = [turn_directory / "runtime.jsonl", turn_directory / "repair-1" / "runtime.jsonl"]
                signature = tuple((str(path), path.stat().st_mtime_ns, path.stat().st_size)
                                  for path in [turn_directory / "turn.json", *paths] if path.is_file())
                cached = self._activity_cache.get(turn_directory)
                if cached is None or cached[0] != signature:
                    rows = recover_activity(paths, turn.get("progress", []))
                    if len(self._activity_cache) >= 32:
                        self._activity_cache.clear()
                    cached = (signature, rows)
                    self._activity_cache[turn_directory] = cached
                turn["progress"] = cached[1]
            if turn["status"] in {"queued", "preparing", "running"} and (
                self._active is None or self._active[:2] != (conversation_id, turn["id"])
            ):
                # A different worker is never silently adopted or killed. The
                # persisted lock must be inspected after an abnormal shutdown.
                turn["status"] = "interrupted"
                turn["error"] = {
                    "type": "WorkerUnavailable",
                    "message": "This server does not own the previous worker. Start a new investigation, or inspect the saved worker PID and conversation lock before resuming.",
                }
            if (turn.get("answer") or {}).get("schema_version") == "scientific_answer.v2":
                from ..views.precedents import answer_step_precedents
                from ..views.step_assessments import answer_step_assessments

                store = ScientificWorkspace(directory).store
                turn["step_precedent_evidence"] = answer_step_precedents(
                    store, turn["answer"],
                )
                turn["step_assessment_evidence"] = answer_step_assessments(store, turn["answer"])
        return {"workspace_mode": "normal", **metadata, "turns": turns}

    def list_conversations(self) -> list[dict[str, Any]]:
        """List saved conversation titles, newest first."""
        items = [json.loads(path.read_text("utf-8")) for path in self.root.glob("*/conversation.json")]
        return sorted(items, key=lambda item: item["created_at"], reverse=True)

    def artifact(self, conversation_id: str, reference: str, *, expanded: bool = False) -> Any:
        """Read verified compact evidence, optionally restoring the full payload."""
        return ScientificWorkspace(self._directory(conversation_id)).store.read_artifact(reference, expanded=expanded)

    def debug_log(self, conversation_id: str, turn_id: str) -> bytes:
        """Snapshot complete debug-log lines while the worker may still append."""
        directory = self._directory(conversation_id)
        path = (directory / "turns" / _identifier(turn_id) / "progress.jsonl").resolve()
        if not path.is_relative_to(directory):
            raise ValueError("Debug log path escapes this conversation")
        data = path.read_bytes()
        return data[:data.rfind(b"\n") + 1]

    def cancel(self, conversation_id: str) -> bool:
        """Signal the active runtime; partial scientific artifacts are preserved."""
        _identifier(conversation_id)
        with self._mutex:
            if self._active and self._active[0] == conversation_id:
                self._active[2].set()
                return True
        return False

    def close(self) -> None:
        """Cancel an active runtime and wait for its worker to preserve terminal state."""
        with self._mutex:
            if self._active:
                self._active[2].set()
        self._pool.shutdown(wait=True)
