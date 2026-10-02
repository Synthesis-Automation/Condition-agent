"""Native process answer-file transport and rejected-output boundary regressions."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from threading import Event

import pytest

from chem_coworker.scientific_workspace.answers.answer_handoff import AnswerSubmissionError
from chem_coworker.scientific_workspace.runtime.agent_runtime import CodexRuntime
from chem_coworker.scientific_workspace.answers.answer_contracts import ANSWER_SCHEMA
from chem_coworker.scientific_workspace.answers.answer_handoff import (
    ANSWER_HANDOFF_SCHEMA,
    ANSWER_HANDOFF_VERSION,
)


HANDOFF = {"schema_version": ANSWER_HANDOFF_VERSION, "answer_file": "answer-draft.json"}
ANSWER = {"answer_markdown": "A source-linked route proposal", "evidence_refs": [],
          "uncertainties": ["Not experimentally verified"], "needs_user_input": False}


def fake_runtime(tmp_path: Path, script_text: str) -> tuple[CodexRuntime, Path]:
    """Use a real local child process, never a model or network service."""
    attempt = tmp_path / "turns" / "attempt"
    attempt.mkdir(parents=True)
    script = tmp_path / "fake_runtime.py"
    script.write_text(
        "import json, os, pathlib, sys\n"
        "attempt = pathlib.Path(sys.argv[1])\n"
        "print(json.dumps({'type':'thread.started', 'thread_id':'saved-thread'}), flush=True)\n"
        + script_text + "\n"
        "print(json.dumps({'type':'turn.completed', 'usage':{'input_tokens':7}}), flush=True)\n",
        "utf-8",
    )
    runtime = object.__new__(CodexRuntime)
    runtime.timeout_seconds = 10
    runtime.model = None
    runtime.command = lambda *_: [sys.executable, str(script), str(attempt)]
    return runtime, attempt


def submit(runtime: CodexRuntime, workspace: Path, attempt: Path):
    return runtime.run(prompt="Investigate the target.", workspace=workspace, turn_directory=attempt,
                       thread_id=None, cancel=Event(), on_event=lambda event: None)


def write_output(name: str, text: str) -> str:
    return f"(attempt / {name!r}).write_text({text!r}, encoding='utf-8')\n"


def test_success_reads_draft_once_and_preserves_full_schema_for_agent(tmp_path):
    runtime, attempt = fake_runtime(
        tmp_path, write_output("answer-draft.json", json.dumps(ANSWER))
        + write_output("agent-final.json", json.dumps(HANDOFF)),
    )
    result = submit(runtime, tmp_path, attempt)
    assert result.answer == ANSWER
    assert result.thread_id == "saved-thread" and result.usage == {"input_tokens": 7}
    assert json.loads((attempt / "answer-schema.json").read_text("utf-8")) == ANSWER_SCHEMA
    assert json.loads((attempt / "answer-handoff-schema.json").read_text("utf-8")) == ANSWER_HANDOFF_SCHEMA
    prompt = (attempt / "prompt.txt").read_text("utf-8")
    assert prompt.startswith("Investigate the target.")
    assert json.dumps(str(attempt / "answer-draft.json")) in prompt
    assert json.dumps(HANDOFF) in prompt
    assert runtime.describe()["answer_transport"]["schema_version"] == ANSWER_HANDOFF_VERSION
    # Later file changes cannot change the already loaded result used by validation.
    (attempt / "answer-draft.json").write_text("{}", "utf-8")
    assert result.answer == ANSWER


def test_cli_uses_small_receipt_schema_for_both_new_and_resumed_threads(tmp_path):
    runtime = object.__new__(CodexRuntime)
    runtime.executable, runtime.model = "codex", None
    for thread in (None, "saved-thread"):
        command = runtime.command(tmp_path, tmp_path, thread)
        assert command[command.index("--output-schema") + 1] == str(tmp_path / "answer-handoff-schema.json")
        assert command[command.index("--output-last-message") + 1] == str(tmp_path / "agent-final.json")


@pytest.mark.parametrize("handoff_text", [
    "{", "[]", json.dumps({**HANDOFF, "answer_file": "../other.json"}),
    json.dumps({**HANDOFF, "answer_file": "/outside/answer-draft.json"}),
    json.dumps({**HANDOFF, "answer_file": "C:\\outside\\answer-draft.json"}),
    json.dumps({**HANDOFF, "schema_version": "scientific_answer.v2"}),
    json.dumps({**HANDOFF, "extra": "value"}),
    '{"schema_version":"scientific_answer_handoff.v1","answer_file":"outside","answer_file":"answer-draft.json"}',
])
def test_rejects_malformed_or_redirected_receipt_without_following_its_path(tmp_path, handoff_text):
    runtime, attempt = fake_runtime(
        tmp_path, write_output("answer-draft.json", json.dumps(ANSWER))
        + write_output("agent-final.json", handoff_text),
    )
    with pytest.raises(AnswerSubmissionError) as caught:
        submit(runtime, tmp_path, attempt)
    failure = caught.value
    assert failure.thread_id == "saved-thread" and failure.usage == {"input_tokens": 7}
    assert failure.payload["stage"] == "handoff"
    assert (attempt / "agent-final.json").read_text("utf-8") == handoff_text
    observed = json.loads((attempt / "runtime-observations.json").read_text("utf-8"))
    assert observed["turn_completed_event"] and observed["process_exit_code"] == 0


@pytest.mark.parametrize("draft_text", ["{", "[]", "null", '{"value":NaN}', '{"x":1,"x":2}'])
def test_rejects_invalid_or_nonobject_draft_for_service_correction(tmp_path, draft_text):
    runtime, attempt = fake_runtime(
        tmp_path, write_output("answer-draft.json", draft_text)
        + write_output("agent-final.json", json.dumps(HANDOFF)),
    )
    with pytest.raises(AnswerSubmissionError) as caught:
        submit(runtime, tmp_path, attempt)
    assert caught.value.payload["stage"] == "draft"
    assert caught.value.thread_id == "saved-thread"
    assert (attempt / "answer-draft.json").read_text("utf-8") == draft_text


@pytest.mark.parametrize("missing,stage", [("agent-final.json", "handoff"), ("answer-draft.json", "draft")])
def test_missing_output_is_a_submission_error_with_runtime_provenance(tmp_path, missing, stage):
    script = "" if missing == "agent-final.json" else write_output("agent-final.json", json.dumps(HANDOFF))
    runtime, attempt = fake_runtime(tmp_path, script)
    with pytest.raises(AnswerSubmissionError) as caught:
        submit(runtime, tmp_path, attempt)
    assert caught.value.payload["stage"] == stage
    assert caught.value.usage == {"input_tokens": 7}


@pytest.mark.parametrize("stale", ["answer-draft.json", "agent-final.json", "answer-schema.json"])
def test_preexisting_outputs_are_preserved_and_never_reused(tmp_path, stale):
    runtime, attempt = fake_runtime(tmp_path, "(attempt / 'executed').touch()")
    (attempt / stale).write_text("previous output", "utf-8")
    with pytest.raises(RuntimeError, match="already contains") as caught:
        submit(runtime, tmp_path, attempt)
    assert not isinstance(caught.value, AnswerSubmissionError)
    assert not (attempt / "executed").exists()
    assert (attempt / stale).read_text("utf-8") == "previous output"


def test_attempt_cannot_be_outside_the_workspace(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    runtime, attempt = fake_runtime(tmp_path, "")
    with pytest.raises(RuntimeError, match="inside the investigation"):
        submit(runtime, workspace, attempt)
    assert not (attempt / "answer-schema.json").exists()


@pytest.mark.parametrize("name", ["answer-draft.json", "agent-final.json"])
def test_rejects_output_hardlinks_without_reading_external_content(tmp_path, name):
    outside = tmp_path / "outside.json"
    outside.write_text("External content must not be used", "utf-8")
    script = f"os.link({str(outside)!r}, attempt / {name!r})\n"
    if name == "answer-draft.json":
        script += write_output("agent-final.json", json.dumps(HANDOFF))
    runtime, attempt = fake_runtime(tmp_path, script)
    with pytest.raises(AnswerSubmissionError, match="not a link") as caught:
        submit(runtime, tmp_path, attempt)
    assert "External content" not in json.dumps(caught.value.payload)
    assert outside.read_text("utf-8") == "External content must not be used"


def test_rejects_output_symlink(tmp_path):
    outside = tmp_path / "outside.json"
    outside.write_text(json.dumps(ANSWER), "utf-8")
    probe = tmp_path / "symlink-probe"
    try:
        probe.symlink_to(outside)
    except OSError:
        pytest.skip("Symlink creation is unavailable in this test environment")
    probe.unlink()
    runtime, attempt = fake_runtime(
        tmp_path, f"(attempt / 'answer-draft.json').symlink_to({str(outside)!r})\n"
        + write_output("agent-final.json", json.dumps(HANDOFF)),
    )
    with pytest.raises(AnswerSubmissionError, match="not a link"):
        submit(runtime, tmp_path, attempt)


def test_real_process_failure_is_not_reclassified_as_answer_correction(tmp_path):
    runtime, attempt = fake_runtime(tmp_path, "sys.exit(3)")
    with pytest.raises(RuntimeError, match="exit 3") as caught:
        submit(runtime, tmp_path, attempt)
    assert not isinstance(caught.value, AnswerSubmissionError)


@pytest.mark.parametrize("encoded,exit_code", [(False, 1), (True, 1), (True, 0)])
def test_failed_turn_surfaces_provider_reason_when_stderr_is_empty(tmp_path, encoded, exit_code):
    reason = "The 'requested-model' model is not supported when using Codex with a ChatGPT account."
    provider_error = {"type": "error", "status": 400,
                      "error": {"type": "invalid_request_error", "message": reason}}
    detail = json.dumps(provider_error) if encoded else provider_error
    events = [
        {"type": "item.completed", "item": {"id": "warning", "type": "error",
                                             "message": "Model metadata not found."}},
        {"type": "error", "message": "Earlier connection error"},
        {"type": "turn.failed", "error": detail},
    ]
    # No trailing newline: the final error must still be consumed on process exit.
    runtime, attempt = fake_runtime(
        tmp_path, f"sys.stdout.write({''.join(json.dumps(e) + chr(10) for e in events).rstrip()!r})\n"
        f"sys.exit({exit_code})",
    )
    with pytest.raises(RuntimeError) as caught:
        submit(runtime, tmp_path, attempt)
    assert not isinstance(caught.value, AnswerSubmissionError)
    assert reason in str(caught.value)
    assert "--agent-model MODEL_ID" in str(caught.value)
    assert "Earlier connection error" not in str(caught.value)
    assert "Model metadata not found" not in str(caught.value)
    assert (attempt / "runtime.stderr.txt").read_text("utf-8") == ""
    observed = json.loads((attempt / "runtime-observations.json").read_text("utf-8"))
    assert observed["runtime_error"] == reason


@pytest.mark.parametrize("event,reason", [
    ({"type": "error", "message": "Authentication required"}, "Authentication required"),
    ({"type": "item.completed", "item": {"id": "failure", "type": "error",
                                          "message": "Invalid configuration"}}, "Invalid configuration"),
])
def test_exit_without_failed_turn_uses_jsonl_error(tmp_path, event, reason):
    runtime, attempt = fake_runtime(
        tmp_path, f"print({json.dumps(event)!r})\nsys.exit(1)",
    )
    with pytest.raises(RuntimeError, match=reason) as caught:
        submit(runtime, tmp_path, attempt)
    assert "--agent-model" not in str(caught.value)


def test_failure_without_jsonl_reason_keeps_stderr_diagnostic(tmp_path):
    runtime, attempt = fake_runtime(tmp_path, "sys.stderr.write('Unable to load configuration')\nsys.exit(2)")
    with pytest.raises(RuntimeError, match="Unable to load configuration"):
        submit(runtime, tmp_path, attempt)


def test_success_with_metadata_warning_remains_successful(tmp_path):
    warning = {"type": "item.completed", "item": {"id": "warning", "type": "error",
                                                "message": "Model metadata not found"}}
    runtime, attempt = fake_runtime(
        tmp_path, f"print({json.dumps(warning)!r})\n"
        + write_output("answer-draft.json", json.dumps(ANSWER))
        + write_output("agent-final.json", json.dumps(HANDOFF)),
    )
    assert submit(runtime, tmp_path, attempt).answer == ANSWER
    observed = json.loads((attempt / "runtime-observations.json").read_text("utf-8"))
    assert observed["runtime_error"] is None
