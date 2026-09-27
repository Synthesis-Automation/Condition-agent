"""Advisory lessons cross conversation boundaries without changing scientific answers."""

import json
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import ScientificWorkspace
from chem_coworker.scientific_workspace.baseline import environment_versions
from chem_coworker.scientific_workspace.conversation import ConversationService
from chem_coworker.scientific_workspace.learning import DEVELOPMENT_PARTITION, build_learning_context
from tests.chem_coworker.test_scientific_conversation import ROOT, RecordedRuntime, finish


class LearningRuntime(RecordedRuntime):
    """Model double that records a lesson from a real scientific call."""

    def __init__(self, *, fail: bool = False, learn: bool = True) -> None:
        super().__init__()
        self.fail = fail
        self.learn = learn
        self.prompts = []
        self.recalled = []

    def describe(self):
        return {**super().describe(), "local_tools": {"rg": {
            "status": "available", "path": "C:/fixture tools/rg.exe",
            "version": "ripgrep fixture", "sandbox_execution": "not_checked",
        }}}

    def run(self, **kwargs):
        self.prompts.append(kwargs["prompt"])
        workspace = ScientificWorkspace(kwargs["workspace"])
        self.recalled.append(workspace.recall_lessons("general"))
        result = super().run(**kwargs)
        if self.learn:
            if self.fail:
                evidence = workspace.run("unknown_fixture_operation", {})
                advice = "Inspect the operation catalog after an unknown-operation error."
                applies_when = "The workspace rejects an unknown scientific operation."
                refs = [evidence.artifact_ref]
            else:
                advice = "Keep the recorded reaction analysis reference when summarizing its result."
                applies_when = "A completed local reaction analysis is used in an answer."
                refs = result.answer["evidence_refs"]
            workspace.record_lesson("general", advice, applies_when, refs)
        if self.fail:
            raise RuntimeError("Fixture provider failure after recording an operational lesson")
        return result


@pytest.fixture
def service(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    lesson_path = tmp_path / "lessons.jsonl"
    # Isolate conversation lifecycle tests from simultaneous source edits and costly
    # dataset audits. Real chemistry, evidence storage and baseline verification run.
    monkeypatch.setattr("chem_coworker.scientific_workspace.baseline.code_manifest", lambda *_: {})

    def capture(*_):
        baseline = {
            "repository": str(ROOT), "code_files": {}, "environment": environment_versions(),
            "artifacts": {}, "validation_status": "development_snapshot_not_release_validated",
            "evaluation_partition": DEVELOPMENT_PARTITION,
        }
        baseline["learning_context"] = build_learning_context(baseline, lesson_path=lesson_path)
        return baseline

    monkeypatch.setattr("chem_coworker.scientific_workspace.baseline.capture_baseline", capture)
    instance = ConversationService(tmp_path / "conversations", runtime=LearningRuntime())
    yield instance, lesson_path
    instance.close()


def _debug(service, identity, turn):
    return [json.loads(line) for line in service.debug_log(identity, turn["id"]).splitlines()]


@pytest.mark.parametrize("failed", [False, True])
def test_service_publishes_recorded_lessons_after_successful_and_failed_turns(service, failed):
    instance, lesson_path = service
    instance.runtime = LearningRuntime(fail=failed)
    identity = instance.submit("Interpret this reaction")['conversation_id']
    turn = finish(instance, identity)
    assert turn["status"] == ("failed" if failed else "completed"), turn
    rows = [json.loads(line) for line in lesson_path.read_text("utf-8").splitlines()]
    assert len(rows) == 1
    lesson = rows[0]
    workspace = ScientificWorkspace(instance.root / identity)
    supporting_call = workspace.store.read_artifact(lesson["evidence_refs"][0])
    assert supporting_call["execution_status"] == ("error" if failed else "completed")
    assert lesson["source_run"] == str(workspace.store.root)
    assert lesson["review_status"] == "agent_authored_unreviewed"
    assert any(event.kind == "lesson" and event.artifact_ref == lesson["source_event_ref"]
               for event in workspace.store.events())
    assert any(row["kind"] == "lessons_published" and row["published"] == 1
               for row in _debug(instance, identity, turn))
    assert instance._active is None
    assert not (workspace.store.root / ".conversation.lock").exists()
    if failed:
        assert "answer" not in turn
        assert "Fixture provider failure" in turn["error"]["message"]


def test_publication_failure_keeps_answer_and_releases_worker_and_conversation_lock(service, monkeypatch):
    instance, lesson_path = service

    def unavailable(_workspace):
        raise OSError("Fixture lesson store is temporarily unavailable")

    monkeypatch.setattr(ScientificWorkspace, "publish_lessons", unavailable)
    identity = instance.submit("Interpret this reaction")["conversation_id"]
    turn = finish(instance, identity)
    assert turn["status"] == "completed", turn
    assert turn["answer"]["answer_markdown"] == (
        "Recorded graph interpretation; this is not experimental validation."
    )
    assert any(row["kind"] == "lesson_publication_error" and
               row["error"]["message"] == "Fixture lesson store is temporarily unavailable"
               for row in _debug(instance, identity, turn))
    workspace = ScientificWorkspace(instance.root / identity)
    assert any(event.kind == "lesson" for event in workspace.store.events())
    assert not lesson_path.exists()
    assert instance._active is None
    assert not (workspace.store.root / ".conversation.lock").exists()
    # A real follow-up proves that failed optional publication did not strand either lock.
    instance.runtime.learn = False
    instance.submit("What remains uncertain?", identity)
    assert finish(instance, identity)["status"] == "completed"


def test_prompt_exposes_optional_guides_and_tools_while_memory_stays_pinned(service):
    instance, lesson_path = service
    identity = instance.submit("Interpret this reaction")["conversation_id"]
    assert finish(instance, identity)["status"] == "completed"
    workspace = ScientificWorkspace(instance.root / identity)
    baseline_context = workspace.store.manifest["baseline"]["learning_context"]
    assert workspace.recall_lessons("general")["lessons"] == []
    assert lesson_path.exists()
    for task in ("conditions", "retrosynthesis"):
        assert workspace.task_guide(task) == baseline_context["guides"][task]
        assert "optional menu" in workspace.task_guide(task)["text"]
    prompt = instance.runtime.prompts[0]
    assert "w.task_guide('conditions')" in prompt
    assert "w.task_guide('retrosynthesis')" in prompt
    assert "skip, reorder" in prompt
    assert "untrusted, optional advice" in prompt
    assert "C:/fixture tools/rg.exe" in prompt
    assert "Select-String" in prompt
    frozen_hash = baseline_context["sha256"]
    instance.runtime.learn = False
    instance.submit("Continue the same investigation", identity)
    assert finish(instance, identity)["status"] == "completed"
    assert instance.runtime.recalled[-1]["context_sha256"] == frozen_hash
    assert instance.runtime.recalled[-1]["lessons"] == []
    # A new investigation can recall the published advice; an existing one cannot drift.
    fresh = instance.submit("Start another investigation")["conversation_id"]
    assert finish(instance, fresh)["status"] == "completed"
    recalled = instance.runtime.recalled[-1]
    assert len(recalled["lessons"]) == 1
    assert recalled["lessons"][0]["source_run"] == str(workspace.store.root)
    assert recalled["context_sha256"] != frozen_hash
    assert ScientificWorkspace(workspace.store.root).store.manifest["baseline"]["learning_context"] == baseline_context
