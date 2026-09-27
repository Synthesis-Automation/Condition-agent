"""Optional advice is pinned, evidence-linked, and isolated from scientific truth."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.baseline import code_manifest
from chem_coworker.scientific_workspace.learning import (
    DEVELOPMENT_PARTITION, build_learning_context,
)


@pytest.fixture
def investigations(tmp_path, monkeypatch):
    repository = tmp_path / "project"
    source = repository / "chem_coworker" / "example.py"
    source.parent.mkdir(parents=True)
    source.write_text("VERSION = 1\n", encoding="utf-8")
    environment = {"python": "fixture-1"}
    monkeypatch.setattr("chem_coworker.scientific_workspace.baseline.environment_versions", lambda: environment.copy())
    stores = []

    def create(*, partition=DEVELOPMENT_PARTITION, repo=repository):
        baseline = {
            "repository": str(repo), "code_files": code_manifest(repo),
            "environment": environment.copy(), "artifacts": {}, "evaluation_partition": partition,
        }
        baseline["learning_context"] = build_learning_context(baseline)
        directory = tmp_path / f"investigation-{len(stores)}"
        InvestigationStore.create(directory, objective="Test procedural memory", baseline=baseline)
        workspace = ScientificWorkspace(directory)
        stores.append(workspace)
        return workspace

    create.repository = repository
    create.source = source
    create.environment = environment
    create.lesson_path = repository / "results" / "ai_native" / "lessons.jsonl"
    return create


def evidence(workspace, *, failed=True):
    return workspace.store.append("call", {
        "operation": "fixture_operation", "execution_status": "error" if failed else "completed",
        "error": {"message": "Fixture network access denied"} if failed else None,
    }).artifact_ref


def lesson(workspace, reference=None, **overrides):
    arguments = {
        "task": "general", "advice": "Inspect captured browser text after a denied download.",
        "applies_when": "Direct download is denied and a browser can access the source.",
        "evidence_refs": [reference or evidence(workspace)],
    }
    arguments.update(overrides)
    return workspace.record_lesson(**arguments)


def test_advice_is_published_once_and_only_recalled_by_future_investigations(investigations):
    current = investigations()
    event = lesson(current)
    stored = current.store.read_artifact(event.artifact_ref)
    assert stored["authority"] == "optional_procedural_advice"
    assert stored["review_status"] == "agent_authored_unreviewed"
    assert not investigations.lesson_path.exists()
    assert current.recall_lessons("retrosynthesis")["lessons"] == []
    assert current.publish_lessons()["published"] == 1
    assert current.publish_lessons()["published"] == 0
    assert current.recall_lessons("retrosynthesis")["lessons"] == []
    future = investigations()
    recalled = future.recall_lessons("retrosynthesis")
    assert [row["lesson_id"] for row in recalled["lessons"]] == [stored["lesson_id"]]
    assert recalled["authority"] == "optional_advice_not_scientific_evidence"
    assert recalled["lessons"][0]["source_event_ref"] == event.artifact_ref
    reopened = ScientificWorkspace(current.store.root)
    assert reopened.recall_lessons("retrosynthesis")["lessons"] == []


def test_zero_lessons_creates_no_store_or_lock(investigations):
    workspace = investigations()
    assert workspace.publish_lessons()["published"] == 0
    assert not investigations.lesson_path.parent.exists()


@pytest.mark.parametrize("overrides", [
    {"task": "new-task"}, {"scope": "global"}, {"advice": " "},
    {"advice": "x" * 1201}, {"applies_when": ""}, {"evidence_refs": []},
    {"evidence_refs": ["../../other.json"]},
])
def test_invalid_lesson_is_not_recorded(investigations, overrides):
    workspace = investigations()
    with pytest.raises(ValueError):
        lesson(workspace, **overrides)
    assert not any(event.kind == "lesson" for event in workspace.store.events())


def test_missing_or_asserted_evidence_is_rejected(investigations):
    workspace = investigations()
    note = workspace.store.append("note", {"text": "I think this succeeded"})
    with pytest.raises(ValueError, match="recorded evidence"):
        lesson(workspace, reference=note.artifact_ref)
    with pytest.raises(FileNotFoundError):
        lesson(workspace, reference="sha256:" + "0" * 64)
    first = lesson(workspace)
    with pytest.raises(ValueError, match="recorded evidence"):
        lesson(workspace, reference=first.artifact_ref)


def test_three_lessons_per_turn_then_new_turn_can_record(investigations):
    workspace = investigations()
    reference = evidence(workspace, failed=False)
    for index in range(3):
        lesson(workspace, reference=reference, advice=f"Procedural lesson {index}")
    with pytest.raises(ValueError, match="at most three"):
        lesson(workspace, reference=reference)
    workspace.store.append("user_message", {"text": "Investigate another question"})
    lesson(workspace, reference=reference)
    assert workspace.publish_lessons()["published"] == 4


def test_retrieval_prefers_exact_task_deduplicates_and_bounds_context(investigations):
    first = investigations()
    lesson(first, advice="General recovery")
    lesson(first, advice="Condition-specific recovery", task="conditions")
    lesson(first, advice="Condition-specific recovery", task="conditions")
    first.publish_lessons()
    second = investigations()
    lesson(second, advice="New general recovery")
    lesson(second, advice="New condition recovery", task="conditions")
    second.publish_lessons()
    future = investigations()
    assert [row["advice"] for row in future.recall_lessons("conditions")["lessons"]] == [
        "New condition recovery", "Condition-specific recovery", "New general recovery",
    ]
    assert [row["advice"] for row in future.recall_lessons("retrosynthesis")["lessons"]] == [
        "New general recovery", "General recovery",
    ]
    assert len(future.recall_lessons("conditions", limit=1)["lessons"]) == 1
    with pytest.raises(ValueError):
        future.recall_lessons("conditions", limit=4)


def test_version_scopes_and_repository_isolation(investigations, tmp_path):
    current = investigations()
    lesson(current, advice="Code-specific recovery")
    lesson(current, advice="Environment-specific recovery", scope="environment")
    current.publish_lessons()
    investigations.source.write_text("VERSION = 2\n", encoding="utf-8")
    future = investigations()
    assert [row["advice"] for row in future.recall_lessons("general")["lessons"]] == ["Environment-specific recovery"]
    investigations.environment["python"] = "fixture-2"
    assert investigations().recall_lessons("general")["lessons"] == []
    investigations.environment["python"] = "fixture-1"
    foreign = investigations(repo=tmp_path / "other-project")
    assert foreign.recall_lessons("general")["lessons"] == []
    # Copying a project's advice into another project cannot confer local provenance.
    context = build_learning_context(foreign.store.manifest["baseline"], lesson_path=investigations.lesson_path)
    assert context["lessons"]["general"] == []
    assert any("another project" in warning for warning in context["warnings"])


@pytest.mark.parametrize("partition", ["untouched_evaluation", "blind_review", None])
def test_non_development_runs_neither_learn_nor_recall(investigations, partition):
    development = investigations()
    lesson(development)
    development.publish_lessons()
    evaluation = investigations(partition=partition)
    assert evaluation.recall_lessons("conditions")["enabled"] is False
    assert evaluation.recall_lessons("conditions")["lessons"] == []
    assert evaluation.task_guide("conditions")["text"]
    with pytest.raises(ValueError, match="only for development"):
        lesson(evaluation)
    assert evaluation.publish_lessons() == {"published": 0, "enabled": False}


def test_retirement_is_evidence_linked_append_only_and_future_only(investigations):
    origin = investigations()
    event = lesson(origin)
    identity = origin.store.read_artifact(event.artifact_ref)["lesson_id"]
    origin.publish_lessons()
    current = investigations()
    current.retire_lesson(identity, "Later evidence contradicts this procedure", [evidence(current)])
    assert current.publish_lessons()["published"] == 1
    assert current.recall_lessons("general")["lessons"]
    assert investigations().recall_lessons("general")["lessons"] == []
    records = [json.loads(line) for line in investigations.lesson_path.read_text("utf-8").splitlines()]
    assert [record["status"] for record in records] == ["active", "retired"]
    with pytest.raises(ValueError, match="frozen context"):
        current.retire_lesson("lesson:unknown", "Does not exist", [evidence(current)])


@pytest.mark.parametrize("damage", ["changed_advice", "missing_evidence", "malformed_store", "invalid_shape"])
def test_damaged_memory_is_skipped_with_warning_not_promoted(investigations, damage):
    origin = investigations()
    reference = evidence(origin)
    lesson(origin, reference=reference)
    origin.publish_lessons()
    if damage in {"changed_advice", "invalid_shape"}:
        record = json.loads(investigations.lesson_path.read_text("utf-8"))
        record["advice"] = "Unrecorded claim" if damage == "changed_advice" else None
        investigations.lesson_path.write_text(json.dumps(record) + "\n", encoding="utf-8")
    elif damage == "missing_evidence":
        (origin.store.root / "artifacts" / (reference.removeprefix("sha256:") + ".json")).unlink()
    else:
        investigations.lesson_path.write_text('{"unfinished":', encoding="utf-8")
    recalled = investigations().recall_lessons("general")
    assert recalled["lessons"] == []
    assert recalled["warnings"]


def test_unverified_retirement_cannot_suppress_recorded_advice(investigations):
    origin = investigations()
    lesson(origin)
    origin.publish_lessons()
    original = json.loads(investigations.lesson_path.read_text("utf-8"))
    forged = {key: original[key] for key in (
        "schema_version", "lesson_id", "source_run", "evidence_refs",
        "evaluation_partition", "source_event_ref",
    )}
    forged.update(status="retired", reason="Unrecorded retirement")
    with investigations.lesson_path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(forged) + "\n")
    recalled = investigations().recall_lessons("general")
    assert [row["lesson_id"] for row in recalled["lessons"]] == [original["lesson_id"]]
    assert recalled["warnings"]


def test_guides_are_hashed_and_context_corruption_is_detected(investigations):
    workspace = investigations()
    guide = workspace.task_guide("retrosynthesis")
    assert "optional menu" in guide["text"]
    assert "disconnect_target" in guide["text"]
    assert "multistep planner" in " ".join(guide["text"].split())
    path = investigations.repository / "chem_coworker" / "scientific_workspace" / "guides" / "conditions.md"
    path.parent.mkdir(parents=True)
    path.write_text("First guide", encoding="utf-8")
    before = code_manifest(investigations.repository)
    path.write_text("Updated guide", encoding="utf-8")
    assert before != code_manifest(investigations.repository)
    corrupted = deepcopy(workspace.store.manifest["baseline"]["learning_context"])
    corrupted["guides"]["retrosynthesis"]["text"] = "Unrecorded guide"
    workspace.store.manifest["baseline"]["learning_context"] = corrupted
    with pytest.raises(ValueError, match="checksum"):
        workspace.task_guide("retrosynthesis")


def test_old_investigation_remains_usable_without_learning_context(investigations):
    workspace = investigations()
    del workspace.store.manifest["baseline"]["learning_context"]
    with pytest.raises(ValueError, match="no frozen guidance"):
        workspace.recall_lessons("general")
