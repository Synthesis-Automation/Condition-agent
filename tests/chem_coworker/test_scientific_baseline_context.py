"""Scientific replay stays strict while recorded application layers can evolve."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from chem_coworker.scientific_workspace.core.baseline import (
    BASELINE_SCHEMA,
    LEGACY_BASELINE_SCHEMA,
    code_manifest,
    legacy_code_manifest,
    scientific_identity,
    verify_baseline,
)
from chem_coworker.scientific_workspace.agent_context.context import (
    capture_application_context,
    current_application_context,
    record_application_context,
)
from chem_coworker.scientific_workspace.agent_context.learning import (
    DEVELOPMENT_PARTITION,
    available_task_guides,
    build_learning_context,
    publish_lessons,
    recall_lessons,
    record_lesson,
    task_guide,
)
from chem_coworker.scientific_workspace.core.store import InvestigationStore


@pytest.fixture
def project(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    repository = tmp_path / "project"
    files = {
        "reactive_taxonomy/chemistry/example.py": "CHEMISTRY = 1\n",
        "reactive_taxonomy/definitions/example.json": '{"version": "1"}\n',
        "chem_coworker/scientific_workspace/workspace.py": "EXECUTION = 1\n",
        "chem_coworker/scientific_workspace/adapters/operations.py": "OPERATIONS = 1\n",
        "chem_coworker/scientific_workspace/core/process_utils.py": "PROCESS_CONTROL = 1\n",
        "chem_coworker/scientific_workspace/answers/answer_contracts.py": "CONTRACT = 1\n",
        "chem_coworker/scientific_workspace/runtime/conversation.py": "RUNTIME = 1\n",
        "chem_coworker/scientific_workspace/runtime/agent_runtime.py": "HARNESS = 1\n",
        "chem_coworker/scientific_workspace/agent_context/prompts.py": "PROMPT_COMPOSITION = 1\n",
        "chem_coworker/scientific_workspace/agent_instructions/core.md": "Core principles\n",
        "chem_coworker/scientific_workspace/presentation/default.md": "Brief output\n",
        "chem_coworker/scientific_workspace/task_playbooks/conditions.md": "Original guide\n",
        "app/web_api/scientific_chat.js": "const render = 1;\n",
    }
    for relative, text in files.items():
        path = repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
    monkeypatch.setattr(
        "chem_coworker.scientific_workspace.core.baseline.environment_versions",
        lambda: {"python": "fixture"},
    )
    return repository


def baseline(repository: Path, *, legacy: bool = False) -> dict:
    result = {
        "schema_version": LEGACY_BASELINE_SCHEMA if legacy else BASELINE_SCHEMA,
        "repository": str(repository),
        "code_files": (legacy_code_manifest if legacy else code_manifest)(repository),
        "environment": {"python": "fixture"}, "artifacts": {},
        "evaluation_partition": DEVELOPMENT_PARTITION,
    }
    result["learning_context"] = build_learning_context(result)
    if not legacy:
        result["scientific_identity"] = scientific_identity(result)
    return result


@pytest.mark.parametrize("build_platform,architecture", [
    ("win-amd64", "AMD64"), ("win32", "x86"),
    ("win-arm64", "ARM64"), ("win-arm32", "ARM"),
])
def test_windows_identity_does_not_depend_on_wmi_or_processor_environment(
    monkeypatch: pytest.MonkeyPatch, build_platform: str, architecture: str,
) -> None:
    from chem_coworker.scientific_workspace.core import baseline as module

    monkeypatch.setattr(module.sys, "platform", "win32")
    monkeypatch.setattr(module.sys, "getwindowsversion", lambda: SimpleNamespace(
        major=10, minor=0, build=26300,
    ), raising=False)
    monkeypatch.setattr(module.sysconfig, "get_platform", lambda: build_platform)

    def forbidden_probe() -> str:
        raise AssertionError("Windows identity must not probe WMI")

    monkeypatch.setattr(module.platform, "machine", forbidden_probe)
    monkeypatch.setattr(module.platform, "platform", forbidden_probe)
    original = module.environment_versions()
    monkeypatch.delenv("PROCESSOR_ARCHITECTURE", raising=False)
    monkeypatch.delenv("PROCESSOR_ARCHITEW6432", raising=False)
    assert module.environment_versions() == original
    assert original["platform"] == f"Windows-10.0.26300-{architecture}"


def test_runtime_mismatch_identifies_all_changed_fields(
    project: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    frozen = baseline(project)
    monkeypatch.setattr(
        "chem_coworker.scientific_workspace.core.baseline.environment_versions",
        lambda: {"python": "different", "rdkit": "new"},
    )
    with pytest.raises(ValueError, match="Scientific runtime changed") as failure:
        verify_baseline(frozen)
    assert "python: recorded='fixture', current='different'" in str(failure.value)
    assert "rdkit: recorded=None, current='new'" in str(failure.value)
    assert frozen["environment"] == {"python": "fixture"}


@pytest.mark.parametrize("relative", [
    "chem_coworker/scientific_workspace/runtime/conversation.py",
    "chem_coworker/scientific_workspace/runtime/agent_runtime.py",
    "chem_coworker/scientific_workspace/agent_context/prompts.py",
    "chem_coworker/scientific_workspace/agent_instructions/core.md",
    "chem_coworker/scientific_workspace/task_playbooks/conditions.md",
    "chem_coworker/scientific_workspace/presentation/default.md",
    "app/web_api/scientific_chat.js",
])
def test_application_changes_leave_scientific_baseline_valid(project: Path, relative: str) -> None:
    frozen = baseline(project)
    (project / relative).write_text("Updated application layer\n", encoding="utf-8")
    verify_baseline(frozen, full_hash=True)
    assert scientific_identity(frozen) == scientific_identity(baseline(project))


@pytest.mark.parametrize("relative", [
    "reactive_taxonomy/chemistry/example.py",
    "reactive_taxonomy/definitions/example.json",
    "chem_coworker/scientific_workspace/workspace.py",
    "chem_coworker/scientific_workspace/adapters/operations.py",
    "chem_coworker/scientific_workspace/core/process_utils.py",
    "chem_coworker/scientific_workspace/answers/answer_contracts.py",
])
def test_scientific_or_integrity_changes_require_new_investigation(project: Path, relative: str) -> None:
    frozen = baseline(project)
    (project / relative).write_text("Updated scientific behavior\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Scientific code or definitions changed"):
        verify_baseline(frozen)
    assert scientific_identity(frozen) != scientific_identity(baseline(project))


@pytest.mark.parametrize("relative", [
    "chem_coworker/scientific_workspace/new_adapter.py",
    "chem_coworker/scientific_workspace/runtime/new_adapter.py",
    "chem_coworker/scientific_workspace/agent_context/new_adapter.py",
    "chem_coworker/scientific_workspace/new_contract.json",
    "chem_coworker/scientific_workspace/task_playbooks/new_adapter.py",
])
def test_unclassified_new_adapters_and_contracts_are_scientific(project: Path, relative: str) -> None:
    frozen = baseline(project)
    (project / relative).write_text("New tool contract\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Scientific code or definitions changed"):
        verify_baseline(frozen)


def test_legacy_scope_remains_strict_and_is_never_silently_migrated(project: Path) -> None:
    legacy = baseline(project, legacy=True)
    frozen_v2 = baseline(project)
    guide = project / "chem_coworker/scientific_workspace/task_playbooks/conditions.md"
    guide.write_text("New task guidance\n", encoding="utf-8")
    verify_baseline(frozen_v2)
    with pytest.raises(ValueError, match="Scientific code or definitions changed"):
        verify_baseline(legacy)
    assert legacy["schema_version"] == LEGACY_BASELINE_SCHEMA
    assert legacy["learning_context"]["guides"]["conditions"]["text"] == "Original guide\n"


def test_legacy_scope_still_verifies_original_guide_directory(project: Path) -> None:
    relative = "chem_coworker/scientific_workspace/guides/conditions.md"
    path = project / relative
    path.parent.mkdir()
    path.write_text("Original directory snapshot\n", encoding="utf-8")
    frozen = baseline(project, legacy=True)
    assert relative in frozen["code_files"]
    path.write_text("Changed old resource\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Scientific code or definitions changed"):
        verify_baseline(frozen)


def test_unknown_baseline_schema_is_rejected(project: Path) -> None:
    frozen = baseline(project)
    frozen["schema_version"] = "scientific_baseline.v99"
    with pytest.raises(ValueError, match="Unsupported scientific baseline schema"):
        verify_baseline(frozen)


def test_turn_context_refresh_is_explicit_and_keeps_earlier_resources(project: Path, tmp_path: Path) -> None:
    frozen = baseline(project)
    original = deepcopy(frozen)
    store = InvestigationStore.create(tmp_path / "investigation", objective="Conditions", baseline=frozen)
    runtime = {"model": "first", "settings": {"profile": "quick"}}
    first = record_application_context(store, turn_id="first", runtime=runtime)
    recorded_first = store.read_artifact(first.artifact_ref)
    runtime["settings"]["profile"] = "mutated after recording"
    assert recorded_first["runtime"]["settings"]["profile"] == "quick"
    assert current_application_context(store)["turn_id"] == "first"
    guide = project / "chem_coworker/scientific_workspace/task_playbooks/conditions.md"
    guide.write_text("Improved guide\n", encoding="utf-8", newline="\n")
    (guide.parent / "new_task.md").write_text("New task advice\n", encoding="utf-8", newline="\n")
    profile = project / "chem_coworker/scientific_workspace/presentation/default.md"
    profile.write_text("Expanded output\n", encoding="utf-8", newline="\n")
    # Neither filesystem edits nor merely capturing the next context changes the
    # guide currently visible to the agent.
    pending = capture_application_context(store, runtime={"model": "second"})
    assert task_guide(store, "conditions")["text"] == "Original guide\n"
    assert available_task_guides(store) == ("conditions",)
    second = record_application_context(store, turn_id="second", runtime={"model": "second"})
    recorded_second = store.read_artifact(second.artifact_ref)
    assert recorded_second["resources"] == pending["resources"]
    assert task_guide(store, "conditions")["text"] == "Improved guide\n"
    assert available_task_guides(store) == ("conditions", "new_task")
    assert recorded_first["resources"]["presentation/default.md"]["text"] == "Brief output\n"
    assert recorded_second["resources"]["presentation/default.md"]["text"] == "Expanded output\n"
    assert recorded_first["scientific_identity"] == recorded_second["scientific_identity"]
    assert recorded_first["sha256"] != recorded_second["sha256"]
    assert store.manifest["baseline"] == original
    reopened = InvestigationStore(store.root)
    assert current_application_context(reopened)["turn_id"] == "second"
    assert task_guide(reopened, "conditions")["text"] == "Improved guide\n"
    assert reopened.read_artifact(first.artifact_ref) == recorded_first


def test_legacy_headless_guide_access_and_context_capture_still_work(project: Path, tmp_path: Path) -> None:
    frozen = baseline(project, legacy=True)
    store = InvestigationStore.create(tmp_path / "headless", objective="Conditions", baseline=frozen)
    assert current_application_context(store) is None
    assert task_guide(store, "conditions")["text"] == "Original guide\n"
    record_application_context(store, turn_id="first", runtime={})
    assert task_guide(store, "conditions")["text"] == "Original guide\n"


def test_baseline_identity_rejects_modified_scientific_inputs(project: Path) -> None:
    frozen = baseline(project)
    frozen["artifacts"]["condition_index"] = {"path": "invented", "status": "missing"}
    with pytest.raises(ValueError, match="Scientific baseline identity mismatch"):
        verify_baseline(frozen)


def test_turn_refresh_keeps_lessons_pinned_but_future_investigations_can_recall(
    project: Path, tmp_path: Path,
) -> None:
    store = InvestigationStore.create(
        tmp_path / "pinned", objective="Conditions", baseline=baseline(project),
    )
    first = record_application_context(store, turn_id="first", runtime={"model": "first"})
    evidence = store.append("call", {
        "operation": "fixture", "execution_status": "error",
        "error": {"message": "Unknown fixture operation"},
    })
    record_lesson(
        store, "conditions", "Inspect the catalog after an unknown operation error.",
        "An unfamiliar operation name failed.", [evidence.artifact_ref],
    )
    assert publish_lessons(store)["published"] == 1
    assert recall_lessons(store, "conditions")["lessons"] == []
    second = record_application_context(store, turn_id="second", runtime={"model": "second"})
    assert recall_lessons(store, "conditions")["lessons"] == []
    original = store.read_artifact(first.artifact_ref)
    refreshed = store.read_artifact(second.artifact_ref)
    assert original["learning_context"] == refreshed["learning_context"]
    assert original["guidance_identity"] == refreshed["guidance_identity"]
    assert original["sha256"] != refreshed["sha256"]
    future = InvestigationStore.create(
        tmp_path / "future", objective="Conditions", baseline=baseline(project),
    )
    assert len(recall_lessons(future, "conditions")["lessons"]) == 1


def test_crlf_guides_have_same_identity_at_creation_and_turn_capture(project: Path, tmp_path: Path) -> None:
    path = project / "chem_coworker/scientific_workspace/task_playbooks/conditions.md"
    path.write_bytes(b"Windows guide\r\n")
    frozen = baseline(project)
    store = InvestigationStore.create(tmp_path / "crlf", objective="Conditions", baseline=frozen)
    record_application_context(store, turn_id="first", runtime={})
    current = current_application_context(store)
    assert current["learning_context"] == frozen["learning_context"]
    assert current["resources"]["task_playbooks/conditions.md"]["text"] == "Windows guide\r\n"


def test_captured_context_is_independent_of_nested_baseline_memory(project: Path, tmp_path: Path) -> None:
    frozen = baseline(project)
    original = deepcopy(frozen)
    store = InvestigationStore.create(tmp_path / "isolated", objective="Conditions", baseline=frozen)
    captured = capture_application_context(store, runtime={})
    captured["learning_context"]["lessons"]["conditions"].append({"advice": "Caller mutation"})
    captured["learning_context"]["warnings"].append("Caller warning")
    assert store.manifest["baseline"] == original
    assert recall_lessons(store, "conditions")["lessons"] == []
    later = capture_application_context(store, runtime={})
    assert later["learning_context"] == original["learning_context"]


def test_ui_renderer_evolution_is_recorded_without_changing_science(project: Path, tmp_path: Path) -> None:
    store = InvestigationStore.create(tmp_path / "renderer", objective="Conditions", baseline=baseline(project))
    before = capture_application_context(store, runtime={})
    renderer = project / "app/web_api/scientific_chat.js"
    renderer.write_text("const render = 2;\n", encoding="utf-8")
    after = capture_application_context(store, runtime={})
    assert before["scientific_identity"] == after["scientific_identity"]
    assert before["guidance_identity"] != after["guidance_identity"]
    before_files = before["layers"]["presentation"]["code_files"]
    after_files = after["layers"]["presentation"]["code_files"]
    assert before_files["app/web_api/scientific_chat.js"] != after_files["app/web_api/scientific_chat.js"]
