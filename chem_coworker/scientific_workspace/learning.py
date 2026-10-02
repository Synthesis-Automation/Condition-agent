"""Small, evidence-linked advisory memory; never a chemistry rule registry."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
from typing import Any, Iterator, Literal, Mapping

from .store import InvestigationEvent, InvestigationStore, canonical_bytes


Task = Literal["general", "conditions", "retrosynthesis"]
Scope = Literal["code", "environment"]
TASKS = ("general", "conditions", "retrosynthesis")
DEVELOPMENT_PARTITION = "development_only_not_an_untouched_evaluation"
LESSON_SCHEMA = "scientific_lesson.v1"
_EVIDENCE_KINDS = {"call", "replay", "custom_execution", "derived_file",
                   "literature_source", "literature_excerpt", "agent_turn_error"}


def _digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _scope_versions(baseline: Mapping[str, Any], scope: str) -> dict[str, Any]:
    versions = {"environment": baseline["environment"]}
    if scope == "code":
        versions["code_sha256"] = _digest(baseline["code_files"])
    return versions


def _bounded_text(value: str, name: str, limit: int) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ValueError(f"{name} must be a nonempty string of at most {limit} characters")
    return value.strip()


def _validate_record(value: Any, *, published: bool = True) -> None:
    """Reject malformed optional records with a recoverable validation error."""
    if not isinstance(value, dict) or value.get("schema_version") != LESSON_SCHEMA:
        raise ValueError("Unsupported lesson store record")
    if not isinstance(value.get("lesson_id"), str) or not re.fullmatch(r"lesson:[0-9a-f]{64}", value["lesson_id"]):
        raise ValueError("Invalid lesson identity")
    _bounded_text(value.get("source_run"), "source_run", 4096)
    refs = value.get("evidence_refs")
    if not isinstance(refs, list) or not 1 <= len(refs) <= 10 or not all(
        isinstance(ref, str) and re.fullmatch(r"sha256:[0-9a-f]{64}", ref) for ref in refs
    ):
        raise ValueError("Invalid lesson evidence references")
    if value.get("evaluation_partition") != DEVELOPMENT_PARTITION:
        raise ValueError("Only development lessons are eligible")
    if value.get("status") == "active":
        if value.get("task") not in TASKS or value.get("scope") not in {"code", "environment"}:
            raise ValueError("Invalid lesson task or scope")
        _bounded_text(value.get("advice"), "advice", 1200)
        _bounded_text(value.get("applies_when"), "applies_when", 800)
        if not isinstance(value.get("versions"), dict):
            raise ValueError("Missing lesson versions")
        if (value.get("authority") != "optional_procedural_advice"
                or value.get("review_status") != "agent_authored_unreviewed"):
            raise ValueError("Lessons must remain unreviewed procedural advice")
    elif value.get("status") == "retired":
        _bounded_text(value.get("reason"), "reason", 800)
    else:
        raise ValueError("Unknown lesson status")
    if published and (not isinstance(value.get("source_event_ref"), str)
                      or not re.fullmatch(r"sha256:[0-9a-f]{64}", value["source_event_ref"])):
        raise ValueError("Missing originating lesson event")


def _records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    result = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            value = json.loads(line)
            _validate_record(value)
            result.append(value)
    return result


def _active(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    active: dict[str, dict[str, Any]] = {}
    for record in records:
        if record.get("status") == "active":
            active[record["lesson_id"]] = record
        elif record.get("status") == "retired":
            active.pop(record["lesson_id"], None)
        else:
            raise ValueError("Unknown lesson status")
    return list(active.values())


def _verify_published_record(record: Mapping[str, Any], baseline: Mapping[str, Any]) -> None:
    source = InvestigationStore(record["source_run"])
    origin = source.manifest["baseline"]
    if (Path(origin["repository"]).resolve() != Path(baseline["repository"]).resolve()
            or origin.get("evaluation_partition") != DEVELOPMENT_PARTITION):
        raise ValueError("Lesson source belongs to another project or evaluation partition")
    original = source.read_artifact(record["source_event_ref"])
    saved = {key: value for key, value in record.items() if key not in {"source_event_ref", "published_at"}}
    if canonical_bytes(original) != canonical_bytes(saved):
        raise ValueError("Published advice differs from its recorded lesson")
    expected_kind = "lesson" if record["status"] == "active" else "lesson_retirement"
    if not any(event.kind == expected_kind and event.artifact_ref == record["source_event_ref"] for event in source.events()):
        raise ValueError("Lesson must reference its originating lesson event")
    if record["status"] == "active" and record["versions"] != _scope_versions(origin, record["scope"]):
        raise ValueError("Lesson versions differ from the originating baseline")
    _evidence(source, record["evidence_refs"])


def build_learning_context(
    baseline: Mapping[str, Any], *, lesson_path: Path | None = None,
    guide_snapshots: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    """Freeze available guides and up to three applicable advisory lessons per task."""
    if guide_snapshots is None:
        directory = Path(baseline["repository"]) / "chem_coworker" / "scientific_workspace" / "guides"
        if not directory.is_dir():
            directory = Path(__file__).with_name("guides")
        guides = {}
        for path in sorted(directory.glob("*.md")):
            data = path.read_bytes()
            guides[path.stem] = {"task": path.stem,
                                 "sha256": hashlib.sha256(data).hexdigest(),
                                 "text": data.decode("utf-8")}
    else:
        guides = {task: dict(guide) for task, guide in guide_snapshots.items()}
    path = lesson_path or Path(baseline["repository"]) / "results" / "ai_native" / "lessons.jsonl"
    enabled = baseline.get("evaluation_partition") == DEVELOPMENT_PARTITION
    selected: dict[str, list[dict[str, Any]]] = {task: [] for task in TASKS}
    warnings = []
    if enabled:
        try:
            verified = []
            for record in _records(path):
                try:
                    _verify_published_record(record, baseline)
                except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
                    warnings.append(f"Skipped unverifiable lesson {record.get('lesson_id')}: {exc}")
                    continue
                verified.append(record)
            eligible = [record for record in _active(verified)
                        if record["evaluation_partition"] == DEVELOPMENT_PARTITION
                        and record["scope"] in {"code", "environment"}
                        and record["versions"] == _scope_versions(baseline, record["scope"])]
            for task in TASKS:
                seen = set()
                # Exact task before general advice; newest within each category.
                candidates = [r for r in reversed(eligible) if r["task"] == task]
                if task != "general":
                    candidates += [r for r in reversed(eligible) if r["task"] == "general"]
                for record in candidates:
                    identity = (record["advice"].casefold(), record["applies_when"].casefold())
                    if identity not in seen:
                        selected[task].append(record)
                        seen.add(identity)
                    if len(selected[task]) == 3:
                        break
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            selected = {task: [] for task in TASKS}
            warnings.append(f"Advisory memory unavailable: {type(exc).__name__}: {exc}")
    payload = {"schema_version": "scientific_learning_context.v1", "enabled": enabled,
               "lesson_store": str(path.resolve()), "guides": guides, "lessons": selected,
               "warnings": warnings, "authority": "optional_advice_not_scientific_evidence"}
    return {**payload, "sha256": _digest(payload)}


def _context(store: InvestigationStore) -> dict[str, Any]:
    from .context import current_application_context

    application = current_application_context(store)
    context = (application["learning_context"] if application is not None
               else store.manifest["baseline"].get("learning_context"))
    if context is None:
        raise ValueError("This investigation has no frozen guidance; start a new investigation")
    payload = {key: value for key, value in context.items() if key != "sha256"}
    if _digest(payload) != context.get("sha256"):
        raise ValueError("Learning context checksum mismatch")
    return context


def task_guide(store: InvestigationStore, task: str) -> dict[str, Any]:
    """Read a frozen optional guide; reading it imposes no required tool sequence."""
    guides = _context(store)["guides"]
    if task not in guides:
        raise ValueError(f"Unknown task guide; available: {', '.join(sorted(guides))}")
    return guides[task]


def available_task_guides(store: InvestigationStore) -> tuple[str, ...]:
    """List frozen guide names without putting task-specific names in orchestration."""
    return tuple(sorted(_context(store)["guides"]))


def recall_lessons(store: InvestigationStore, task: Task, limit: int = 3) -> dict[str, Any]:
    """Return pinned advisory lessons, never newly published lessons mid-run."""
    if task not in TASKS or type(limit) is not int or not 1 <= limit <= 3:
        raise ValueError("Use task general/conditions/retrosynthesis and limit 1..3")
    context = _context(store)
    return {"context_sha256": context["sha256"], "task": task,
            "lessons": context["lessons"][task][:limit], "warnings": context["warnings"],
            "enabled": context["enabled"], "authority": context["authority"]}


def _evidence(store: InvestigationStore, refs: list[str]) -> tuple[str, ...]:
    if not isinstance(refs, list) or not 1 <= len(refs) <= 10 or not all(isinstance(ref, str) for ref in refs):
        raise ValueError("Supply 1..10 evidence_refs from this investigation")
    kinds = {event.artifact_ref: event.kind for event in store.events()}
    for ref in refs:
        store.read_artifact(ref)
        if kinds.get(ref) not in _EVIDENCE_KINDS:
            raise ValueError("Lessons require recorded evidence, not another lesson or agent assertion")
    return tuple(dict.fromkeys(refs))


def record_lesson(
    store: InvestigationStore, task: Task, advice: str, applies_when: str,
    evidence_refs: list[str], scope: Scope = "code",
) -> InvestigationEvent:
    """Record one unreviewed procedural lesson; publication occurs after the turn."""
    from .baseline import verify_baseline

    context = _context(store)
    baseline = store.manifest["baseline"]
    if not context["enabled"] or baseline.get("evaluation_partition") != DEVELOPMENT_PARTITION:
        raise ValueError("Lesson extraction is enabled only for development investigations")
    verify_baseline(baseline)
    if task not in TASKS or scope not in {"code", "environment"}:
        raise ValueError("Invalid lesson task or scope")
    refs = _evidence(store, evidence_refs)
    events = store.events()
    boundary = max((event.sequence for event in events if event.kind == "user_message"), default=0)
    if sum(event.kind == "lesson" and event.sequence > boundary for event in events) >= 3:
        raise ValueError("Record at most three lessons per turn; zero is fine")
    lesson = {"schema_version": LESSON_SCHEMA, "task": task,
              "advice": _bounded_text(advice, "advice", 1200),
              "applies_when": _bounded_text(applies_when, "applies_when", 800),
              "evidence_refs": list(refs), "source_run": str(store.root), "scope": scope,
              "versions": _scope_versions(baseline, scope), "status": "active",
              "evaluation_partition": DEVELOPMENT_PARTITION,
              "authority": "optional_procedural_advice", "review_status": "agent_authored_unreviewed"}
    lesson["lesson_id"] = "lesson:" + _digest(lesson)
    return store.append("lesson", lesson, evidence_refs=refs)


def retire_lesson(store: InvestigationStore, lesson_id: str, reason: str, evidence_refs: list[str]) -> InvestigationEvent:
    """Queue retirement of a lesson pinned in this run, preserving its prior records."""
    from .baseline import verify_baseline

    context = _context(store)
    if not context["enabled"]:
        raise ValueError("Lesson changes are enabled only for development investigations")
    verify_baseline(store.manifest["baseline"])
    known = {r["lesson_id"] for rows in context["lessons"].values() for r in rows}
    if lesson_id not in known:
        raise ValueError("Select a lesson from this investigation's frozen context")
    refs = _evidence(store, evidence_refs)
    return store.append("lesson_retirement", {
        "schema_version": LESSON_SCHEMA, "lesson_id": lesson_id, "status": "retired",
        "reason": _bounded_text(reason, "reason", 800), "source_run": str(store.root),
        "evidence_refs": list(refs), "evaluation_partition": DEVELOPMENT_PARTITION,
    }, evidence_refs=refs)


@contextmanager
def _writer(path: Path) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.with_suffix(".lock")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError("Lesson writer active; retry after it finishes or inspect a stale lock") from exc
    try:
        os.close(descriptor)
        yield
    finally:
        lock.unlink()


def publish_lessons(store: InvestigationStore) -> dict[str, Any]:
    """Idempotently publish recorded lessons from a finished worker, including failed runs."""
    context = _context(store)
    if not context["enabled"]:
        return {"published": 0, "enabled": False}
    candidates = [event for event in store.events() if event.kind in {"lesson", "lesson_retirement"}]
    if not candidates:
        return {"published": 0, "enabled": True}
    path = Path(context["lesson_store"])
    published = 0
    with _writer(path):
        previous = _records(path)
        seen = {(record.get("source_run"), record.get("source_event_ref")) for record in previous}
        with path.open("a", encoding="utf-8") as stream:
            for event in candidates:
                if (str(store.root), event.artifact_ref) in seen:
                    continue
                lesson = store.read_artifact(event.artifact_ref)
                _validate_record(lesson, published=False)
                expected_kind = "lesson" if lesson["status"] == "active" else "lesson_retirement"
                if event.kind != expected_kind or lesson["source_run"] != str(store.root):
                    raise ValueError("Lesson publication must match its originating event")
                # Verify the supporting snapshots still exist unchanged.
                _evidence(store, lesson["evidence_refs"])
                record = {**lesson, "source_event_ref": event.artifact_ref,
                          "published_at": datetime.now(timezone.utc).isoformat()}
                stream.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                stream.flush()
                seen.add((str(store.root), event.artifact_ref))
                published += 1
    return {"published": published, "enabled": True, "lesson_store": str(path)}
