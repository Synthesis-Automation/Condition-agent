"""Immutable JSON artifacts and append-only investigation events."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from time import sleep
from typing import Any, Iterator, Mapping


SCHEMA_VERSION = "scientific_investigation.v1"


def _json_default(value: Any) -> Any:
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if is_dataclass(value) and not isinstance(value, type):
        return asdict(value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Unsupported artifact type: {type(value).__name__}")


def canonical_bytes(value: Any) -> bytes:
    """Serialize scientific results without silently accepting NaN or unknown objects."""
    return json.dumps(
        value, default=_json_default, sort_keys=True, ensure_ascii=False,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")


def _write_json(path: Path, value: Any) -> None:
    data = canonical_bytes(value)
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".writing-")
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        # Windows readers briefly hold handles that can deny replacement. Retry
        # only this transient error; never expose a partially written JSON file.
        for attempt in range(8):
            try:
                os.replace(temporary, path)
                break
            except PermissionError:
                if attempt == 7:
                    raise
                sleep(0.01 * (attempt + 1))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@dataclass(frozen=True)
class InvestigationEvent:
    """One immutable call, note, or lifecycle transition."""

    sequence: int
    kind: str
    artifact_ref: str
    evidence_refs: tuple[str, ...]
    created_at: str
    schema_version: str = SCHEMA_VERSION


class InvestigationStore:
    """Local single-writer store; readers can reopen it in any agent session."""

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root).expanduser().resolve()
        self.manifest = json.loads((self.root / "investigation.json").read_text("utf-8"))
        if self.manifest.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported investigation schema")

    @classmethod
    def create(
        cls, root: str | Path, *, objective: str, baseline: Mapping[str, Any],
        constraints: tuple[str, ...] = (), agent_metadata: Mapping[str, Any] | None = None,
    ) -> "InvestigationStore":
        """Create a new directory; existing investigations are never overwritten."""
        if not objective.strip():
            raise ValueError("An investigation requires an objective")
        path = Path(root).expanduser().resolve()
        path.mkdir(parents=True, exist_ok=False)
        (path / "artifacts").mkdir()
        (path / "events").mkdir()
        _write_json(path / "investigation.json", {
            "schema_version": SCHEMA_VERSION, "objective": objective,
            "constraints": constraints, "baseline": baseline,
            "agent_metadata": dict(agent_metadata or {}),
            "created_at": datetime.now(timezone.utc).isoformat(),
        })
        return cls(path)

    @contextmanager
    def _writer(self) -> Iterator[None]:
        lock = self.root / ".writer.lock"
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as exc:
            raise RuntimeError("Investigation writer active; inspect stale lock after a crash") from exc
        try:
            os.close(descriptor)
            yield
        finally:
            lock.unlink()

    def read_artifact(self, reference: str) -> Any:
        """Read and verify a content-addressed artifact; reject path-like references."""
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", reference):
            raise ValueError("Invalid artifact reference")
        digest = reference.split(":", 1)[1]
        data = (self.root / "artifacts" / f"{digest}.json").read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"Artifact checksum mismatch: {reference}")
        return json.loads(data)

    def events(self) -> tuple[InvestigationEvent, ...]:
        """Load history in sequence order, detecting missing or invalid events."""
        result = []
        for sequence, path in enumerate(sorted((self.root / "events").glob("*.json")), 1):
            value = json.loads(path.read_text("utf-8"))
            if value["sequence"] != sequence or path.name != f"{sequence:08d}.json":
                raise ValueError("Investigation history is not contiguous")
            if value["schema_version"] != SCHEMA_VERSION:
                raise ValueError("Unsupported event schema")
            value["evidence_refs"] = tuple(value["evidence_refs"])
            self.read_artifact(value["artifact_ref"])
            for reference in value["evidence_refs"]:
                self.read_artifact(reference)
            result.append(InvestigationEvent(**value))
        return tuple(result)

    def append(
        self, kind: str, value: Any, *, evidence_refs: tuple[str, ...] = (),
    ) -> InvestigationEvent:
        """Atomically append an event whose evidence and complete payload are retrievable."""
        data = canonical_bytes(value)
        digest = hashlib.sha256(data).hexdigest()
        with self._writer():
            for reference in evidence_refs:
                self.read_artifact(reference)
            events = self.events()
            artifact_path = self.root / "artifacts" / f"{digest}.json"
            if not artifact_path.exists():
                _write_json(artifact_path, value)
            else:
                self.read_artifact(f"sha256:{digest}")
            event = InvestigationEvent(
                sequence=len(events) + 1, kind=kind, artifact_ref=f"sha256:{digest}",
                evidence_refs=evidence_refs, created_at=datetime.now(timezone.utc).isoformat(),
            )
            _write_json(self.root / "events" / f"{event.sequence:08d}.json", event)
            return event

    def note(
        self, kind: str, text: str, *, evidence_refs: tuple[str, ...] = (),
    ) -> InvestigationEvent:
        """Record agent hypotheses and decisions without labeling them observations."""
        if kind not in {"hypothesis", "decision", "question", "limitation", "review"}:
            raise ValueError("Unsupported note kind")
        if not text.strip():
            raise ValueError("Note text is required")
        return self.append(kind, {
            "text": text, "origin": "agent_authored", "review_status": "unreviewed",
        }, evidence_refs=evidence_refs)

    def attach_file(
        self, path: str | Path, *, description: str, evidence_refs: tuple[str, ...] = (),
    ) -> InvestigationEvent:
        """Preserve a UTF-8 custom script or report as a derived artifact, without executing it."""
        source = Path(path).resolve()
        data = source.read_bytes()
        return self.append("derived_file", {
            "source_path": str(source), "description": description,
            "content": data.decode("utf-8"), "source_sha256": hashlib.sha256(data).hexdigest(),
            "origin": "derived_analysis", "review_status": "unreviewed",
        }, evidence_refs=evidence_refs)

    def set_status(self, status: str, reason: str) -> InvestigationEvent:
        """Append a lifecycle transition, retaining prior completion and resumption history."""
        if status not in {"active", "completed", "insufficient_evidence", "failed", "cancelled", "budget_exhausted"}:
            raise ValueError("Unsupported investigation status")
        if not reason.strip():
            raise ValueError("A status transition requires a reason")
        return self.append("status", {"status": status, "reason": reason})

    def summary(self) -> dict[str, Any]:
        """Return enough context for another session to resume without the original chat."""
        events = self.events()
        states = [event for event in events if event.kind == "status"]
        return {
            "objective": self.manifest["objective"], "constraints": self.manifest["constraints"],
            "status": self.read_artifact(states[-1].artifact_ref)["status"] if states else "active",
            "baseline_validation": self.manifest["baseline"].get("validation_status"),
            "events": [asdict(event) for event in events],
            "notes": [{"kind": event.kind, "artifact_ref": event.artifact_ref,
                       **self.read_artifact(event.artifact_ref)}
                      for event in events if event.kind in {"hypothesis", "decision", "question", "limitation", "review"}],
        }
