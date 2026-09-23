"""Recorded execution and reproducible replay for external scientific agents."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from time import monotonic
from typing import Any, Mapping

from .baseline import capture_baseline, verify_baseline
from .operations import ScientificOperations
from .store import canonical_bytes, InvestigationEvent, InvestigationStore


class ScientificWorkspace:
    """A resumable scientific workspace usable directly from Python or the CLI."""

    def __init__(self, root: str | Path) -> None:
        self.store = InvestigationStore(root)
        self.operations = ScientificOperations(self.store)

    @classmethod
    def create(
        cls, root: str | Path, *, objective: str, repository: str | Path,
        artifacts: Mapping[str, str | Path] | None = None,
        constraints: tuple[str, ...] = (), agent_metadata: Mapping[str, Any] | None = None,
    ) -> "ScientificWorkspace":
        """Identify selected inputs and initialize an unvalidated development investigation."""
        repository = Path(repository).resolve()
        paths = {name: (repository / Path(path)).resolve() for name, path in (artifacts or {}).items()}
        baseline = capture_baseline(repository, paths)
        InvestigationStore.create(
            root, objective=objective, baseline=baseline,
            constraints=constraints, agent_metadata=agent_metadata,
        )
        return cls(root)

    def run(
        self, operation: str, arguments: Mapping[str, Any], *, evidence_refs: tuple[str, ...] = (),
    ) -> InvestigationEvent:
        """Record full inputs and results, including failures, without remapping domain statuses."""
        # Snapshot caller-owned input containers before invoking mutable third-party code.
        inputs = json.loads(canonical_bytes(arguments))
        started = monotonic()
        payload: dict[str, Any] = {
            "operation": operation, "arguments": inputs,
            "origin": "deterministic_computation", "review_status": "unreviewed",
        }
        try:
            for reference in evidence_refs:
                self.store.read_artifact(reference)
            if self.store.summary()["status"] != "active":
                raise ValueError("Investigation is stopped; append an active status to resume")
            verify_baseline(self.store.manifest["baseline"])
            payload["result"] = json.loads(canonical_bytes(
                self.operations.invoke(operation, json.loads(canonical_bytes(inputs)))
            ))
            payload["execution_status"] = "completed"
        except KeyboardInterrupt:
            payload.update(execution_status="cancelled", error={"type": "KeyboardInterrupt", "message": "Execution interrupted"})
        except Exception as exc:
            payload.update(execution_status="error", error={"type": type(exc).__name__, "message": str(exc)})
        payload["duration_seconds"] = round(monotonic() - started, 6)
        if operation == "revise_routes" and isinstance(inputs.get("source_ref"), str):
            try:
                self.store.read_artifact(inputs["source_ref"])
            except (OSError, ValueError):
                pass
            else:
                evidence_refs = tuple(dict.fromkeys((*evidence_refs, inputs["source_ref"])))
        # Invalid caller references remain in the saved request/error, not in verified links.
        valid_refs = []
        for reference in evidence_refs:
            try:
                self.store.read_artifact(reference)
            except (OSError, ValueError):
                continue
            valid_refs.append(reference)
        event = self.store.append("call", payload, evidence_refs=tuple(valid_refs))
        if payload["execution_status"] == "cancelled":
            self.store.set_status("cancelled", "Scientific operation interrupted")
        return event

    def replay(self, reference: str) -> InvestigationEvent:
        """Recompute a saved successful call against full baseline hashes and compare results."""
        source = self.store.read_artifact(reference)
        if source.get("execution_status") != "completed" or "operation" not in source:
            raise ValueError("Only completed scientific calls can be replayed")
        verify_baseline(self.store.manifest["baseline"], full_hash=True)
        actual = self.operations.invoke(source["operation"], source["arguments"])
        return self.store.append("replay", {
            "source_ref": reference,
            "matches": canonical_bytes(actual) == canonical_bytes(source["result"]),
            "result": actual,
        }, evidence_refs=(reference,))

    def call_summary(self, event: InvestigationEvent) -> dict[str, Any]:
        """Return a compact pointer while retaining complete scientific output on disk."""
        value = self.store.read_artifact(event.artifact_ref)
        result = value.get("result")
        summary = {"event": asdict(event), "operation": value.get("operation"),
                   "execution_status": value.get("execution_status"), "error": value.get("error")}
        if isinstance(result, dict):
            summary["result_summary"] = {key: result[key] for key in (
                "valid", "error", "status", "retrieval_level", "candidate_count", "warnings", "total",
            ) if key in result}
        return summary
