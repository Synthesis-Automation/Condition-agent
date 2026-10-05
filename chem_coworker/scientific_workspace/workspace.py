"""Recorded execution and reproducible replay for external scientific agents."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from threading import Event
from time import monotonic
from typing import Any, Mapping

from .adapters.operations import ScientificOperations
from .core.baseline import capture_baseline, verify_baseline
from .core.operation_contracts import OperationDefinition, OperationProvider
from .core.store import InvestigationEvent, InvestigationStore, canonical_bytes


class ScientificWorkspace:
    """A resumable scientific workspace usable directly from Python or the CLI."""

    def __init__(self, root: str | Path, *, operations: OperationProvider | None = None) -> None:
        self.store = InvestigationStore(root)
        self.operations: OperationProvider = operations if operations is not None else ScientificOperations(self.store)

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
        timings: dict[str, float] = {}
        definition: OperationDefinition | None = None
        try:
            definition = self.operations.definition(operation)
            payload["operation_contract_version"] = definition.contract_version
            payload["scientific_identity"] = self.store.manifest["baseline"].get("scientific_identity")
            phase_started = monotonic()
            for reference in evidence_refs:
                self.store.read_artifact(reference)
            if self.store.summary()["status"] != "active":
                raise ValueError("Investigation is stopped; append an active status to resume")
            verify_baseline(self.store.manifest["baseline"])
            timings["baseline_and_evidence_seconds"] = round(monotonic() - phase_started, 6)
            phase_started = monotonic()
            result = self.operations.invoke(operation, json.loads(canonical_bytes(inputs)))
            timings["operation_seconds"] = round(monotonic() - phase_started, 6)
            phase_started = monotonic()
            serialized = canonical_bytes(result)
            payload["result"] = json.loads(serialized)
            payload["result_bytes"] = len(serialized)
            timings["serialization_seconds"] = round(monotonic() - phase_started, 6)
            payload["execution_status"] = "completed"
            if definition.execution_status_field:
                operation_status = (
                    payload["result"].get(definition.execution_status_field)
                    if isinstance(payload["result"], dict) else None
                )
                if operation_status not in ("completed", "timed_out", "error", "cancelled"):
                    raise ValueError(
                        f"Operation returned invalid {definition.execution_status_field}; "
                        "expected completed, error, timed_out or cancelled"
                    )
                if operation_status in {"timed_out", "error", "cancelled"}:
                    payload["execution_status"] = operation_status
                    payload["error"] = payload["result"].get("error", {
                        "type": "OperationIncomplete", "message": "Scientific operation did not finish",
                    })
        except KeyboardInterrupt:
            payload.update(execution_status="cancelled", error={"type": "KeyboardInterrupt", "message": "Execution interrupted"})
        except Exception as exc:
            payload.update(execution_status="error", error={"type": type(exc).__name__, "message": str(exc)})
        payload["duration_seconds"] = round(monotonic() - started, 6)
        payload["timings"] = timings
        if definition is not None:
            evidence_refs = tuple(dict.fromkeys((*evidence_refs, *definition.evidence_references(inputs))))
            if definition.result_evidence_field and isinstance(payload.get("result"), dict):
                result_refs = payload["result"].get(definition.result_evidence_field, [])
                if isinstance(result_refs, list) and all(isinstance(ref, str) for ref in result_refs):
                    evidence_refs = tuple(dict.fromkeys((*evidence_refs, *result_refs)))
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
        definition = self.operations.definition(source["operation"])
        saved_version = source.get("operation_contract_version")
        if saved_version is not None and saved_version != definition.contract_version:
            raise ValueError("Operation contract changed; replay requires the recorded contract version")
        actual = self.operations.invoke(source["operation"], source["arguments"])
        return self.store.append("replay", {
            "source_ref": reference,
            "operation_contract_version": definition.contract_version,
            "scientific_identity": self.store.manifest["baseline"].get("scientific_identity"),
            "matches": canonical_bytes(definition.replay_projection(actual))
            == canonical_bytes(definition.replay_projection(source["result"])),
            "comparison_scope": definition.replay_comparison,
            "result": actual,
        }, evidence_refs=(reference,))

    def run_python(
        self, script: str, parameters: dict[str, Any], *,
        evidence_refs: tuple[str, ...] = (), timeout_seconds: int = 60, cancel: Event | None = None,
    ) -> InvestigationEvent:
        """Record a local script's inputs/output and execution; no implicit chemistry validation."""
        from .core.execution import run_python

        return run_python(self.store, script, parameters, evidence_refs, timeout_seconds, cancel)

    def call_summary(self, event: InvestigationEvent | str, *, detailed: bool = False) -> dict[str, Any]:
        """Summarize a new event or saved reference without rerunning science."""
        from .views.brief_summaries import bound_call_summary, summarize_call_brief
        from .views.call_summaries import summarize_call

        if isinstance(event, str):
            reference = event
            event = next((item for item in self.store.events() if item.artifact_ref == reference), None)
            if event is None:
                raise ValueError("Summary requires a recorded artifact reference")
        value = self.store.read_artifact(event.artifact_ref)
        project = summarize_call if detailed else summarize_call_brief
        return bound_call_summary({"event": asdict(event), **project(value)}, value)

    def run_summary(
        self, operation: str, arguments: Mapping[str, Any], *, evidence_refs: tuple[str, ...] = (),
    ) -> dict[str, Any]:
        """Record complete scientific evidence and return only its bounded decision view."""
        return self.call_summary(self.run(operation, arguments, evidence_refs=evidence_refs))

    def describe(self) -> dict[str, Any]:
        """Show baseline identity and configured inputs without dumping code manifests."""
        baseline = self.store.manifest["baseline"]
        events = self.store.events()
        states = [event for event in events if event.kind == "status"]
        return {
            "schema_version": "scientific_workspace_description.v1",
            "objective": self.store.manifest["objective"],
            "constraints": self.store.manifest["constraints"],
            "scientific_identity": baseline.get("scientific_identity"),
            "baseline_validation": baseline.get("validation_status"),
            "repository": baseline.get("repository"),
            "artifacts": {name: {key: value.get(key) for key in ("path", "status", "sha256")}
                          for name, value in baseline.get("artifacts", {}).items()},
            "status": self.store.read_artifact(states[-1].artifact_ref)["status"] if states else "active",
            "event_count": len(events),
            "recent_events": [asdict(event) for event in events[-5:]],
            "history_omitted": max(0, len(events) - 5),
        }

    def help(self, names: str | list[str]) -> list[dict[str, Any]]:
        """Return signatures/examples for 1..10 named public helpers or scientific operations."""
        import inspect
        from .views.helper_help import HELPER_EXAMPLES, nested_input_help

        selected = [names] if isinstance(names, str) else names
        if not isinstance(selected, list) or not 1 <= len(selected) <= 10:
            raise ValueError("Select 1..10 public names")
        operations = {entry["name"]: entry for entry in self.operations.catalog()}
        entries = []
        for name in selected:
            if name in HELPER_EXAMPLES:
                method = getattr(self, name)
                entries.append({"name": name, "signature": str(inspect.signature(method)),
                                "description": inspect.getdoc(method), "example": HELPER_EXAMPLES[name]})
            elif name in operations:
                entries.append(operations[name])
            else:
                raise ValueError(f"Unknown public helper or operation: {name}")
            entries[-1].update(nested_input_help(name))
        return entries

    def batch_summary(self, references: list[str]) -> dict[str, Any]:
        """Read 1..30 saved events under one shared 16 KiB console budget; never rerun calls."""
        if not isinstance(references, list) or not 1 <= len(references) <= 30:
            raise ValueError("Supply 1..30 saved references")
        summaries = [self.call_summary(ref) for ref in references]
        items = [{"artifact_ref": item["event"]["artifact_ref"], "operation": item.get("operation"),
                  "execution_status": item.get("execution_status"), "summary_omitted": True}
                 for item in summaries]
        result = {"schema_version": "scientific_batch_summary.v1", "items": items,
                  "serialized_byte_limit": 16000,
                  "hint": "Inspect omitted calls by reference with call_summary or inspect_artifact; full evidence remains saved."}
        for index, summary in enumerate(summaries):
            original = items[index]
            items[index] = summary
            if len(canonical_bytes(result)) > 16000:
                items[index] = original
        return result

    def answer_template(self, answer_markdown: str) -> dict[str, Any]:
        """Return empty answer boilerplate; supply all chemistry, claims and attribution yourself."""
        from .answers.answer_finalization import _complete_empty_fields

        return _complete_empty_fields({"answer_markdown": answer_markdown})

    def attach_literature_reaction(
        self, draft: Mapping[str, Any], step_id: str, preparation_ref: str,
    ) -> dict[str, Any]:
        """Copy a draft and attach the immutable source drawing plus its exact required citation."""
        from .answers.authoring import attach_literature_reaction

        return attach_literature_reaction(self.store, draft, step_id, preparation_ref)

    def attach_recipe_check(
        self, draft: Mapping[str, Any], step_id: str, reference: str,
    ) -> dict[str, Any]:
        """Copy a draft and bind a completed proposed-recipe check to the exact step graphs."""
        from .answers.authoring import attach_recipe_check

        return attach_recipe_check(self.store, draft, step_id, reference)

    def answer_preflight(self, draft: Mapping[str, Any], *, offset: int = 0, limit: int = 10) -> dict[str, Any]:
        """Validate without rerunning science; page recipe/input warnings via next_offset, limit 1..10."""
        from .answers.authoring import answer_preflight

        return answer_preflight(self.store, draft, offset, limit)

    def inspect_artifact(
        self, artifact_ref: str, path: tuple[str | int, ...] | list[str | int] = (), *,
        offset: int = 0, limit: int = 5,
    ) -> dict[str, Any]:
        """Inspect a bounded saved JSON field without recomputation or mutation.

        Use literal keys/indices, for example ("result", "recommendations", 0).
        Lists/mappings support pages of 1..20 items; omitted text, nested values,
        and ancestor context are explicitly disclosed. Full artifacts stay on disk.
        """
        from .views.artifact_inspection import inspect_artifact_payload

        return inspect_artifact_payload(
            self.store.read_artifact(artifact_ref), artifact_ref=artifact_ref,
            path=path, offset=offset, limit=limit,
        )

    def capabilities(self) -> dict[str, Any]:
        """Probe local imports and configured file presence; do not assert web access."""
        from .core.capabilities import local_capabilities

        return local_capabilities(self.store.manifest["baseline"])

    def fetch_source(
        self, url: str, *, title: str | None = None, retry_network: bool = False,
        reference_id: str | None = None,
    ) -> InvestigationEvent:
        """Save a bounded public source snapshot with extraction and retrieval provenance."""
        from .adapters.literature import fetch_source

        return fetch_source(self.store, url, title=title, retry_network=retry_network, reference_id=reference_id)

    def capture_source(
        self, text: str, *, url: str, title: str | None = None, locator: str | None = None,
        reference_id: str | None = None,
    ) -> InvestigationEvent:
        """Save an agent-supplied passage without claiming independent retrieval."""
        from .adapters.literature import capture_source

        return capture_source(self.store, text, url=url, title=title, locator=locator, reference_id=reference_id)

    def capture_source_file(
        self, file: str, *, url: str, title: str | None = None, locator: str | None = None,
        reference_id: str | None = None, text_path: list[str | int] | None = None,
    ) -> InvestigationEvent:
        """Import UTF-8 text/JSON inside the investigation with exact byte/path provenance.

        JSON strings unwrap automatically; JSON objects need text_path=['text']
        or their actual literal path. reference_id is optional REF1:<64 lowercase
        hex> from the corpus, never a DOI. Put a DOI in url/title or omit reference_id.
        Captures remain agent-supplied, not independently verified.
        """
        from .adapters.literature import capture_source_file

        return capture_source_file(self.store, file, url=url, title=title, locator=locator,
                                   reference_id=reference_id, text_path=text_path)

    def capture_source_image(self, file: str, *, source_ref: str, locator: str) -> InvestigationEvent:
        """Save a PNG/JPEG source page or scheme inside the investigation, linked to captured text.

        Maximum 2 MiB/20 million pixels. This preserves bytes and page attribution;
        it does not verify chemical assignments. Pass its artifact_ref as scheme_refs
        to prepare_literature_reaction.
        """
        from .adapters.source_images import capture_source_image

        return capture_source_image(self.store, file, source_ref, locator)

    def prepared_literature_reaction(self, reference: str) -> dict[str, Any]:
        """Read an answer-ready block from a completed recorded preparation; never recompute."""
        from .adapters.literature_reactions import prepared_literature_reaction

        return prepared_literature_reaction(self.store, reference)

    def inspect_source(
        self, source_ref: str, *, query: str | None = None, offset: int = 0, limit: int = 4000,
    ) -> dict[str, Any]:
        """Read a bounded source passage with exact snapshot locations."""
        from .adapters.literature import inspect_source

        return inspect_source(self.store, source_ref, query=query, offset=offset, limit=limit)

    def record_source_excerpt(
        self, source_ref: str, *, start: int | None = None, end: int | None = None,
        excerpt: str | None = None, locator: str | None = None,
    ) -> InvestigationEvent:
        """Save an exact passage from an existing snapshot, not a generated quotation."""
        from .adapters.literature import record_source_excerpt

        return record_source_excerpt(
            self.store, source_ref, start=start, end=end, excerpt=excerpt, locator=locator,
        )

    def record_evidence_review(
        self, draft: Mapping[str, Any], findings: list[dict[str, Any]],
    ) -> InvestigationEvent:
        """Record the agent's explicit challenge of a draft; this is not independent review."""
        from .answers.evidence_review import record_evidence_review

        return record_evidence_review(self.store, draft, findings)

    def finalize_answer(
        self, draft_path: str | Path, draft: Mapping[str, Any], *,
        findings: list[dict[str, Any]] | None = None,
    ) -> dict[str, str]:
        """Save a validated answer and optional explicit self-review; return its handoff.

        Empty lists, null yields/URLs and schema boilerplate may be omitted from
        the input. Scientific basis, citations and content must remain explicit.
        """
        from .answers.answer_finalization import finalize_answer

        return finalize_answer(self.store, draft_path, draft, findings=findings)

    def task_guide(self, task: str) -> dict[str, Any]:
        """Read an available optional guide from recorded application context."""
        from .agent_context.learning import task_guide

        return task_guide(self.store, task)

    def recall_lessons(self, task: str, limit: int = 3) -> dict[str, Any]:
        """Read at most three frozen procedural lessons; these are not scientific evidence."""
        from .agent_context.learning import recall_lessons

        return recall_lessons(self.store, task, limit)

    def record_lesson(
        self, task: str, advice: str, applies_when: str, evidence_refs: list[str], scope: str = "code",
    ) -> InvestigationEvent:
        """Save evidence-linked procedural advice for later runs; at most three per turn."""
        from .agent_context.learning import record_lesson

        return record_lesson(self.store, task, advice, applies_when, evidence_refs, scope)

    def retire_lesson(self, lesson_id: str, reason: str, evidence_refs: list[str]) -> InvestigationEvent:
        """Record a correction that retires a recalled lesson for subsequent investigations."""
        from .agent_context.learning import retire_lesson

        return retire_lesson(self.store, lesson_id, reason, evidence_refs)

    def publish_lessons(self) -> dict[str, Any]:
        """Publish pending lessons after a CLI investigation; the conversation service does this itself."""
        from .agent_context.learning import publish_lessons

        return publish_lessons(self.store)
