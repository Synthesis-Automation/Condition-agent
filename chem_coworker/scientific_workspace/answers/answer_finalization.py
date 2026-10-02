"""Small answer-authoring helper over the existing scientific answer contract."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping

from ..core.store import InvestigationStore, _write_json
from .answer_contracts import ScientificAnswer, validate_answer_evidence
from .answer_handoff import ANSWER_FILENAME, ANSWER_HANDOFF_VERSION, _attempt_directory
from .evidence_review import record_evidence_review


def _complete_empty_fields(draft: Mapping[str, Any]) -> dict[str, Any]:
    """Fill authoring boilerplate only; never infer chemistry, attribution or support."""
    value = deepcopy(dict(draft))
    value.setdefault("schema_version", "scientific_answer.v2")
    value.setdefault("needs_user_input", False)
    for name in ("evidence_refs", "uncertainties", "sources", "molecules",
                 "target_molecule_ids", "steps", "routes", "claims"):
        value.setdefault(name, [])

    def defaults(item: Any, **fields: Any) -> None:
        if isinstance(item, dict):
            for key, default in fields.items():
                item.setdefault(key, default)

    def entries(container: Any, name: str) -> list[Any]:
        items = container.get(name) if isinstance(container, dict) else None
        return items if isinstance(items, list) else []

    for source in entries(value, "sources"):
        defaults(source, url=None)
    for name in ("molecules", "claims", "steps"):
        for item in entries(value, name):
            defaults(item, source_ids=[], limitations=[])
    for step in entries(value, "steps"):
        defaults(step, after_step_ids=[], conditions=[], yield_info=None)
        for condition in entries(step, "conditions"):
            defaults(condition, source_ids=[], limitations=[])
        if isinstance(step, dict):
            defaults(step.get("yield_info"), source_ids=[], limitations=[])
            defaults(step.get("rationale"), source_ids=[], limitations=[])
    for route in entries(value, "routes"):
        defaults(route, limitations=[])
    return value


def finalize_answer(
    store: InvestigationStore, draft_path: str | Path, draft: Mapping[str, Any], *,
    findings: list[dict[str, Any]] | None = None,
) -> dict[str, str]:
    """Validate, optionally record explicit self-review, and atomically save a draft.

    Missing empty collections/nulls may be omitted while authoring. Basis, IDs,
    source locators, scientific text and review findings remain the author's
    responsibility. The stored file is the unchanged scientific_answer.v2 contract;
    the service still validates it before publication. No scientific tools run here.
    """
    path = Path(draft_path)
    if not path.is_absolute():
        path = store.root / path
    if path.name != ANSWER_FILENAME:
        raise ValueError(f"Use the current attempt's {ANSWER_FILENAME} path")
    directory = _attempt_directory(store.root, path.parent)
    answer = ScientificAnswer.model_validate(_complete_empty_fields(draft))
    answer.evidence_refs = validate_answer_evidence(answer, store)
    payload = answer.model_dump()
    if findings is not None:
        record_evidence_review(store, payload, findings)
    _write_json(directory / ANSWER_FILENAME, payload)
    return {"schema_version": ANSWER_HANDOFF_VERSION, "answer_file": ANSWER_FILENAME}
