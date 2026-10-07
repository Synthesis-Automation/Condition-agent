"""Fixed-file answer transport; scientific validation remains in the service."""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path
from typing import Any

ANSWER_HANDOFF_VERSION = "scientific_answer_handoff.v1"
ANSWER_FILENAME = "answer-draft.json"
HANDOFF_FILENAME = "agent-final.json"
ANSWER_HANDOFF_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "properties": {
        "schema_version": {"type": "string", "enum": [ANSWER_HANDOFF_VERSION]},
        "answer_file": {"type": "string", "enum": [ANSWER_FILENAME]},
    },
    "required": ["schema_version", "answer_file"],
}
_EXPECTED_HANDOFF = {"schema_version": ANSWER_HANDOFF_VERSION, "answer_file": ANSWER_FILENAME}
_MAX_ANSWER_BYTES = 8 * 1024 * 1024


class AnswerSubmissionError(ValueError):
    """A completed agent submitted unusable output, eligible for one correction."""

    def __init__(
        self, message: str, *, thread_id: str, usage: dict[str, Any],
        payload: dict[str, Any],
    ) -> None:
        super().__init__(message)
        self.thread_id = thread_id
        self.usage = dict(usage)
        self.payload = payload


def _is_link(path: Path, information: os.stat_result) -> bool:
    return path.is_symlink() or bool(
        getattr(information, "st_file_attributes", 0)
        & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
    )


def _attempt_directory(workspace: Path, turn_directory: Path) -> Path:
    workspace, directory = workspace.absolute(), turn_directory.absolute()
    if not directory.is_relative_to(workspace):
        raise ValueError("Answer attempt directory must be inside the investigation")
    current = directory
    while True:
        information = current.lstat()
        if _is_link(current, information) or not stat.S_ISDIR(information.st_mode):
            raise ValueError("Answer attempt directories must be ordinary directories, not links")
        if current == workspace:
            break
        current = current.parent
    if not directory.resolve().is_relative_to(workspace.resolve()):
        raise ValueError("Answer attempt directory resolves outside the investigation")
    return directory


def prepare_answer_handoff(workspace: Path, turn_directory: Path) -> Path:
    """Require a fresh attempt directory without deleting or accepting prior output."""
    try:
        directory = _attempt_directory(workspace, turn_directory)
        stale = [name for name in (
            ANSWER_FILENAME, HANDOFF_FILENAME, "agent-final.txt", "answer-schema.json", "answer-handoff-schema.json",
            "runtime-request.json", "runtime.jsonl", "runtime.stderr.txt", "runtime-observations.json", "prompt.txt",
        )
                 if os.path.lexists(directory / name)]
        if stale:
            raise ValueError("Answer attempt already contains output; use a fresh correction attempt directory")
        return directory
    except (OSError, ValueError) as exc:
        # A preflight failure is an adapter/setup error, not a completed agent
        # submission. It must not spend an answer-correction attempt.
        raise RuntimeError(str(exc)) from exc


def answer_handoff_prompt(prompt: str, turn_directory: Path) -> str:
    """Tell the runtime where to save this attempt and return only its small receipt."""
    return prompt + "\n\nANSWER FILE HANDOFF FOR THIS ATTEMPT:\n" + (
        f"Save the complete ScientificAnswer JSON at {json.dumps(str(turn_directory / ANSWER_FILENAME))}.\n"
        f"The full answer schema is available at {json.dumps(str(turn_directory / 'answer-schema.json'))}; "
        "read it only if needed.\n"
        "Validate that exact draft with finalize_answer; a separate self-review is optional. "
        "Then finish with only this JSON object:\n"
        + json.dumps(_EXPECTED_HANDOFF)
        + "\nDo not repeat the full answer in your final message or print it through a shell. "
        "The service reads the saved draft once and validates it before publication. "
        "Each correction attempt has its own output path; use the path above, not an earlier draft.\n"
    )


def _object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON field: {key[:80]}")
        result[key] = value
    return result


def _invalid_constant(value: str) -> None:
    raise ValueError(f"Invalid JSON numeric constant: {value}")


def _read_file_bytes(workspace: Path, directory: Path, name: str, limit: int) -> bytes:
    path = directory / name
    information = path.lstat()
    if (_is_link(path, information) or not stat.S_ISREG(information.st_mode)
            or information.st_nlink != 1):
        raise ValueError(f"{name} must be an ordinary file, not a link")
    if information.st_size > limit:
        raise ValueError(f"{name} exceeds its {limit}-byte transport limit")
    with path.open("rb") as stream:
        opened = os.fstat(stream.fileno())
        if (opened.st_dev, opened.st_ino) != (information.st_dev, information.st_ino):
            raise ValueError(f"{name} changed while opening the answer")
        _attempt_directory(workspace, directory)
        raw = stream.read(limit + 1)
        after = os.fstat(stream.fileno())
    if len(raw) > limit or (after.st_size, after.st_mtime_ns) != (information.st_size, information.st_mtime_ns):
        raise ValueError(f"{name} changed while reading the answer")
    return raw


def _read_json_file(workspace: Path, directory: Path, name: str, limit: int) -> Any:
    raw = _read_file_bytes(workspace, directory, name, limit)
    # UTF-8 BOM is accepted because Windows text writers commonly emit one.
    return json.loads(raw.decode("utf-8-sig"), object_pairs_hook=_object_pairs, parse_constant=_invalid_constant)


def load_text_answer(workspace: Path, turn_directory: Path) -> dict[str, Any]:
    """Read free-form output with file integrity checks and no answer schema."""
    directory = _attempt_directory(workspace, turn_directory)
    text = _read_file_bytes(workspace, directory, "agent-final.txt", _MAX_ANSWER_BYTES).decode("utf-8-sig")
    if not text.strip():
        raise ValueError("Agent returned an empty final answer")
    return {"schema_version": "agent_text.v1", "answer_markdown": text}


def load_answer_handoff(
    workspace: Path, turn_directory: Path, *, thread_id: str, usage: dict[str, Any],
) -> dict[str, Any]:
    """Read a strict receipt and its fixed local draft once, preserving rejected files."""
    stage = "handoff"
    metadata: dict[str, Any] = {"schema_version": ANSWER_HANDOFF_VERSION}
    try:
        directory = _attempt_directory(workspace, turn_directory)
        handoff = _read_json_file(workspace, directory, HANDOFF_FILENAME, 4096)
        metadata["json_type"] = type(handoff).__name__
        if isinstance(handoff, dict):
            metadata["fields"] = list(handoff)[:12]
            for key in ("schema_version", "answer_file"):
                if isinstance(handoff.get(key), str):
                    metadata["received_" + key] = handoff[key][:200]
        if handoff != _EXPECTED_HANDOFF:
            raise ValueError("Expected the exact scientific_answer_handoff.v1 receipt for answer-draft.json")
        stage = "draft"
        answer = _read_json_file(workspace, directory, ANSWER_FILENAME, _MAX_ANSWER_BYTES)
        if not isinstance(answer, dict):
            metadata["draft_json_type"] = type(answer).__name__
            raise ValueError("answer-draft.json must contain a ScientificAnswer JSON object")
        return answer
    except (OSError, UnicodeError, ValueError) as exc:
        raise AnswerSubmissionError(str(exc), thread_id=thread_id, usage=usage,
                                    payload={**metadata, "stage": stage, "issue": type(exc).__name__}) from exc
