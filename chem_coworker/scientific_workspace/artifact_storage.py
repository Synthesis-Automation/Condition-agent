"""Lossless, versioned storage projection for large route-assessment evidence.

Scientific contracts stay unchanged. Only explicitly named diagnostic fields in
route results are separated; structures, gates, warnings and precedents stay inline.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
from typing import Any, Callable


STORAGE_VERSION = "scientific_artifact_storage.v1"
SECTION_VERSION = "scientific_evidence_section.v1"
REFERENCE_VERSION = "scientific_evidence_reference.v1"
MIN_SECTION_BYTES = 2048
_ROUTE_OPERATIONS = frozenset({
    "assess_route_step", "assess_route_proposal", "revise_route_branch", "inspect_route_step",
})
_RESULT_SCHEMAS = frozenset({
    "route_step_investigation.v1", "route_investigation.v1", "route_step_inspection.v1",
})


def compact_route_artifact(kind: str, payload: Any) -> tuple[Any, dict[str, dict[str, Any]]]:
    """Return a compact storage envelope and deduplicated content-addressed sections.

    Small results and all other event kinds preserve their original representation.
    No input is mutated, and this function does not run scientific calculations.
    """
    from .store import canonical_bytes

    if not isinstance(payload, dict):
        return payload, {}
    result = payload.get("result")
    if not isinstance(result, dict) or result.get("schema_version") not in _RESULT_SCHEMAS:
        return payload, {}
    if kind != "replay" and not (kind == "call" and payload.get("execution_status") == "completed"
                                  and payload.get("operation") in _ROUTE_OPERATIONS):
        return payload, {}
    projected = deepcopy(payload)
    sections: dict[str, dict[str, Any]] = {}
    locations = []

    def separate(parent: dict[str, Any], key: str, path: list[str | int]) -> None:
        value = parent[key]
        if not isinstance(value, (dict, list)) or len(canonical_bytes(value)) < MIN_SECTION_BYTES:
            return
        section = {"schema_version": SECTION_VERSION, "value": value}
        reference = "sha256:" + hashlib.sha256(canonical_bytes(section)).hexdigest()
        sections[reference] = section
        preview = {name: value[name] for name in (
            "signature_id", "schema_version", "evidence_quality", "transformation_class",
        ) if isinstance(value, dict) and name in value}
        parent[key] = {"schema_version": REFERENCE_VERSION, "artifact_ref": reference,
                       "value_type": "object" if isinstance(value, dict) else "array",
                       "item_count": len(value), "summary": preview}
        locations.append(path + [key])

    def walk(value: Any, path: list[str | int]) -> None:
        if isinstance(value, list):
            for index, item in enumerate(value):
                walk(item, path + [index])
        elif isinstance(value, dict):
            for name, item in list(value.items()):
                if name == "template_ids" and len(path) >= 2 and path[-2] == "operator_matches":
                    separate(value, name, path)
                elif name == "reaction_signature" and "operator_matches" in value:
                    separate(value, name, path)
                elif name == "molecule_audits" and isinstance(item, dict):
                    for smiles in item:
                        separate(item, smiles, path + [name])
                else:
                    walk(item, path + [name])

    walk(projected["result"], ["result"])
    if not sections:
        return payload, {}
    return {
        "storage_schema_version": STORAGE_VERSION,
        "expanded_sha256": hashlib.sha256(canonical_bytes(payload)).hexdigest(),
        "section_paths": locations,
        "payload": projected,
    }, sections


def expand_route_artifact(stored: Any, read_section: Callable[[str], Any]) -> Any:
    """Restore the exact original payload, validating every linked section and hash.

    Section contents are terminal values, never recursively executed or expanded.
    That keeps cycles and arbitrary path resolution out of the storage protocol.
    """
    from .store import canonical_bytes

    if not isinstance(stored, dict) or "storage_schema_version" not in stored:
        return stored
    if stored["storage_schema_version"] != STORAGE_VERSION:
        raise ValueError("Unsupported scientific artifact storage schema")
    if not isinstance(stored.get("payload"), dict) or not isinstance(stored.get("section_paths"), list):
        raise ValueError("Malformed compact scientific artifact")
    payload = deepcopy(stored["payload"])
    seen = set()
    sections = {}
    for path in stored["section_paths"]:
        if (not isinstance(path, list) or not path or path[0] != "result"
                or not all(type(key) in (str, int) for key in path)):
            raise ValueError("Invalid compact evidence path")
        key = tuple(path)
        if key in seen:
            raise ValueError("Duplicate compact evidence path")
        seen.add(key)
        try:
            parent = payload
            for part in path[:-1]:
                parent = parent[part]
            marker = parent[path[-1]]
            if not isinstance(marker, dict) or marker.get("schema_version") != REFERENCE_VERSION:
                raise ValueError("Invalid compact evidence reference")
            reference = marker["artifact_ref"]
            if reference not in sections:
                sections[reference] = read_section(reference)
            section = sections[reference]
            if not isinstance(section, dict) or section.get("schema_version") != SECTION_VERSION:
                raise ValueError("Invalid scientific evidence section")
            value = section["value"]
            expected = dict if marker["value_type"] == "object" else list if marker["value_type"] == "array" else None
            if expected is None or not isinstance(value, expected) or len(value) != marker["item_count"]:
                raise ValueError("Scientific evidence section type or count mismatch")
            parent[path[-1]] = deepcopy(value)
        except (KeyError, IndexError, TypeError) as exc:
            raise ValueError("Malformed compact evidence link") from exc
    if hashlib.sha256(canonical_bytes(payload)).hexdigest() != stored.get("expanded_sha256"):
        raise ValueError("Expanded scientific artifact checksum mismatch")
    return payload
