"""Lossless JSON shard storage with shared chemistry and evidence objects.

Each shard is self contained. References can only address objects emitted earlier
in that shard, which prevents cycles and keeps streaming readers deterministic.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, TextIO

STORAGE_SCHEMA_VERSION = "reaction_object_shard.v1"
_REFERENCE_KEY = "$reaction_object"


def json_bytes(value: Any) -> bytes:
    """Serialize evidence deterministically without lossy transformations."""
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":")).encode("utf-8")


def hydrate(value: Any, objects: Mapping[str, Any]) -> Any:
    """Restore a JSON value, rejecting dangling or malformed references."""
    if isinstance(value, dict):
        if set(value) == {_REFERENCE_KEY}:
            identity = value[_REFERENCE_KEY]
            if not isinstance(identity, str) or identity not in objects:
                raise ValueError(f"Missing reaction object: {identity}")
            return hydrate(objects[identity], objects)
        return {key: hydrate(item, objects) for key, item in value.items()}
    if isinstance(value, list):
        return [hydrate(item, objects) for item in value]
    return value


def write_object_records(handle: TextIO, records: Iterable[Mapping[str, Any]]) -> int:
    """Write records and deduplicate large nested values within a bounded shard."""
    handle.write(json_bytes({"storage_schema": STORAGE_SCHEMA_VERSION}).decode() + "\n")
    objects: dict[str, bytes] = {}

    def intern(value: Any, *, root: bool = False) -> Any:
        if isinstance(value, dict):
            if _REFERENCE_KEY in value:
                raise ValueError("Source field conflicts with reserved object reference")
            compact = {str(key): intern(item) for key, item in value.items()}
        elif isinstance(value, (list, tuple)):
            compact = [intern(item) for item in value]
        else:
            return value
        encoded = json_bytes(compact)
        if root or len(encoded) < 512:
            return compact
        identity = hashlib.sha256(encoded).hexdigest()
        previous = objects.get(identity)
        if previous is None:
            objects[identity] = encoded
            handle.write(json_bytes({"object_id": identity, "value": compact}).decode() + "\n")
        elif previous != encoded:
            raise ValueError("Reaction object hash collision")
        return {_REFERENCE_KEY: identity}

    count = 0
    for record in records:
        compact = intern(dict(record), root=True)
        handle.write(json_bytes({"record": compact}).decode() + "\n")
        count += 1
    return count


def iter_record_events(path: str | Path) -> Iterator[dict[str, Any]]:
    """Read storage events or wrap ordinary JSONL records as record events."""
    source = Path(path)
    opener = gzip.open if source.suffix == ".gz" else open
    object_format = False
    with opener(source, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"Invalid JSONL at line {number}: {error.msg}") from error
            if not isinstance(value, dict):
                raise ValueError(f"Not a JSON object at {source}:{number}")
            if set(value) == {"storage_schema"}:
                if value["storage_schema"] != STORAGE_SCHEMA_VERSION:
                    raise ValueError("Unsupported reaction object storage schema")
                object_format = True
                yield {"reset_objects": True}
                continue
            if not object_format:
                yield {"record": value}
            elif set(value) == {"object_id", "value"}:
                if hashlib.sha256(json_bytes(value["value"])).hexdigest() != value["object_id"]:
                    raise ValueError(f"Reaction object checksum failed at {source}:{number}")
                yield value
            elif set(value) == {"record"} and isinstance(value["record"], dict):
                yield value
            else:
                raise ValueError(f"Invalid storage event at {source}:{number}")


def iter_record_shard(path: str | Path) -> Iterator[dict[str, Any]]:
    """Hydrate canonical records from an ordinary or object-normalized shard."""
    objects: dict[str, Any] = {}
    for event in iter_record_events(path):
        if "reset_objects" in event:
            objects.clear()
        elif "object_id" in event:
            # References must already resolve before admitting an object.
            value = event["value"]
            _check_references(value, objects)
            objects[event["object_id"]] = value
        else:
            yield hydrate(event["record"], objects)


def _check_references(value: Any, objects: Mapping[str, Any]) -> None:
    if isinstance(value, dict):
        if set(value) == {_REFERENCE_KEY}:
            if value[_REFERENCE_KEY] not in objects:
                raise ValueError("Forward or missing reaction object reference")
        else:
            for item in value.values():
                _check_references(item, objects)
    elif isinstance(value, list):
        for item in value:
            _check_references(item, objects)
