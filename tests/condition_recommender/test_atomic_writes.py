import json
from pathlib import Path

import pytest

from condition_recommender.conversion import atomic as atomic_module


def test_atomic_json_retries_transient_permission_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "shard_manifest.json"
    destination.write_text('{"old": true}\n', encoding="utf-8")
    real_replace = atomic_module.os.replace
    attempts = 0

    def transient_replace(source: Path, target: Path) -> None:
        nonlocal attempts
        attempts += 1
        if attempts < 3:
            raise PermissionError(5, "Access is denied", str(target))
        real_replace(source, target)

    monkeypatch.setattr(atomic_module.os, "replace", transient_replace)
    monkeypatch.setattr(atomic_module.time, "sleep", lambda _delay: None)

    atomic_module.atomic_json(destination, {"new": True})

    assert json.loads(destination.read_text(encoding="utf-8")) == {"new": True}
    assert attempts == 3
    assert not tuple(tmp_path.glob(".shard_manifest.json.*.tmp"))


def test_atomic_json_preserves_destination_after_persistent_permission_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "shard_manifest.json"
    original = '{"old": true}\n'
    destination.write_text(original, encoding="utf-8")

    def denied_replace(_source: Path, target: Path) -> None:
        raise PermissionError(5, "Access is denied", str(target))

    monkeypatch.setattr(atomic_module.os, "replace", denied_replace)
    monkeypatch.setattr(atomic_module.time, "sleep", lambda _delay: None)

    with pytest.raises(PermissionError):
        atomic_module.atomic_json(destination, {"new": True})

    assert destination.read_text(encoding="utf-8") == original
    assert not tuple(tmp_path.glob(".shard_manifest.json.*.tmp"))
