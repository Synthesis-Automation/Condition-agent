"""Persistent checksum reuse preserves hashes, invalidation and full audits."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import pytest

from chem_coworker.scientific_workspace.core import baseline
from chem_coworker.scientific_workspace.core.fingerprint_cache import FingerprintCache


def test_unchanged_file_reuses_hash_across_cache_instances(tmp_path: Path, monkeypatch) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    cache = FingerprintCache(tmp_path / "cache")
    original = baseline.artifact_identity(path, cache=cache)
    monkeypatch.setattr(baseline, "sha256_file", lambda *a, **k: pytest.fail("Unexpected full read"))
    reused = []
    assert baseline.artifact_identity(
        path, cache=FingerprintCache(cache.directory), on_reuse=lambda: reused.append(True),
    ) == original
    assert reused == [True]


@pytest.mark.parametrize("change", ["size", "mtime", "replacement", "missing", "corrupt"])
def test_changed_files_and_invalid_cache_entries_are_not_reused(tmp_path: Path, change: str) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    cache = FingerprintCache(tmp_path / "cache")
    baseline.artifact_identity(path, cache=cache)
    original_stat = path.stat()
    if change == "size":
        path.write_bytes(b"a different length")
    elif change == "mtime":
        path.write_bytes(b"modified")
        os.utime(path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns + 1_000_000_000))
    elif change == "replacement":
        replacement = tmp_path / "replacement"
        replacement.write_bytes(b"replaced")
        os.utime(replacement, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
        replacement.replace(path)
    elif change == "missing":
        path.unlink()
    else:
        next(cache.directory.glob("*.json")).write_text("{broken", encoding="utf-8")
    reads = []
    result = baseline.artifact_identity(path, cache=cache, on_progress=reads.append)
    if change == "missing":
        assert result["status"] == "missing"
    else:
        assert reads[-1] == path.stat().st_size
        assert result["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_cancelled_or_changing_read_does_not_populate_cache(tmp_path: Path) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    cache = FingerprintCache(tmp_path / "cache")

    def cancel(count: int) -> None:
        if count:
            raise InterruptedError("cancel")

    with pytest.raises(InterruptedError):
        baseline.artifact_identity(path, cache=cache, on_progress=cancel)
    assert cache.get(path.resolve(), path.stat()) is None

    def change(count: int) -> None:
        if count:
            stat = path.stat()
            os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))

    with pytest.raises(ValueError, match="changed during fingerprinting"):
        baseline.artifact_identity(path, cache=cache, on_progress=change)
    assert cache.get(path.resolve(), path.stat()) is None


def test_expected_release_checksum_is_enforced_on_cached_files(tmp_path: Path) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    cache = FingerprintCache(tmp_path / "cache")
    baseline.artifact_identity(path, cache=cache)
    with pytest.raises(ValueError, match="checksum mismatch"):
        baseline.artifact_identity(path, cache=cache, expected_sha256="0" * 64)


def test_uncached_verification_reads_content_even_with_a_cached_hash(tmp_path: Path, monkeypatch) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    cache = FingerprintCache(tmp_path / "cache")
    baseline.artifact_identity(path, cache=cache)
    reads = []
    original_hash = baseline.sha256_file

    def tracked_hash(path, **kwargs):
        reads.append(path)
        return original_hash(path, **kwargs)

    monkeypatch.setattr(baseline, "sha256_file", tracked_hash)
    baseline.artifact_identity(path)
    assert reads == [path.resolve()]


def test_cache_unavailable_falls_back_to_full_hash(tmp_path: Path) -> None:
    path = tmp_path / "data"
    path.write_bytes(b"original")
    blocked = tmp_path / "not_a_directory"
    blocked.write_text("file", encoding="utf-8")
    assert baseline.artifact_identity(path, cache=FingerprintCache(blocked)) == baseline.artifact_identity(path)


def test_warm_capture_preserves_scientific_identity(tmp_path: Path) -> None:
    repository = Path(__file__).resolve().parents[2]
    path = tmp_path / "source.bin"
    path.write_bytes(b"source evidence")
    artifacts = {"source": path, "alias": path}
    cold = baseline.capture_baseline(repository, artifacts, fingerprint_cache=tmp_path / "cache")
    events = []
    warm = baseline.capture_baseline(repository, artifacts, fingerprint_cache=tmp_path / "cache",
                                     on_progress=events.append)
    full = baseline.capture_baseline(repository, artifacts)
    assert cold == warm == full
    assert len([event for event in events if event["stage"] == "fingerprint_reused"]) == 1
