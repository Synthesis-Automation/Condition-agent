"""Disposable local checksum cache; full replay verification never consults it."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


def stat_token(stat: os.stat_result) -> dict[str, int]:
    """Identify a file instance and ordinary content/metadata changes."""
    return {key: int(getattr(stat, key)) for key in (
        "st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns",
    )}


class FingerprintCache:
    """Reuse completed hashes while file identity and stat metadata remain equal.

    This is an I/O optimization for local files, not proof against edits that
    preserve all metadata. Explicit full verification must still read the bytes.
    Invalid or inaccessible cache entries are misses; cache writes are atomic.
    """

    def __init__(self, directory: Path) -> None:
        self.directory = directory

    def _entry_path(self, path: Path) -> Path:
        key = hashlib.sha256(str(path).encode("utf-8")).hexdigest()
        return self.directory / f"{key}.json"

    def get(self, path: Path, stat: os.stat_result) -> str | None:
        """Return a validated cache entry for this exact file snapshot, if any."""
        try:
            record = json.loads(self._entry_path(path).read_text(encoding="utf-8"))
            if not isinstance(record, dict):
                return None
            digest = record.get("sha256")
            if (record.get("schema_version") == "artifact_fingerprint_cache.v1"
                    and record.get("path") == str(path)
                    and record.get("stat") == stat_token(stat)
                    and isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest)):
                return digest
        except (OSError, ValueError):
            pass
        return None

    def put(self, path: Path, stat: os.stat_result, digest: str) -> None:
        """Save only a completed hash; cache unavailability does not block work."""
        temporary = None
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            record = {"schema_version": "artifact_fingerprint_cache.v1",
                      "path": str(path), "stat": stat_token(stat), "sha256": digest}
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=self.directory,
                prefix=".fingerprint-", suffix=".tmp", delete=False,
            ) as stream:
                temporary = Path(stream.name)
                json.dump(record, stream, sort_keys=True)
            os.replace(temporary, self._entry_path(path))
        except OSError:
            pass
        finally:
            if temporary is not None:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError:
                    pass
