"""Resolve and pin one complete, coherent processed dataset release."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

from .corpus_io import file_sha256
from .conversion.atomic import atomic_json

DEFAULT_PROCESSED_ROOT = Path(__file__).resolve().parents[1] / "datasets" / "processed_datasets"


@dataclass(frozen=True)
class ProcessedRelease:
    """An immutable manifest and its root, resolved once by a consumer."""

    root: Path
    manifest: dict[str, Any]

    def artifact(self, name: str) -> Path:
        """Resolve a declared artifact within the pinned release."""
        entry = self.manifest["artifacts"].get(name)
        if not entry:
            raise FileNotFoundError(f"Processed release capability is unavailable: {name}")
        path = (self.root / entry["relative_path"]).resolve()
        if not path.is_relative_to(self.root.resolve()):
            raise ValueError("Processed artifact escapes its release")
        if not path.is_file() or path.stat().st_size != entry["size_bytes"]:
            raise ValueError(f"Processed artifact is missing or inconsistent: {name}")
        return path


def resolve_processed_release(root: str | Path = DEFAULT_PROCESSED_ROOT,
                              *, verify_artifacts: bool = False) -> ProcessedRelease:
    """Resolve the active manifest without silently falling back to an old corpus."""
    selected = Path(root).resolve()
    manifest_path = selected if selected.is_file() else selected / "manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if payload.get("schema_version") == "processed_reaction_active.v1":
        target = (manifest_path.parent / payload["release_manifest"]).resolve()
        if not target.is_relative_to(manifest_path.parent.resolve()):
            raise ValueError("Active release reference escapes the dataset root")
        if file_sha256(target) != payload["sha256"]:
            raise ValueError("Active processed manifest checksum mismatch")
        manifest_path = target
        payload = json.loads(target.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "processed_reaction_release.v1" or not payload.get("build_complete"):
        raise ValueError("Processed dataset release is incompatible or unfinished")
    result = ProcessedRelease(manifest_path.parent, payload)
    for name, entry in payload.get("artifacts", {}).items():
        path = result.artifact(name)
        if verify_artifacts and file_sha256(path) != entry["sha256"]:
            raise ValueError(f"Processed artifact checksum mismatch: {name}")
    return result


def processed_artifact(name: str, root: str | Path = DEFAULT_PROCESSED_ROOT) -> Path:
    """Return a declared artifact from the active complete release."""
    return resolve_processed_release(root).artifact(name)


def publish_processed_release(
    release_root: str | Path,
    output_root: str | Path = DEFAULT_PROCESSED_ROOT,
) -> ProcessedRelease:
    """Verify a complete release and atomically select it for app/agent readers."""
    release = resolve_processed_release(release_root, verify_artifacts=True)
    destination = Path(output_root).resolve()
    active = {
        "schema_version": "processed_reaction_active.v1",
        "release_id": release.manifest["release_id"],
        "release_manifest": (release.root / "manifest.json").relative_to(destination).as_posix(),
        "sha256": file_sha256(release.root / "manifest.json"),
    }
    atomic_json(destination / "manifest.json", active)
    return release
