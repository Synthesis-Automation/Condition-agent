"""Publication binding, relocation, and incomplete-release protection."""

import json

import pytest

from condition_recommender.corpus_io import file_sha256
from condition_recommender.processed_release import resolve_processed_release


def _release(tmp_path):
    directory = tmp_path / "releases" / "id"
    directory.mkdir(parents=True)
    artifact = directory / "index.sqlite"
    artifact.write_bytes(b"immutable")
    manifest = {"schema_version": "processed_reaction_release.v1", "release_id": "id",
                "build_complete": True, "artifacts": {"condition_index": {
                    "relative_path": artifact.name, "sha256": file_sha256(artifact),
                    "size_bytes": artifact.stat().st_size}}}
    path = directory / "manifest.json"
    path.write_text(json.dumps(manifest))
    (tmp_path / "manifest.json").write_text(json.dumps({"schema_version": "processed_reaction_active.v1",
        "release_manifest": "releases/id/manifest.json", "sha256": file_sha256(path)}))
    return path, manifest


def test_release_pins_identity_and_detects_same_size_corruption(tmp_path):
    path, _ = _release(tmp_path)
    release = resolve_processed_release(tmp_path, verify_artifacts=True)
    assert release.artifact("condition_index").parent == path.parent
    release.artifact("condition_index").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        resolve_processed_release(tmp_path, verify_artifacts=True)


def test_unfinished_and_escaping_release_artifacts_are_rejected(tmp_path):
    path, manifest = _release(tmp_path)
    manifest["build_complete"] = False
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="unfinished"):
        resolve_processed_release(path)
    manifest["build_complete"] = True
    manifest["artifacts"]["condition_index"]["relative_path"] = "../../outside.sqlite"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="escapes"):
        resolve_processed_release(path)
