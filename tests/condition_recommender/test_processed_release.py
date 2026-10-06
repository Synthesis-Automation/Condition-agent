"""Publication binding, relocation, and incomplete-release protection."""

import json

import pytest

from condition_recommender.corpus_io import file_sha256
from condition_recommender.processed_release import publish_processed_release, resolve_processed_release


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


def test_publication_restores_missing_active_pointer_after_completed_manifest(tmp_path):
    path, _ = _release(tmp_path)
    (tmp_path / "manifest.json").unlink()
    pinned = publish_processed_release(path.parent, tmp_path)
    active = resolve_processed_release(tmp_path, verify_artifacts=True)
    assert active.root == pinned.root == path.parent
    assert active.manifest["release_id"] == "id"


def test_failed_publication_preserves_existing_active_pointer(tmp_path):
    path, _ = _release(tmp_path)
    active = tmp_path / "manifest.json"
    before = active.read_bytes()
    (path.parent / "index.sqlite").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        publish_processed_release(path.parent, tmp_path)
    assert active.read_bytes() == before


def test_builder_resume_publishes_already_complete_release(tmp_path, monkeypatch):
    from app import processed_dataset_builder as builder

    path, manifest = _release(tmp_path)
    manifest["builder_version"] = "processed_pipeline.v4"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    (tmp_path / "manifest.json").unlink()
    monkeypatch.setattr(builder, "prepare_release", lambda *args: (path.parent, manifest))

    def unexpected_conversion(*args, **kwargs):
        raise AssertionError("A completed release must not be converted again")

    monkeypatch.setattr(builder, "convert_release", unexpected_conversion)
    result = builder.build_processed_datasets(tmp_path / "input", tmp_path)
    assert result["release_id"] == "id"
    assert resolve_processed_release(tmp_path, verify_artifacts=True).root == path.parent
