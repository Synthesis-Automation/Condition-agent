"""Prepared forward artifacts must not mask a newer paired retro library."""

import os
import json

import pytest

import app.web_api.runtime as runtime_module
from condition_recommender.corpus_io import file_sha256
from condition_recommender.processed_release import publish_processed_release


def _published_root(root):
    directory = root / "releases" / "test"
    directory.mkdir(parents=True)
    artifacts = {}
    for name in ("condition_index", "fragment_index", "retro_library", "forward_library"):
        path = directory / name
        path.write_bytes(name.encode())
        artifacts[name] = {"relative_path": name, "sha256": file_sha256(path), "size_bytes": path.stat().st_size}
    (directory / "manifest.json").write_text(json.dumps({
        "schema_version": "processed_reaction_release.v1", "release_id": "test",
        "build_complete": True, "artifacts": artifacts,
    }))
    publish_processed_release(directory, root)
    return directory


@pytest.mark.parametrize("use_environment", [False, True])
def test_custom_operator_root_overrides_published_condition_defaults(tmp_path, monkeypatch, use_environment):
    root = tmp_path / "processed"
    published = _published_root(root)
    monkeypatch.setattr(runtime_module, "DEFAULT_LIBRARY_ROOT", root)
    custom = tmp_path / "custom"
    kwargs = {}
    if use_environment:
        monkeypatch.setenv("CORE_RETROSYNTHESIS_LIBRARY_ROOT", str(custom))
    else:
        kwargs["retrosynthesis_library_root"] = custom
    runtime = runtime_module.LocalRecommendationRuntime(**kwargs)
    assert runtime.index_path == published / "condition_index"
    assert runtime.fragment_index_path == published / "fragment_index"
    assert runtime._retrosynthesis_library_path("full") == custom / "full/operator_library_v3.json.gz"
    assert runtime._forward_library_path("full") == custom / "full/forward_operator_library_v1.json.gz"


def test_explicit_processed_operator_root_is_pinned(tmp_path, monkeypatch):
    root = tmp_path / "processed"
    _published_root(root)
    monkeypatch.setattr(runtime_module, "DEFAULT_LIBRARY_ROOT", root)
    custom = tmp_path / "custom"
    published = _published_root(custom)
    runtime = runtime_module.LocalRecommendationRuntime(retrosynthesis_library_root=custom)
    assert runtime._retrosynthesis_library_path("full") == published / "retro_library"
    assert runtime._forward_library_path("full") == published / "forward_library"


def test_stale_forward_artifact_rebuilds_and_tracks_source_changes(
    tmp_path, monkeypatch,
) -> None:
    mode = tmp_path / "compact"
    mode.mkdir()
    retro = mode / "operator_library_v3.json.gz"
    prepared = mode / "forward_operator_library_v1.json.gz"
    retro.touch()
    prepared.touch()
    os.utime(prepared, ns=(10_000_000_000, 10_000_000_000))
    os.utime(retro, ns=(20_000_000_000, 20_000_000_000))
    runtime = runtime_module.LocalRecommendationRuntime(
        retrosynthesis_library_root=tmp_path,
    )
    builds = []

    def build(source):
        result = object()
        builds.append(result)
        return result

    monkeypatch.setattr(runtime, "_get_retrosynthesis_library", lambda mode: object())
    monkeypatch.setattr(runtime_module, "build_forward_library", build)
    monkeypatch.setattr(
        runtime_module, "load_forward_library",
        lambda path: pytest.fail("stale prepared artifact must not be loaded"),
    )
    first = runtime._get_forward_library("compact")
    assert runtime._get_forward_library("compact") is first
    assert len(builds) == 1
    assert runtime.capabilities()["forward_library_modes"]["compact"]["prepared"] is False
    os.utime(retro, ns=(30_000_000_000, 30_000_000_000))
    assert runtime._get_forward_library("compact") is not first
    assert len(builds) == 2


@pytest.mark.parametrize("paired_source_exists", [False, True])
def test_current_or_standalone_prepared_forward_artifact_is_loaded(
    tmp_path, monkeypatch, paired_source_exists,
) -> None:
    mode = tmp_path / "compact"
    mode.mkdir()
    prepared = mode / "forward_operator_library_v1.json.gz"
    prepared.touch()
    if paired_source_exists:
        retro = mode / "operator_library_v3.json.gz"
        retro.touch()
        os.utime(retro, ns=(10_000_000_000, 10_000_000_000))
    os.utime(prepared, ns=(20_000_000_000, 20_000_000_000))
    runtime = runtime_module.LocalRecommendationRuntime(
        retrosynthesis_library_root=tmp_path,
    )
    expected = object()
    monkeypatch.setattr(runtime_module, "load_forward_library", lambda path: expected)
    monkeypatch.setattr(
        runtime_module, "build_forward_library",
        lambda source: pytest.fail("current artifact should be loaded directly"),
    )
    assert runtime._get_forward_library("compact") is expected
    assert runtime._prepared_forward_library_current("compact") is True
