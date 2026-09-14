"""Prepared forward artifacts must not mask a newer paired retro library."""

import os

import pytest

import app.web_api.runtime as runtime_module


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
