"""Batch workers release source analyses without changing compiled artifacts."""

import gzip
import json

import pytest

from core_retrosynthesis import full_scale, generic_compiler


def test_shard_cleanup_preserves_compilation_and_reuse(tmp_path) -> None:
    source = tmp_path / "source.jsonl.gz"
    row = {"reaction_id": "reduction", "reaction_smiles": "CC=O>>CCO"}
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")
    config = full_scale.FullScaleBuildConfig()
    output = tmp_path / "output"
    original = full_scale.compile_operator_shard(source, output, config=config)
    assert generic_compiler._materialized_analysis.cache_info().currsize > 0
    original_payload = gzip.decompress(
        (output / "shards" / (full_scale._shard_stem(source) + ".library.json.gz"))
        .read_bytes()
    )
    rebuilt = full_scale._compile_operator_shard_job((source, output, config, True))
    assert generic_compiler._materialized_analysis.cache_info().currsize == 0
    assert rebuilt == original
    assert gzip.decompress(
        (output / "shards" / (full_scale._shard_stem(source) + ".library.json.gz"))
        .read_bytes()
    ) == original_payload
    generic_compiler._materialized_analysis(row["reaction_smiles"])
    reused = full_scale._compile_operator_shard_job((source, output, config, False))
    assert reused["_reused"] is True
    assert generic_compiler._materialized_analysis.cache_info().currsize == 0


def test_failed_shard_also_releases_cached_analyses(tmp_path, monkeypatch) -> None:
    generic_compiler._materialized_analysis("CC=O>>CCO")
    assert generic_compiler._materialized_analysis.cache_info().currsize > 0

    def fail(*args, **kwargs):
        raise RuntimeError("fixture failure")

    monkeypatch.setattr(full_scale, "compile_operator_shard", fail)
    with pytest.raises(RuntimeError, match="fixture failure"):
        full_scale._compile_operator_shard_job(
            (tmp_path / "source", tmp_path, full_scale.FullScaleBuildConfig(), False)
        )
    assert generic_compiler._materialized_analysis.cache_info().currsize == 0
