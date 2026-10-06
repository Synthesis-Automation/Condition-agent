"""Lossless, deterministic storage and corrupt reference regressions."""

import io
import json

import pytest

from condition_recommender.record_storage import iter_record_shard, write_object_records


def test_shared_objects_round_trip_and_are_not_aliased(tmp_path):
    core = {"atoms": [{"atom_index": i, "environment": "C" * 30} for i in range(40)]}
    rows = [{"reaction_core": core, "reaction_observation": {"core": core},
             "warning": "ambiguous", "raw": {"blank": None}},
            {"reaction_core": core, "reaction_observation": {"core": core}}]
    first, second = io.StringIO(), io.StringIO()
    assert write_object_records(first, rows) == 2
    write_object_records(second, rows)
    assert first.getvalue() == second.getvalue()
    path = tmp_path / "records.jsonl"
    path.write_text(first.getvalue(), encoding="utf-8")
    restored = list(iter_record_shard(path))
    assert restored == rows
    restored[0]["reaction_core"]["atoms"][0]["atom_index"] = 99
    assert restored[0]["reaction_observation"]["core"]["atoms"][0]["atom_index"] == 0
    assert restored[1]["reaction_core"]["atoms"][0]["atom_index"] == 0
    assert len(first.getvalue()) < len(json.dumps(rows))


def test_dangling_and_corrupt_objects_are_rejected(tmp_path):
    path = tmp_path / "records.jsonl"
    path.write_text('{"storage_schema":"reaction_object_shard.v1"}\n'
                    '{"record":{"core":{"$reaction_object":"missing"}}}\n')
    with pytest.raises(ValueError, match="Missing reaction object"):
        list(iter_record_shard(path))
    path.write_text('{"storage_schema":"reaction_object_shard.v1"}\n'
                    '{"object_id":"wrong","value":{"evidence":"unchanged"}}\n')
    with pytest.raises(ValueError, match="checksum"):
        list(iter_record_shard(path))


def test_ordinary_jsonl_remains_readable(tmp_path):
    path = tmp_path / "source.jsonl"
    path.write_text('{"reaction_id":"one","warnings":["unresolved"]}\n')
    assert list(iter_record_shard(path))[0]["warnings"] == ["unresolved"]


def test_selected_fields_keep_object_and_reference_validation(tmp_path, monkeypatch):
    import condition_recommender.record_storage as storage

    path = tmp_path / "records.jsonl"
    rows = [{"observation_id": "one", "evidence": {"atoms": list(range(500))}}]
    handle = io.StringIO()
    write_object_records(handle, rows)
    path.write_text(handle.getvalue(), encoding="utf-8")
    original = storage.hydrate

    def guarded_hydrate(value, objects):
        assert not (isinstance(value, dict) and "evidence" in value)
        return original(value, objects)

    monkeypatch.setattr(storage, "hydrate", guarded_hydrate)
    assert list(iter_record_shard(path, fields=("observation_id",))) == [{"observation_id": "one"}]

    path.write_text('{"storage_schema":"reaction_object_shard.v1"}\n'
                    '{"record":{"observation_id":"one","omitted":{"$reaction_object":"missing"}}}\n')
    with pytest.raises(ValueError, match="missing reaction object reference"):
        list(iter_record_shard(path, fields=("observation_id",)))

    path.write_text('{"storage_schema":"reaction_object_shard.v1"}\n'
                    '{"object_id":"wrong","value":{"evidence":"unchanged"}}\n'
                    '{"record":{"observation_id":"one"}}\n')
    with pytest.raises(ValueError, match="checksum"):
        list(iter_record_shard(path, fields=("observation_id",)))


def test_selected_fields_resolve_references_and_reset_between_members(tmp_path):
    import gzip

    path = tmp_path / "records.jsonl.gz"
    evidence = {"atoms": list(range(500))}
    for mode, identity in (("wt", "one"), ("at", "two")):
        with gzip.open(path, mode, encoding="utf-8") as handle:
            write_object_records(handle, [{"observation_id": identity, "evidence": evidence}])
    assert list(iter_record_shard(path, fields=("evidence",))) == [
        {"evidence": evidence}, {"evidence": evidence},
    ]
