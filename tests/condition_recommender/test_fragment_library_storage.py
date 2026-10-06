"""Large fragment libraries stay bounded and failed packing retains its scan."""

from dataclasses import asdict
import json
import random
import sqlite3

import pytest

from condition_recommender import fragment_index as storage
from condition_recommender.fragment_search import search_fragment_precedents
from reactive_taxonomy import featurize_reaction


def library_database():
    db = sqlite3.connect(":memory:")
    db.execute("CREATE TABLE library (ordinal INTEGER PRIMARY KEY,payload BLOB,raw_size INTEGER,sha256 TEXT)")
    return db


def test_library_chunks_support_payload_larger_than_sqlite_blob_limit():
    payload = random.Random(42).randbytes(8192)
    with library_database() as db:
        db.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 1024)
        with pytest.raises(sqlite3.DataError):
            db.execute("INSERT INTO library VALUES (0,?,0,'')", (payload,))
        manifest = storage._write_serialized_library(db, payload, chunk_bytes=512)
        assert manifest["chunk_count"] == 16
        assert storage._read_serialized_library(db, manifest) == payload


@pytest.mark.parametrize("damage", ["missing", "checksum", "order", "payload"])
def test_library_chunks_reject_incomplete_or_corrupt_evidence(damage):
    with library_database() as db:
        manifest = storage._write_serialized_library(db, b"data" * 1000, chunk_bytes=512)
        if damage == "missing":
            db.execute("DELETE FROM library WHERE ordinal=7")
        elif damage == "checksum":
            db.execute("UPDATE library SET sha256='incorrect' WHERE ordinal=0")
        elif damage == "order":
            db.execute("UPDATE library SET ordinal=99 WHERE ordinal=0")
        else:
            db.execute("UPDATE library SET payload=? WHERE ordinal=0", (b"corrupt",))
        with pytest.raises(ValueError, match="checksum|incomplete|order|corrupt"):
            storage._read_serialized_library(db, manifest)


def test_packing_failure_resumes_completed_scan_without_reading_observations(tmp_path, monkeypatch):
    source = tmp_path / "records.jsonl"
    reaction = "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]"
    source.write_text(json.dumps({
        "observation_id": "observation", "reaction_id": "reaction", "reaction_smiles": reaction,
        "reaction_observation": asdict(featurize_reaction(reaction).observation),
    }) + "\n", encoding="utf-8")
    baseline = storage.build_fragment_index(source, tmp_path / "baseline.sqlite")
    destination = tmp_path / "resume.sqlite"
    writer = storage._write_serialized_library

    def fail_packing(*args, **kwargs):
        raise RuntimeError("packing failure")

    monkeypatch.setattr(storage, "_write_serialized_library", fail_packing)
    with pytest.raises(RuntimeError, match="packing failure"):
        storage.build_fragment_index(source, destination)
    assert not destination.exists()
    assert destination.with_name(destination.name + ".building").is_file()
    monkeypatch.setattr(storage, "_write_serialized_library", writer)

    def unexpected_scan(*args, **kwargs):
        raise AssertionError("Completed observations must not be scanned again")

    monkeypatch.setattr(storage, "iter_canonical_records", unexpected_scan)
    resumed = storage.build_fragment_index(source, destination)
    assert resumed == baseline
    assert not destination.with_name(destination.name + ".building").exists()
    result = search_fragment_precedents(destination, "COC")
    assert result["hits"][0]["observation_id"] == "observation"
    assert result["hits"][0]["matched_sides"] == ["product"]
