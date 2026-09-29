"""Canonical discovery indexing, evidence joins, breadth, and deterministic search."""

from dataclasses import asdict
import json
from pathlib import Path
import sqlite3

import pytest
from rdkit import Chem

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.fragment_search import compile_fragment_query, fragment_embeddings
from condition_recommender.corpus_io import canonical_source_files, file_sha256
from condition_recommender.fragment_index import build_fragment_index, open_fragment_index
from condition_recommender.fragment_search import search_fragment_precedents


REACTIONS = [
    "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]",
    "[CH3:1][O:3][CH3:4].[CH3:5][OH:6]>>[CH3:1][O:3][CH3:4].[CH2:5]=[O:6]",
    "CCBr.CO>>CCOC",
    "c1ccccc1>>Cc1ccccc1",
    "c1ccccc1>>c1ccc2ccccc2c1",
]


def write_records(path: Path, reactions: list[str] = REACTIONS) -> None:
    rows = [{"observation_id": f"obs-{i}", "reaction_id": "rxn-shared" if i < 2 else f"rxn-{i}",
             "reference_id": f"ref-{i}", "reaction_smiles": r, "admission_tier": "review",
             "admission_reasons": ["conditions_unresolved"],
             "reaction_observation": asdict(featurize_reaction(r).observation)}
            for i, r in enumerate(reactions)]
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n", "utf-8")


@pytest.fixture
def index(tmp_path):
    source = tmp_path / "records.jsonl"
    write_records(source)
    procedures = tmp_path / "procedures.jsonl"
    procedures.write_text("\n".join(json.dumps(r) for r in [
        {"observation_id": "obs-0", "reaction_id": "rxn-shared", "text": "First observation. " * 200},
        {"observation_id": "obs-1", "reaction_id": "rxn-shared", "text": "WRONG OBSERVATION"},
        {"reaction_id": "rxn-shared", "text": "Unassigned reaction procedure"},
    ]), "utf-8")
    output = tmp_path / "fragment.sqlite"
    build_fragment_index(source, output, procedure_catalog=procedures)
    return output


def test_search_review_records_with_distinct_evidence_and_exact_procedure_join(index):
    result = search_fragment_precedents(index, "COC")
    assert result["search_status"] == "complete"
    assert result["counts"]["products"] == {"value": 2, "precision": "exact"}
    hit = result["hits"][0]
    assert hit["observation_id"] == "obs-0"
    assert "constructed" in hit["relationships"]
    assert hit["admission_tier"] == "review"
    assert len(hit["procedures"]) == 2
    assert hit["procedure_match_scope"] == "exact_observation"
    assert "WRONG OBSERVATION" not in json.dumps(hit)
    text = hit["procedures"][0]["record"]["text"]
    assert "".join(chunk["text"] for chunk in text["chunks"]) == "First observation. " * 200
    assert result["relationship_groups"]["carried_through"]["value"] == 1


@pytest.mark.parametrize("query,format,topology", [
    ("COC", "smiles", "preserve_rings"), ("c1ccccc1", "smiles", "preserve_rings"),
    ("[#6]-[O,N]", "smarts", "subgraph"), ("*~*", "smarts", "subgraph"),
])
def test_screened_hits_equal_exhaustive_graph_search(index, query, format, topology):
    compiled = compile_fragment_query(query, format, topology)
    with sqlite3.connect(index) as db:
        expected = {s for (s,) in db.execute("SELECT smiles FROM products")
                    if fragment_embeddings(compiled, Chem.MolFromSmiles(s))[0]}
    result = search_fragment_precedents(index, query, format, topology, limit=10)
    actual = {m["product_smiles"] for h in result["hits"] for m in h["matches"]}
    assert expected == actual


def test_zero_hits_is_distinct_from_missing_or_stale_index(index, tmp_path):
    assert search_fragment_precedents(index, "P(=O)(O)O")["counts"]["products"]["value"] == 0
    with pytest.raises(FileNotFoundError):
        search_fragment_precedents(tmp_path / "missing.sqlite", "CO")
    with sqlite3.connect(index) as db:
        metadata = json.loads(db.execute("SELECT payload FROM metadata").fetchone()[0])
        metadata["rdkit_version"] = "incompatible"
        db.execute("UPDATE metadata SET payload=?", (json.dumps(metadata),))
    with pytest.raises(ValueError, match="incompatible"):
        open_fragment_index(index)


def test_broad_query_reports_lower_bound_without_arbitrary_top_list(index, monkeypatch):
    import condition_recommender.fragment_search as search
    policy = search.fragment_search_policy()
    monkeypatch.setattr(search, "fragment_search_policy", lambda: {**policy, "max_matched_products": 1})
    result = search.search_fragment_precedents(index, "C", topology="subgraph")
    assert result["search_status"] == "too_broad"
    assert result["counts"]["products"] == {"value": 2, "precision": "at_least"}
    assert result["hits"] == [] and result["refinement_hints"]


def test_observation_limit_is_partial(index, monkeypatch):
    import condition_recommender.fragment_search as search
    policy = search.fragment_search_policy()
    monkeypatch.setattr(search, "fragment_search_policy", lambda: {**policy, "max_observations": 1})
    result = search.search_fragment_precedents(index, "COC")
    assert result["search_status"] == "partial"
    assert result["ranking_scope"] == "examined_subset"
    assert result["stop_reason"] == "observation_limit"
    assert result["counts"]["observations"]["precision"] == "at_least"


def test_prefix_pilot_and_atomic_failure_are_explicit(tmp_path):
    source, output = tmp_path / "records.jsonl", tmp_path / "index.sqlite"
    write_records(source)
    manifest = build_fragment_index(source, output, max_records=1)
    assert manifest["source_scope"] == "prefix_pilot" and not manifest["source_coverage_complete"]
    original = output.read_bytes()
    source.write_text('{"observation_id":"ok","reaction_smiles":"CC>>CC"}\n{bad', "utf-8")
    with pytest.raises(ValueError):
        build_fragment_index(source, output)
    assert output.read_bytes() == original


def test_strict_manifest_checks_missing_incomplete_and_corrupt_shards(tmp_path):
    shard = tmp_path / "rows.jsonl"
    write_records(shard)
    path = tmp_path / "shard_manifest.json"
    data = {"artifact_type": "generic_sharded_conversion", "shards": [
        {"status": "complete", "output_path": shard.name, "output_sha256": file_sha256(shard)}]}
    path.write_text(json.dumps(data), "utf-8")
    assert canonical_source_files(path, strict=True) == (shard.resolve(),)
    data["shards"][0]["status"] = "partial"
    path.write_text(json.dumps(data), "utf-8")
    with pytest.raises(ValueError, match="Incomplete"):
        canonical_source_files(path, strict=True)
    data["shards"][0]["status"] = "complete"
    data["shards"][0]["output_sha256"] = "wrong"
    path.write_text(json.dumps(data), "utf-8")
    with pytest.raises(ValueError, match="checksum"):
        canonical_source_files(path, strict=True)


def test_multiple_product_copies_keep_constructed_and_retained_embeddings(tmp_path):
    source, output = tmp_path / "records.jsonl", tmp_path / "index.sqlite"
    write_records(source, ["[CH3:1][Br:2].[OH:3][CH3:4].[CH3:5][O:6][CH3:7]>>[CH3:1][O:3][CH3:4].[CH3:5][O:6][CH3:7]"])
    build_fragment_index(source, output)
    result = search_fragment_precedents(output, "COC")
    assert result["counts"]["observations"]["value"] == 1
    assert result["relationship_groups"]["constructed"]["value"] == 1
    assert result["relationship_groups"]["carried_through"]["value"] == 1
    assert {m["product_component_index"] for m in result["hits"][0]["matches"]} == {0, 1}


def test_missing_reaction_id_never_joins_an_unassigned_procedure(tmp_path):
    source, output, catalog = tmp_path / "rows.jsonl", tmp_path / "index.sqlite", tmp_path / "procedures.jsonl"
    source.write_text(json.dumps({"observation_id": "known-observation", "reaction_id": "",
                                  "reaction_smiles": "CCBr.O>>CCO"}), "utf-8")
    catalog.write_text(json.dumps({"reaction_id": "", "text": "UNRELATED PROCEDURE"}), "utf-8")
    build_fragment_index(source, output, procedure_catalog=catalog)
    result = search_fragment_precedents(output, "CCO")
    assert result["hits"][0]["procedures"] == []
    assert result["hits"][0]["procedure_match_scope"] is None
