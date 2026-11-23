"""Experiment evidence retains controls, independent of recipe admission."""

from contextlib import closing
import json
from pathlib import Path
import sqlite3

import pytest

from condition_recommender.experiment_catalog import build_experiment_catalog
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.processed_catalog import ProcessedCatalog, build_processed_catalog


LIU_OBSERVATIONS = Path(__file__).resolve().parents[2] / "examples/common/liu_thioether_controls.jsonl"


def test_source_experiments_round_trip_and_local_analogue(tmp_path):
    destination = tmp_path / "evidence.sqlite"
    report = build_experiment_catalog(LIU_OBSERVATIONS, destination)
    catalog = ProcessedCatalog(destination)
    inputs = [json.loads(line) for line in LIU_OBSERVATIONS.read_text().splitlines()]
    query = inputs[0]["reaction"]["reaction_smiles"]
    assert catalog.supports_experiment_analogues
    result = catalog.experiment_candidates(query)
    assert result.scope == "exact_graph" and result.lookup_complete and not result.truncated
    assert report["counts"]["experiments_with_local_core"] == len(inputs) == 6
    by_entry = {item["source_row_number"]: item for item in result.records}
    for original in inputs:
        entry = by_entry[original["source"]["source_row_number"]]
        assert entry["source"]["raw_fields"]["source_observation"] == original
        assert entry["yield_pct"] == original["outcomes"][0]["value"]
        assert entry["observation_id"] == original["observation_id"]
    assert by_entry[5]["index_eligibility"] == "ineligible"
    assert by_entry[5]["resolved_recipe"]["declared_absences"] == ["base"]
    assert by_entry[2]["resolved_recipe"]["stages"][0]["provenance"]["illumination"] == "no light"
    assert by_entry[2]["yield_pct"] is None  # source reports Trace

    ethyl_analogue = query.replace("CC(=O)c1", "CCC(=O)c1")
    assert not catalog.experiment_candidates(ethyl_analogue).records
    related = catalog.experiment_candidates(ethyl_analogue, include_analogues=True)
    assert related.scope == "exact_and_indexed_local_cores"
    assert related.records == result.records
    capped = catalog.experiment_candidates(query, limit=2, include_analogues=True)
    assert capped.truncated and not capped.lookup_complete
    assert capped.records == result.records[:2]


def test_unknown_outcome_and_unspecified_conditions_survive_json_import(tmp_path):
    original = json.loads(LIU_OBSERVATIONS.read_text().splitlines()[0])
    original.update(outcomes=[])
    original["conditions"].update(stages=[], procedure_text="", declared_absences=[])
    original["raw_fields"] = {}
    source = tmp_path / "source.json"
    source.write_text(json.dumps([original]))
    path = tmp_path / "import.sqlite"
    build_experiment_catalog(source, path)
    row, = ProcessedCatalog(path).experiment_candidates(original["reaction"]["reaction_smiles"]).records
    assert row["yield_pct"] is None
    assert not row["resolved_recipe"]["declared_absences"]
    assert row["source"]["raw_fields"]["source_observation"]["outcomes"] == []
    source.write_text(json.dumps([original, original]))
    before = path.read_bytes()
    with pytest.raises(sqlite3.IntegrityError):
        build_experiment_catalog(source, path)
    assert path.read_bytes() == before


def test_legacy_exact_lookup_uses_bound_fragment_and_both_reaction_sides(tmp_path):
    source = tmp_path / "converted.jsonl"
    smiles = ("CCBr.O>>CCO", "O.CCBr>>CCO", "CCCl.O>>CCO")
    rows = [dict(observation_id=f"obs-{i}", reaction_id=f"step-{i}",
                 reaction_smiles=sm, source_dataset="synthetic", source={}, yield_pct=0)
            for i, sm in enumerate(smiles)]
    source.write_text("\n".join(map(json.dumps, rows)))
    path = tmp_path / "legacy.sqlite"
    report = build_processed_catalog(source, path)
    # A v2 release has all observation payloads and no reaction lookup table.
    with closing(sqlite3.connect(path)) as db:
        report["schema_version"] = "processed_evidence_catalog.v2"
        db.execute("UPDATE metadata SET payload=?", (json.dumps(report),))
        db.execute("DROP TABLE experiment_lookup")
        db.commit()
    fragment = tmp_path / "fragment.sqlite"
    build_fragment_index(source, fragment, evidence_catalog=path)
    legacy = ProcessedCatalog(path)
    assert not legacy.supports_experiment_analogues
    assert not legacy.experiment_candidates(smiles[0]).lookup_complete
    with pytest.raises(ValueError, match="bound fragment index"):
        legacy.experiment_candidates(smiles[0], include_analogues=True)
    result = legacy.experiment_candidates("O.[CH3:1]CBr>O>CCO", fragment_index=fragment)
    assert result.scope == "bound_fragment_exact_graph" and result.lookup_complete
    assert result.candidates_examined == 3
    assert {r["observation_id"] for r in result.records} == {"obs-0", "obs-1"}
    assert legacy.experiment_candidates(smiles[0], fragment_index=fragment, limit=1).truncated
    mismatched = tmp_path / "other.sqlite"
    with closing(sqlite3.connect(path)) as original, closing(sqlite3.connect(mismatched)) as other:
        original.backup(other)
        report["catalog_id"] = "unrelated"
        other.execute("UPDATE metadata SET payload=?", (json.dumps(report),))
        other.commit()
    with pytest.raises(ValueError, match="not bound"):
        ProcessedCatalog(mismatched).experiment_candidates(smiles[0], fragment_index=fragment)


def test_exact_only_unverified_and_replaced_catalog(tmp_path):
    # An incomplete source mapping and a missing target cannot yield local keys.
    source = tmp_path / "records.jsonl"
    source.write_text("\n".join(json.dumps(dict(observation_id=f"sample-{i}", reaction_smiles=sm))
                                for i, sm in enumerate(("[CH3:1]Br>>CO", "CBr>>"))))
    path = tmp_path / "catalog.sqlite"
    build_processed_catalog(source, path, include_local_cores=False)
    exact = ProcessedCatalog(path)
    assert not exact.supports_experiment_analogues
    with pytest.raises(ValueError, match="Local experiment lookup unavailable"):
        exact.experiment_candidates("CBr>>CO", include_analogues=True)
    # Whole graph identity supplies a candidate only, and the source map survives.
    row, = exact.experiment_candidates("CBr>>CO").records
    assert row["reaction_smiles"] == "[CH3:1]Br>>CO"
    invalid = exact.experiment_candidates("CBr>>")
    assert not invalid.lookup_complete and not invalid.records
    report = build_processed_catalog(source, path)
    assert report["counts"]["experiments_with_local_core"] == 0
    assert report["counts"]["experiments_with_graph"] == 1
    with pytest.raises(ValueError, match="changed"):
        exact.observation("sample-0")
    with closing(sqlite3.connect(path)) as db:
        report["experiment_local_core_definition"] = "outdated"
        db.execute("UPDATE metadata SET payload=?", (json.dumps(report),))
        db.commit()
    with pytest.raises(ValueError, match="definition mismatch"):
        ProcessedCatalog(path).experiment_candidates("CBr>>CO", include_analogues=True)
