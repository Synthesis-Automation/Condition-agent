"""Exact observation evidence lookup and bounded indexed catalogs."""

import io
import json
from dataclasses import asdict

from condition_recommender.processed_catalog import ProcessedCatalog, build_processed_catalog
from condition_recommender.record_storage import write_object_records
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.fragment_search import search_fragment_precedents
from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_core_graph import shared_core_availability
from tests.reactive_taxonomy.test_shared_core_availability import dense_reduction


def test_over_budget_graph_preserves_exact_catalog_lookup_and_coverage(tmp_path):
    reaction = dense_reduction()
    analysis = featurize_reaction(reaction)
    row = {"observation_id": "dense", "reaction_id": "rx-dense", "reaction_smiles": reaction,
           "reaction_signature": asdict(analysis.reaction_signature), "reaction_core": asdict(analysis.reaction_core),
           "reaction_completeness": {"status": "verified"}}
    source, destination = tmp_path / "rows.jsonl", tmp_path / "catalog.sqlite"
    with source.open("w", encoding="utf-8") as out:
        write_object_records(out, [row, {**row, "observation_id": "dense-2"}])
    report = build_processed_catalog(source, destination)
    assert report["local_core_limit_reasons"] == {"LOCAL_CORE_CENTER_DISTANCE_PAIR_LIMIT": 2}
    assert report["counts"]["experiments_with_local_core"] == 0
    assert report["counts"]["experiments_with_graph"] == 2
    assert report["experiment_local_core_availability"] == shared_core_availability().fingerprint
    catalog = ProcessedCatalog(destination)
    assert catalog.supports_experiment_analogues
    found = catalog.experiment_candidates(reaction, limit=10, include_analogues=False)
    assert {r["observation_id"] for r in found.records} == {"dense", "dense-2"}


def test_shared_catalog_preserves_and_indexes_source_evidence(tmp_path):
    rows = [{"observation_id": f"obs-{i}", "reaction_id": "shared",
             "reaction_smiles": "CCBr>O>CCO", "source_dataset": "test",
             "reference_id": "ref", "reference_identity": {"reference_id": "ref", "raw": "paper"},
             "source": {"experimental_procedure": f"Procedure {i}", "raw_fields": {"original": "x" * 800}},
             "admission_tier": "review", "admission_reasons": ["ambiguous_mapping"]}
            for i in range(2)]
    source = tmp_path / "records.jsonl"
    with source.open("w", encoding="utf-8") as f:
        write_object_records(f, rows)
    catalog_path = tmp_path / "catalogs.sqlite"
    report = build_processed_catalog(source, catalog_path)
    catalog = ProcessedCatalog(catalog_path)
    assert report["counts"]["observations"] == 2
    assert catalog.observation("obs-0") == rows[0]
    assert catalog.references(["ref", "missing"])["ref"]["raw"] == "paper"
    page = catalog.procedures(["shared"], limit=1)
    assert page["total"] == 2 and page["next_offset"] == 1
    assert page["records"][0]["procedure_text"] == "Procedure 0"
    index = tmp_path / "fragment.sqlite"
    build_fragment_index(source, index, evidence_catalog=catalog_path)
    result = search_fragment_precedents(index, "CCBr", search_side="reactant")
    hit = result["hits"][0]
    assert hit["procedure_match_scope"] == "exact_observation"
    assert hit["procedures"][0]["record"]["procedure_text"] == "Procedure 0"
    assert hit["record"]["source"]["experimental_procedure"] == "Procedure 0"
