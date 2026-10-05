"""Exact observation evidence lookup and bounded indexed catalogs."""

import io
import json

from condition_recommender.processed_catalog import ProcessedCatalog, build_processed_catalog
from condition_recommender.record_storage import write_object_records
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.fragment_search import search_fragment_precedents


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
