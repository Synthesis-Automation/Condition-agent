"""Verify structural evidence, source inspection and provenance across a screened extension."""

from contextlib import closing
import json
import sqlite3

import pytest

from condition_recommender.conversion.source_evidence import structural_search_exclusions
from condition_recommender.fragment_index import build_fragment_index, bound_evidence_catalog_paths
from condition_recommender.fragment_search import lookup_exact_product, search_fragment_precedents
from condition_recommender.processed_catalog import ProcessedCatalog, build_processed_catalog
from condition_recommender.screened_fragment_extension import extend_fragment_index


def _row(identity, *, paired=False, issues=(), signature=True):
    source = {"experimental_procedure": f"Procedure for {identity}."}
    if paired:
        source["source_transfer_assessment"] = {"mapping": {"issues": list(issues)}}
    return {
        "observation_id": identity, "reaction_id": "route:" + identity, "reference_id": "doi:" + identity,
        "reaction_smiles": "[CH3:1][CH2:2][Br:3].[OH2:4]>>[CH3:1][CH2:2][OH:4]",
        "source": source, "source_dataset": "cbyd" if paired else "older",
        "reaction_signature": {"schema_version": "example"} if signature else None,
    }


def test_paired_screen_distinguishes_contradiction_missing_signature_and_missing_assessment():
    issues = ("source_mapped_products_disagree", "source_duplicate_atom_maps")
    assert structural_search_exclusions(_row("conflict", paired=True, issues=issues, signature=False)) == tuple(sorted(issues))
    assert structural_search_exclusions(_row("missing", paired=True, signature=False)) == (
        "source_reaction_signature_unavailable",)
    assert structural_search_exclusions({**_row("old"), "source_dataset": "cbyd"}) == (
        "source_mapping_unavailable_or_unreadable",)
    row = _row("review", paired=True)
    row["admission_reasons"] = ["source_model_mapping_requires_review", "partial_atom_mapping"]
    assert structural_search_exclusions(row) == ()
    row["admission_reasons"].append("source_mapped_inputs_disagree")
    assert structural_search_exclusions(row) == ("source_mapped_inputs_disagree",)
    assert structural_search_exclusions(_row("old", signature=False)) == ()


@pytest.mark.parametrize("historical", [False, True])
def test_combined_index_keeps_exact_evidence_and_joins_each_selected_hit(tmp_path, monkeypatch, historical):
    indexes = []
    for name, rows in [
        ("older", [_row("old", signature=False)]),
        ("cbyd", [_row("usable", paired=True), _row("flagged", paired=True, issues=("source_duplicate_atom_maps",)),
                  _row("no-signature", paired=True, signature=False)]),
    ]:
        root = tmp_path / name
        root.mkdir()
        source = root / "records.jsonl"
        source.write_text("".join(json.dumps(row) + "\n" for row in rows), "utf-8")
        catalog = root / "catalogs.sqlite"
        index = root / "fragment_index.sqlite"
        with monkeypatch.context() as patch:
            if historical:
                # Before the screen shipped, both indexes included every projected source row.
                patch.setattr("condition_recommender.fragment_index.structural_search_exclusions", lambda _: ())
                patch.setattr("condition_recommender.processed_catalog.structural_search_exclusions", lambda _: ())
            build_processed_catalog(source, catalog)
            build_fragment_index(source, index, evidence_catalog=catalog)
        assert len(lookup_exact_product(index, "CCO")["precedents"]) == (3 if name == "cbyd" and historical else 1)
        indexes.append(index)
    directory = tmp_path / "combined"
    manifest = extend_fragment_index(*indexes, directory)
    index = directory / "fragment_index.sqlite"
    assert manifest["supplemental_observations"] == 1
    assert manifest["counts"]["source_screen_excluded"] == 2
    assert bound_evidence_catalog_paths(index) == (indexes[0].with_name("catalogs.sqlite"), directory / "catalogs.sqlite")
    hits = search_fragment_precedents(index, "CCO")["hits"]
    assert {hit["observation_id"] for hit in hits} == {"old", "usable"}
    for hit in hits:
        assert hit["record"]["source"]["experimental_procedure"] == f"Procedure for {hit['observation_id']}."
        assert hit["procedures"][0]["link_scope"] == "exact_observation"
        assert not structural_search_exclusions(hit["record"])
    assert {p["observation_id"] for p in lookup_exact_product(index, "CCO")["precedents"]} == {"old", "usable"}
    evidence = ProcessedCatalog(directory / "catalogs.sqlite")
    assert evidence.observation("flagged")["source_dataset"] == "cbyd"
    candidate = evidence.experiment_candidates(_row("usable")["reaction_smiles"])
    assert {r["observation_id"] for r in candidate.records} == {"usable"}
    with pytest.raises(FileExistsError):
        extend_fragment_index(*indexes, directory)
    with closing(sqlite3.connect(directory / "catalogs.sqlite")) as db:
        db.execute("UPDATE metadata SET payload=json_set(payload, '$.catalog_id', 'unexpected')")
        db.commit()
    with pytest.raises(ValueError, match="binding mismatch"):
        bound_evidence_catalog_paths(index)
