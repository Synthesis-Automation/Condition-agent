"""Regression tests for lossless SciFinder RDF extraction and combination."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


PROCESSOR_DIR = Path(__file__).resolve().parents[2] / "data-processor"
if str(PROCESSOR_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSOR_DIR))

from process_reactions import parse_rdf  # noqa: E402
from Scifinder_rdf_processer import (  # noqa: E402
    RDFWorker,
    ReactionMarkdownGenerator,
    merge_rdf_records,
)


def test_parse_rdf_preserves_all_fields_and_decimal_yield(tmp_path: Path) -> None:
    rdf_path = tmp_path / "sample.rdf"
    rdf_path.write_text(
        """$RDFILE 1
$RFMT
$RXN

      TEST


  1  1
$MOL
reactant
  test

  0  0  0  0  0  0  0  0  0  0  1 V2000
M  END
$MOL
product
  test

  0  0  0  0  0  0  0  0  0  0  1 V2000
M  END
$DTYPE RXN:RCT(1):CAS_RN
$DATUM 64-17-5
$DTYPE RXN:RCT(1):AMD
$DATUM reactant-amd
$DTYPE RXN:PRO(1):CAS_RN
$DATUM 67-56-1
$DTYPE RXN:VAR(1):CAS_Reaction_Number
$DATUM 12345
$DTYPE RXN:VAR(1):PRO(1):YIELD
$DATUM 87.5
$DTYPE RXN:VAR(1):EXP_PROC
$DATUM Heat for 2 h
then cool to room temperature
$DTYPE RXN:VAR(1):UNRECOGNIZED_CHEMISTRY_FIELD
$DATUM retained value
""",
        encoding="utf-8",
    )

    record = parse_rdf(str(rdf_path))["12345"]

    assert record["yield_pct"] == 87.5
    assert record["pro_yields"] == {1: 87.5}
    assert record["rct_cas"] == ["64-17-5"]
    assert record["rct_amd"] == ["reactant-amd"]
    assert record["exp_proc"] == ["Heat for 2 h then cool to room temperature"]
    assert record["raw_fields"][
        "RXN:VAR(1):UNRECOGNIZED_CHEMISTRY_FIELD"
    ] == ["retained value"]
    assert len(record["rct_mol"]) == 1
    assert len(record["pro_mol"]) == 1


def test_duplicate_records_merge_values_and_provenance() -> None:
    primary = {
        "rct_cas": ["1-00-0"],
        "source_files": ["a.rdf"],
        "reaction_types": ["Type A"],
        "duplicate_reaction_ids": ["10"],
        "pro_yields": {1: 80},
        "title": "First title",
        "raw_fields": {"FIELD": ["one"]},
    }
    duplicate = {
        "rct_cas": ["2-00-0"],
        "source_files": ["b.rdf"],
        "reaction_types": ["Type B"],
        "duplicate_reaction_ids": ["11"],
        "pro_yields": {1: 81.5, 2: 50},
        "title": "Second title",
        "raw_fields": {"FIELD": ["two"]},
    }

    merged = merge_rdf_records(primary, duplicate)

    assert merged["rct_cas"] == ["1-00-0", "2-00-0"]
    assert merged["source_files"] == ["a.rdf", "b.rdf"]
    assert merged["duplicate_reaction_ids"] == ["10", "11"]
    assert merged["pro_yields"] == {1: 80, 2: 50}
    assert merged["product_yield_values"] == {"1": [80, 81.5]}
    assert merged["field_values"]["title"] == ["First title", "Second title"]
    assert merged["raw_fields"]["FIELD"] == ["one", "two"]


def test_exact_cross_id_duplicates_are_removed_but_provenance_is_merged() -> None:
    worker = RDFWorker(".", "", "", combine_all=True)
    rows = [
        {
            "ReactionID": "10",
            "ReactionType": "Type A",
            "ReactantSMILES": "CC.O",
            "ProductSMILES": "CCO",
            "Temperature_C": 25,
            "Time_h": 2,
        },
        {
            "ReactionID": "11",
            "ReactionType": "Type B",
            "ReactantSMILES": "O.CC",
            "ProductSMILES": "CCO",
            "Temperature_C": 25,
            "Time_h": 2,
        },
    ]
    common = {
        "rct_cas": ["64-17-5", "7732-18-5"],
        "pro_cas": ["64-17-5"],
        "pro_yields": {1: 80},
        "notes": ["same"],
    }
    rdf_map = {
        "10": {
            **common,
            "source_files": ["a.rdf"],
            "reaction_types": ["Type A"],
            "duplicate_reaction_ids": ["10"],
        },
        "11": {
            **common,
            "source_files": ["b.rdf"],
            "reaction_types": ["Type B"],
            "duplicate_reaction_ids": ["11"],
        },
    }

    unique, removed = worker._deduplicate_rows(rows, rdf_map)

    assert removed == 1
    assert [row["ReactionID"] for row in unique] == ["10"]
    assert unique[0]["ReactionType"] == "Type A | Type B"
    assert rdf_map["10"]["duplicate_reaction_ids"] == ["10", "11"]
    assert rdf_map["10"]["source_files"] == ["a.rdf", "b.rdf"]
    assert "11" not in rdf_map


def test_csv_export_keeps_long_reactions_and_complete_metadata(tmp_path: Path) -> None:
    long_smiles = "C" * 170
    rows = [{
        "ReactionID": "10",
        "ReactionType": "Test",
        "Yield_%": 75.5,
        "Temperature_C": 20,
        "Time_h": 1,
        "ReactantSMILES": long_smiles,
        "ProductSMILES": "C",
        "Reference": "Title | Author | Citation",
    }]
    rdf_map = {"10": {
        "rct_cas": ["1-00-0"],
        "pro_cas": ["2-00-0"],
        "pro_yields": {1: 75.5},
        "title": "Title",
        "authors": "Author",
        "citation": "Citation",
        "source_files": ["source.rdf"],
        "duplicate_reaction_ids": ["10"],
        "raw_fields": {"CUSTOM": ["value"]},
    }}
    output_path = tmp_path / "out.csv"

    ReactionMarkdownGenerator().generate_csv_export(
        rows, str(output_path), "test", rdf_map=rdf_map
    )

    with output_path.open(encoding="utf-8", newline="") as input_file:
        exported = list(csv.DictReader(input_file))
    assert len(exported) == 1
    assert exported[0]["reaction_smiles"] == f"{long_smiles}>>C"
    assert exported[0]["yield_pct"] == "75.5"
    assert exported[0]["title"] == "Title"
    assert exported[0]["product_yields_json"] == '{"1": 75.5}'
    assert "raw_fields_json" not in exported[0]
    assert "conflicting_values_json" not in exported[0]


def test_worker_output_generation_does_not_create_markdown(tmp_path: Path) -> None:
    markdown_path = tmp_path / "must_not_exist.md"
    csv_path = tmp_path / "output.csv"
    worker = RDFWorker(
        str(tmp_path),
        str(markdown_path),
        str(csv_path),
    )
    rows = [{
        "ReactionID": "10",
        "ReactionType": "Test",
        "Yield_%": 80,
        "Temperature_C": "",
        "Time_h": "",
        "ReactantSMILES": "C",
        "ProductSMILES": "CO",
        "Reference": "",
    }]

    worker._generate_outputs(rows, {"10": {"pro_yields": {1: 80}}})

    assert csv_path.exists()
    assert not markdown_path.exists()
