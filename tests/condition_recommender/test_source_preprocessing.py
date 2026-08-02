import csv
import gzip
import json
from pathlib import Path

import pytest

from condition_recommender.ingestion import (
    INTERMEDIATE_OBSERVATION_SCHEMA_VERSION,
    detect_adapter,
    preprocess_file,
    preprocess_files,
)
from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.conversion.input_schema import iter_conversion_records
from condition_recommender.models import IndexEligibility


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_records(path: str | Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def _literature_row() -> dict[str, str]:
    return {
        "reaction_id": "lit-1",
        "reaction_type": "Suzuki source label",
        "yield_pct": "81.5",
        "temperature_c": "80",
        "time_h": "2",
        "reaction_smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "reference": "Example Journal (2024), 1, 1-2",
        "reactant_cas": "",
        "product_cas": "",
        "reagent_cas": "584-08-7",
        "catalyst_cas": "14221-01-3",
        "solvent_cas": "108-88-3",
        "experimental_procedure": "Stirred at 80 C for 2 h.",
        "stages": "1",
        "steps": "1",
        "product_yield_1": "81.5",
        "product_yield_2": "",
        "product_yield_3": "",
        "product_yield_4": "",
        "product_yield_5": "",
        "product_yield_6": "",
        "product_yield_7": "",
        "notes": "source note",
    }


def _hitea_row() -> dict[str, str]:
    return {
        "Dataset entry number": "15408",
        "ReactionClass": "Pd coupling",
        "SCREEN_ID": "SCRN_11",
        "NOTEBOOK_ID": "NB-1",
        "REACTION_ID": "RXN-1",
        "KeyWord_STD": "SUZUKI",
        "ReactionGroup": "Metal Mediated",
        "RXN_SMILES": "[CH3:1][Br:2]>>[CH3:1][OH:2]",
        "PCAT_CMPD_ID": "P-1",
        "PRODUCT_STRUCTURE": "CO",
        "Product_Yield_PCT_Area_UV": "20",
        "Product_UV_Analysis_Wavelength_nm": "280",
        "Product_Yield_Mass_Ion_Count": "547137.877",
        "Product_Selectivity": "18.2",
        "Solvent_1_Name": "MeCN",
        "Reaction_T": "25",
        "Reaction_Time_hrs": "18",
        "Catalyst_1_eq": "0.0625",
        "Catalyst_1_ID[1]": "MFCD00012453",
        "catalyst_1_ID_1_SMILES": "CC(=O)[O-].[Pd+2]",
        "Catalyst_1_ID[2]": "",
        "catalyst_1_ID_2_SMILES": "",
        "Catalyst_1_Short_Hand": "Pd catalyst",
        "Catalyst_2_eq": "",
        "Catalyst_2_ID[1]": "",
        "catalyst_2_ID_1_SMILES": "",
        "Catalyst_2_ID[2]": "",
        "catalyst_2_ID_2_SMILES": "",
        "Catalyst_2_Short_Hand": "",
        "Reactant_1_eq": "1",
        "Reactant_1_ID": "MFCD1",
        "Reactant_1_SMILES": "CBr",
        "Reactant_2_eq": "",
        "Reactant_2_ID": "",
        "reactant_2_SMILES": "",
        "Reactant_3_eq": "",
        "Reactant_3_ID": "",
        "reactant_3_SMILES": "",
        "Reagent_1_eq": "3",
        "Reagent_1_ID": "MFCD00006493",
        "Reagent_1_Short_Hand": "Base",
        "Reagent_2_eq": "",
        "Reagent_2_ID": "",
        "Reagent_2_Short_Hand": "",
        "yie_react": "",
    }


def _weak_row(*, fg_a: str = "ArH", fg_b: str = "ArBr") -> dict[str, str]:
    return {
        "yield%": "48.01",
        "Base": "K2CO3",
        "Catalyst": "XantPhos Pd(allyl)Cl",
        "Solvent": "DMAc",
        "Ligand": "XantPhos",
        "Additive": "KOPiv",
        "Coupling Reagent": "",
        "Secondary Solvent": "water",
        "Tertiary Solvent": "",
        "Reaction Type": "CH-Activation",
        "FG A": fg_a,
        "FG B": fg_b,
        "z-Score": "0.42",
        "conditions": "17 h at 25 °C",
    }


def test_literature_preprocessing_is_source_faithful(tmp_path: Path) -> None:
    source = tmp_path / "literature.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_literature_row()])

    report = preprocess_file(source, output)
    record = _read_records(report["output_path"])[0]

    assert report["adapter_id"] == "literature_csv.v1"
    assert report["input_row_count"] == report["output_row_count"] == 1
    assert record["schema_version"] == INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
    assert record["observation_kind"] == "structure_backed"
    assert record["reaction"]["supplied_mapping_status"] == "not_supplied"
    assert record["reaction"]["source_reaction_type"] == "Suzuki source label"
    assert record["source"]["reference"].startswith("Example Journal")
    assert {
        item["source_role_hint"] for item in record["conditions"]["components"]
    } == {
        "catalyst",
        "reagent",
        "solvent",
    }
    assert record["conditions"]["stages"][0]["temperature_c"] == 80.0
    assert record["raw_fields"]["notes"] == "source note"


def test_hitea_preprocessing_groups_identifiers_and_preserves_outcome_basis(
    tmp_path: Path,
) -> None:
    source = tmp_path / "hitea.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_hitea_row()])

    report = preprocess_file(source, output)
    record = _read_records(report["output_path"])[0]

    assert report["adapter_id"] == "hitea_approved_csv.v1"
    assert record["reaction"]["supplied_mapping_status"] == "supplied_unvalidated"
    assert record["source"]["source_groups"]["screen_id"] == "SCRN_11"
    catalyst = next(
        item
        for item in record["conditions"]["components"]
        if item["component_key"] == "catalyst_1_component_1"
    )
    assert {item["identifier_type"] for item in catalyst["identifiers"]} == {
        "mfcd",
        "smiles",
    }
    group = record["conditions"]["component_groups"][0]
    assert group["display_name"] == "Pd catalyst"
    assert group["amount"]["value"] == 0.0625
    assert record["outcomes"][0]["outcome_type"] == "uv_area_yield_pct"
    assert record["outcomes"][0]["metadata"]["wavelength_nm"] == 280.0


def test_weak_label_preprocessing_retains_formerly_filtered_rows(
    tmp_path: Path,
) -> None:
    source = tmp_path / "weak.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_weak_row(fg_a="ArH", fg_b="ArH")])

    report = preprocess_file(source, output)
    record = _read_records(report["output_path"])[0]

    assert report["adapter_id"] == "weak_label_v2_1.v1"
    assert report["output_row_count"] == 1
    assert record["observation_kind"] == "label_only"
    assert record["reaction"]["reaction_smiles"] is None
    assert not record["reaction"]["structure_available"]
    assert "WEAK_LABEL_IDENTICAL_SITE_LABELS" in record["warnings"]
    assert record["conditions"]["stages"][0]["time_h"] == 17.0
    assert record["conditions"]["stages"][0]["temperature_c"] == 25.0
    assert record["outcomes"][1]["outcome_type"] == "source_z_score"


def test_preprocessing_cache_is_checksum_and_adapter_bound(tmp_path: Path) -> None:
    source = tmp_path / "literature.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_literature_row()])

    first = preprocess_file(source, output)
    first_bytes = Path(first["output_path"]).read_bytes()
    reused = preprocess_file(source, output)

    assert not first["reused"]
    assert reused["reused"]
    assert Path(reused["output_path"]).read_bytes() == first_bytes
    assert not (output / "literature.csv.preprocessing.json").exists()

    changed = _literature_row()
    changed["yield_pct"] = "82"
    _write_csv(source, [changed])
    rebuilt = preprocess_file(source, output)

    assert not rebuilt["reused"]
    assert rebuilt["source_sha256"] != first["source_sha256"]


def test_batch_rejects_colliding_source_filenames(tmp_path: Path) -> None:
    first = tmp_path / "one" / "source.csv"
    second = tmp_path / "two" / "source.csv"
    first.parent.mkdir()
    second.parent.mkdir()
    _write_csv(first, [_literature_row()])
    _write_csv(second, [_literature_row()])

    with pytest.raises(ValueError, match="colliding names"):
        preprocess_files((first, second), tmp_path / "output")


def test_adapter_detection_rejects_unknown_schema(tmp_path: Path) -> None:
    source = tmp_path / "unknown.csv"
    _write_csv(source, [{"unknown": "value"}])

    with pytest.raises(ValueError, match="No registered source adapter"):
        detect_adapter(source)


def test_generic_converter_consumes_structure_backed_intermediate_file(
    tmp_path: Path,
) -> None:
    source = tmp_path / "literature.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_literature_row()])
    preprocessing = preprocess_file(source, output)

    raw_record = next(iter_conversion_records(preprocessing["output_path"]))
    converted = convert_record(raw_record)

    assert raw_record.primary_outcome_type == "reported_yield_pct"
    assert len(raw_record.condition_component_inputs) == 3
    assert converted.reaction_signature is not None
    assert converted.resolved_recipe_id.startswith("RCR1:")
    assert converted.source["input_schema_version"] == (
        INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
    )
    assert converted.source["primary_outcome_type"] == "reported_yield_pct"


def test_generic_converter_retains_but_never_indexes_label_only_observation(
    tmp_path: Path,
) -> None:
    source = tmp_path / "weak.csv"
    output = tmp_path / "intermediate"
    _write_csv(source, [_weak_row()])
    preprocessing = preprocess_file(source, output)

    raw_record = next(iter_conversion_records(preprocessing["output_path"]))
    converted = convert_record(raw_record)

    assert raw_record.reaction_smiles == ""
    assert raw_record.condition_component_inputs
    assert converted.reaction_signature is None
    assert converted.canonical_reaction_id is None
    assert converted.index_eligibility == IndexEligibility.INELIGIBLE
    assert converted.resolved_recipe is not None


def test_conversion_engine_accepts_intermediate_artifact(tmp_path: Path) -> None:
    source = tmp_path / "literature.csv"
    intermediate = tmp_path / "intermediate"
    converted = tmp_path / "converted"
    _write_csv(source, [_literature_row()])
    preprocessing = preprocess_file(source, intermediate)

    report = convert_datasets(preprocessing["output_path"], converted)

    assert report["file_count"] == 1
    assert report["row_count"] == 1
    assert report["source_row_counts"] == {"literature_reaction_dataset": 1}
    with (converted / "records.jsonl").open(encoding="utf-8") as handle:
        record = json.loads(next(handle))
    assert record["source"]["input_schema_version"] == (
        INTERMEDIATE_OBSERVATION_SCHEMA_VERSION
    )
