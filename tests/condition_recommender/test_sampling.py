import csv
import json
from pathlib import Path

from condition_recommender.conversion.input_schema import iter_csv_records
from condition_recommender.conversion.sampling import build_reference_safe_samples


_FIELDS = (
    "reaction_id",
    "reaction_type",
    "yield_pct",
    "reaction_smiles",
    "reference",
    "reagent_cas",
    "catalyst_cas",
    "solvent_cas",
    "experimental_procedure",
    "stages",
)


def _write_source(path: Path, *, offset: int, count: int) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS)
        writer.writeheader()
        for local_index in range(count):
            index = offset + local_index
            reference = (
                "Journal of Shared Conditions (2022), 10, 1-10"
                if local_index in {0, 1}
                else f"Journal {index} (2020), {index}, 1-2"
            )
            reaction = (
                "CCBr.N>>CCN"
                if local_index == 2
                else f"[CH3:{index + 1}]Br>>[CH3:{index + 1}]"
            )
            writer.writerow(
                {
                    "reaction_id": f"reaction-{index}",
                    "reaction_type": "source-label",
                    "yield_pct": str(index % 101),
                    "reaction_smiles": reaction,
                    "reference": reference,
                    "reagent_cas": "584-08-7",
                    "catalyst_cas": "14221-01-3" if index % 2 else "",
                    "solvent_cas": "108-88-3",
                    "experimental_procedure": "procedure" if index % 3 == 0 else "",
                    "stages": "2" if index % 4 == 0 else "1",
                }
            )


def _build_source_directory(tmp_path: Path) -> Path:
    source = tmp_path / "reaction_dataset"
    source.mkdir()
    for file_index in range(3):
        _write_source(
            source / f"source_{file_index}.csv",
            offset=file_index * 40,
            count=40,
        )
    return source


def test_reference_safe_sampling_is_deterministic_and_leakage_free(
    tmp_path: Path,
) -> None:
    source = _build_source_directory(tmp_path)
    first = tmp_path / "first"
    second = tmp_path / "second"
    parameters = {
        "smoke_size": 6,
        "development_size": 24,
        "validation_size": 12,
        "test_size": 12,
        "seed": 23,
        "max_rows_per_group": 4,
    }

    first_report = build_reference_safe_samples(source, first, **parameters)
    second_report = build_reference_safe_samples(source, second, **parameters)

    assert first_report == second_report
    assert first_report["actual_sizes"] == {
        "smoke": 6,
        "development": 24,
        "validation": 12,
        "untouched_test": 12,
    }
    assert first_report["reference_leakage_count"] == 0
    assert first_report["reaction_text_leakage_count"] == 0
    assert (first / "sample_manifest.v1.json").read_bytes() == (
        second / "sample_manifest.v1.json"
    ).read_bytes()
    for name in ("smoke", "development", "validation", "untouched_test"):
        assert (first / f"{name}.csv").read_bytes() == (
            second / f"{name}.csv"
        ).read_bytes()


def test_smoke_is_a_development_subset_and_source_provenance_survives(
    tmp_path: Path,
) -> None:
    source = _build_source_directory(tmp_path)
    output = tmp_path / "samples"
    build_reference_safe_samples(
        source,
        output,
        smoke_size=5,
        development_size=20,
        validation_size=10,
        test_size=10,
        seed=17,
        max_rows_per_group=4,
    )
    manifest = json.loads(
        (output / "sample_manifest.v1.json").read_text(encoding="utf-8")
    )
    development = {
        (entry["source_path"], entry["source_row_number"])
        for entry in manifest["entries"]
        if entry["partition"] == "development"
    }
    smoke = {
        (entry["source_path"], entry["source_row_number"])
        for entry in manifest["entries"]
        if entry["smoke"]
    }

    assert smoke
    assert smoke.issubset(development)
    sampled_records = tuple(iter_csv_records(output / "development.csv"))
    assert sampled_records
    assert {record.source_dataset for record in sampled_records}.issubset(
        {"source_0", "source_1", "source_2"}
    )
    assert all(
        Path(record.source_path).parent == source for record in sampled_records
    )
    assert all(record.source_row_number >= 2 for record in sampled_records)


def test_sampling_rejects_smoke_larger_than_development(tmp_path: Path) -> None:
    source = _build_source_directory(tmp_path)

    try:
        build_reference_safe_samples(
            source,
            tmp_path / "samples",
            smoke_size=11,
            development_size=10,
            validation_size=5,
            test_size=5,
        )
    except ValueError as exc:
        assert "subset of development" in str(exc)
    else:
        raise AssertionError("Expected invalid sample targets to fail")
