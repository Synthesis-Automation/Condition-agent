import csv
import json

from condition_recommender.conversion.audit import audit_datasets
from condition_recommender.conversion.identities import (
    canonical_reaction_identity,
    observation_id,
    raw_recipe_id,
)
from condition_recommender.conversion.input_schema import adapt_row, iter_csv_records


def _row(
    reaction_id: str,
    reaction_smiles: str,
    *,
    reaction_type: str = "source-label",
    yield_pct: str = "80",
    catalyst_cas: str = "14221-01-3",
    reagent_cas: str = "584-08-7",
    solvent_cas: str = "108-88-3",
) -> dict[str, str]:
    return {
        "reaction_id": reaction_id,
        "reaction_type": reaction_type,
        "reaction_smiles": reaction_smiles,
        "yield_pct": yield_pct,
        "temperature_c": "",
        "time_h": "",
        "reference": "reference",
        "reactant_cas": "",
        "product_cas": "",
        "reagent_cas": reagent_cas,
        "catalyst_cas": catalyst_cas,
        "solvent_cas": solvent_cas,
        "experimental_procedure": "",
        "stages": "1",
        "steps": "1",
        "notes": "",
    }


def _write(path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_common_adapter_preserves_source_and_supports_aliases() -> None:
    record = adapt_row(
        {
            "ReactionID": "source-1",
            "ReactionType": "Suzuki",
            "ReactionSmiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            "Yield": "81.5",
            "reagents_cas": "584-08-7; 584-08-7",
            "catalysts_cas": "14221-01-3",
            "solvents_cas": "108-88-3",
            "custom_field": "preserved",
        },
        source_dataset="pilot",
        source_path="pilot.csv",
        source_row_number=2,
    )

    assert record.reaction_id == "source-1"
    assert record.source_declared_family == "Suzuki"
    assert record.yield_pct == 81.5
    assert record.reagent_cas == ("584-08-7",)
    assert record.raw_fields["custom_field"] == "preserved"
    assert observation_id(record).startswith("OBS1:")
    assert raw_recipe_id(record).startswith("RAWCOND1:")


def test_canonical_reaction_identity_is_component_and_agent_invariant() -> None:
    first = canonical_reaction_identity(
        "Brc1ccccc1.OB(O)c1ccccc1>CCN>c1ccc(-c2ccccc2)cc1"
    )
    reordered = canonical_reaction_identity(
        "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert first is not None
    assert reordered is not None
    assert first.reaction_id == reordered.reaction_id
    assert canonical_reaction_identity("not-a-reaction") is None


def test_streaming_audit_reports_duplicates_without_modifying_sources(tmp_path) -> None:
    dataset_dir = tmp_path / "datasets"
    dataset_dir.mkdir()
    suzuki = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    suzuki_reordered = "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    first_path = dataset_dir / "first.csv"
    second_path = dataset_dir / "second.csv"
    _write(
        first_path,
        [
            _row("a-1", suzuki),
            _row("a-2", suzuki_reordered, solvent_cas="68-12-2"),
        ],
    )
    _write(
        second_path,
        [
            _row("b-1", suzuki),
            _row(
                "b-2",
                "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
                reaction_type="hydrogenation",
            ),
            _row("b-3", "not-a-reaction", yield_pct="not-numeric"),
        ],
    )
    source_bytes = {
        path.name: path.read_bytes() for path in (first_path, second_path)
    }

    output_dir = tmp_path / "audit"
    report = audit_datasets(
        dataset_dir, output_dir, chemistry_sample_per_file=10
    )

    assert report["file_count"] == 2
    assert report["row_count"] == 5
    assert report["canonical_reaction_count"] == 4
    duplicates = report["duplicate_groups"]
    assert duplicates["unique_canonical_reactions"] == 2
    assert duplicates["repeated_reaction_groups"] == 1
    assert duplicates["multi_recipe_reaction_groups"] == 1
    assert duplicates["cross_dataset_reaction_groups"] == 1
    assert report["chemistry_sample"]["sampled_rows"] == 5
    assert report["chemistry_sample"]["signature_count"] == 4
    assert report["yield_counts"] == {"missing_or_non_numeric": 1, "valid": 4}
    assert json.loads((output_dir / "audit_report.json").read_text()) == report
    assert (output_dir / "audit_report.md").exists()
    assert {
        path.name: path.read_bytes() for path in (first_path, second_path)
    } == source_bytes


def test_iter_csv_records_is_streaming_and_source_faithful(tmp_path) -> None:
    path = tmp_path / "pilot.csv"
    _write(path, [_row("row-1", "CCBr.CN>>CCNC")])

    iterator = iter_csv_records(path)
    record = next(iterator)

    assert record.source_dataset == "pilot"
    assert record.source_row_number == 2
    assert record.raw_fields["reaction_id"] == "row-1"
