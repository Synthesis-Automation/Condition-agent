import csv
import json

from condition_registry.coverage import analyze_coverage


def test_dataset_coverage_routes_identity_and_role_categories(tmp_path) -> None:
    source = tmp_path / "source.csv"
    rows = [{
        "reaction_id": "r1",
        "catalyst_cas": "14221-01-3, 584-08-7",
        "reagent_cas": "496-15-1, 999999-99-4",
        "solvent_cas": "7732-18-5, 7732-18-4",
    }]
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = analyze_coverage({"test_family": source}, tmp_path / "out")
    assert report["source_rows"] == 1
    assert report["unique_identifier_source_pairs"] == 6
    assert report["category_counts"] == {
        "resolved": 1,
        "multiple_roles": 1,
        "generic_only": 1,
        "role_conflict": 1,
        "missing_substance": 1,
        "invalid_identifier": 1,
    }
    assert json.loads((tmp_path / "out" / "coverage_report.json").read_text()) == report
    for filename in (
        "resolved.csv", "ambiguous_roles.csv", "generic_only.csv",
        "role_conflicts.csv", "missing_substances.csv", "invalid_identifiers.csv",
    ):
        assert (tmp_path / "out" / filename).exists()
