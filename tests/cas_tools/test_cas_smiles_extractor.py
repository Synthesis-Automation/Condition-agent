import csv
import json
from pathlib import Path

from cas_tools.cas_smiles_extractor import (
    CASSmilesPair,
    discover_csv_files,
    extract_cas_smiles_pairs_from_csv,
    extract_cas_smiles_pairs_from_folder,
    write_cas_smiles_pairs,
)


def _write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_extracts_pairs_recursively_from_all_json_columns(tmp_path: Path) -> None:
    source = tmp_path / "literature.csv"
    fieldnames = [
        "reactant_cas",
        "reactant_smiles",
        "reactants_json",
        "products_json",
        "reagents_json",
        "other_json",
    ]
    _write_rows(
        source,
        fieldnames,
        [
            {
                "reactant_cas": "64-17-5",
                "reactant_smiles": "this flat value must not override JSON",
                "reactants_json": json.dumps(
                    [{"index": 1, "cas_rn": "64-17-5", "smiles": "CCO"}]
                ),
                "products_json": json.dumps(
                    [
                        {
                            "cas_rn": ["67-56-1", "123-39-7"],
                            "smiles": "CO",
                        }
                    ]
                ),
                "reagents_json": json.dumps(
                    [{"cas_number": "67-64-1", "canonical_smiles": "CC(C)=O"}]
                ),
                "other_json": json.dumps(
                    {"wrapper": {"cas_no": "71-43-2", "compound_smiles": "c1ccccc1"}}
                ),
            }
        ],
    )

    result = extract_cas_smiles_pairs_from_csv(source)

    assert set(result.pairs) == {
        CASSmilesPair("64-17-5", "CCO"),
        CASSmilesPair("67-56-1", "CO"),
        CASSmilesPair("123-39-7", "CO"),
        CASSmilesPair("67-64-1", "CC(C)=O"),
        CASSmilesPair("71-43-2", "c1ccccc1"),
    }
    assert result.rows_read == 1
    assert not result.warnings


def test_preserves_conflicts_and_deduplicates_exact_pairs(tmp_path: Path) -> None:
    source = tmp_path / "conflicts.csv"
    _write_rows(
        source,
        ["items_json"],
        [
            {
                "items_json": json.dumps(
                    [
                        {"cas_rn": "64-17-5", "smiles": "CCO"},
                        {"cas_rn": "64-17-5", "smiles": "OCC"},
                        {"cas_rn": "64-17-5", "smiles": "CCO"},
                        {"cas_rn": "not-a-cas", "smiles": "N"},
                    ]
                )
            }
        ],
    )

    result = extract_cas_smiles_pairs_from_csv(source)

    assert result.pair_occurrences == 3
    assert result.pairs == (
        CASSmilesPair("64-17-5", "CCO"),
        CASSmilesPair("64-17-5", "OCC"),
    )


def test_carries_reaction_id_and_citation_and_preserves_each_reaction(tmp_path: Path) -> None:
    source = tmp_path / "provenance.csv"
    _write_rows(
        source,
        ["reaction_id", "citation", "reactants_json"],
        [
            {
                "reaction_id": "RXN-1",
                "citation": "Journal A (2024), 1, 10-12",
                "reactants_json": json.dumps(
                    [{"cas_rn": "64-17-5", "smiles": "CCO"}]
                ),
            },
            {
                "reaction_id": "RXN-2",
                "citation": "Journal B (2025), 2, 20-22",
                "reactants_json": json.dumps(
                    [{"cas_rn": "64-17-5", "smiles": "CCO"}]
                ),
            },
        ],
    )

    result = extract_cas_smiles_pairs_from_csv(source)

    assert result.pairs == (
        CASSmilesPair("64-17-5", "CCO", "RXN-1", "Journal A (2024), 1, 10-12"),
        CASSmilesPair("64-17-5", "CCO", "RXN-2", "Journal B (2025), 2, 20-22"),
    )


def test_uses_matching_flat_columns_only_without_structured_role(tmp_path: Path) -> None:
    source = tmp_path / "flat.csv"
    _write_rows(
        source,
        ["compound_cas", "compound_smiles"],
        [
            {"compound_cas": "64-17-5", "compound_smiles": "CCO"},
            {
                "compound_cas": "67-56-1, 67-64-1",
                "compound_smiles": "CO.CC(C)=O",
            },
        ],
    )

    result = extract_cas_smiles_pairs_from_csv(source)

    assert set(result.pairs) == {
        CASSmilesPair("64-17-5", "CCO"),
        CASSmilesPair("67-56-1", "CO"),
        CASSmilesPair("67-64-1", "CC(C)=O"),
    }


def test_malformed_json_warns_without_scanning_unrelated_text(tmp_path: Path) -> None:
    source = tmp_path / "malformed.csv"
    _write_rows(
        source,
        ["reactants_json", "notes"],
        [{"reactants_json": '[{"cas_rn": "64-17-5"', "notes": "64-17-5 and CCO"}],
    )

    result = extract_cas_smiles_pairs_from_csv(source)

    assert result.pairs == ()
    assert len(result.warnings) == 1
    assert "invalid JSON" in result.warnings[0]


def test_folder_discovery_excludes_output_and_writer_has_exact_columns(tmp_path: Path) -> None:
    source = tmp_path / "nested" / "source.csv"
    output = tmp_path / "pairs.csv"
    _write_rows(
        source,
        ["items_json"],
        [{"items_json": json.dumps({"cas": "64-17-5", "smiles": "CCO"})}],
    )
    output.write_text("old output", encoding="utf-8")
    (tmp_path / "ignored.txt").write_text("64-17-5", encoding="utf-8")

    assert discover_csv_files(tmp_path, excluded_paths=[output]) == [source]
    result = extract_cas_smiles_pairs_from_folder(tmp_path, excluded_paths=[output])
    count = write_cas_smiles_pairs(result.pairs, output)

    assert count == 1
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows == [
        ["cas_no", "compound_smiles", "reaction_id", "citation"],
        ["64-17-5", "CCO", "", ""],
    ]


def test_writer_emits_one_row_per_cas_and_selects_best_supported_smiles(
    tmp_path: Path,
) -> None:
    output = tmp_path / "deduplicated.csv"
    pairs = [
        CASSmilesPair("64-17-5", "OCC", "RXN-3", "Citation C"),
        CASSmilesPair("64-17-5", "CCO", "RXN-2", "Citation B"),
        CASSmilesPair("64-17-5", "CCO", "RXN-1", "Citation A"),
        CASSmilesPair("67-56-1", "CO", "RXN-4", "Citation D"),
    ]

    count = write_cas_smiles_pairs(pairs, output)

    assert count == 2
    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows == [
        ["cas_no", "compound_smiles", "reaction_id", "citation"],
        ["64-17-5", "CCO", "RXN-1", "Citation A"],
        ["67-56-1", "CO", "RXN-4", "Citation D"],
    ]
