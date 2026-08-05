import csv
from pathlib import Path

from condition_registry import ConditionRegistry
from condition_registry import curation
from condition_registry import weak_label_reconciliation as reconciliation


SOURCE_FIELDS = (
    "Base",
    "Catalyst",
    "Solvent",
    "Ligand",
    "Additive",
    "Coupling Reagent",
    "Secondary Solvent",
    "Tertiary Solvent",
)


def _write_csv(
    path: Path,
    fieldnames: tuple[str, ...],
    rows: tuple[dict[str, str], ...] = (),
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _registry_paths(tmp_path: Path) -> tuple[Path, Path, Path]:
    substances = tmp_path / "substances.csv"
    additions = tmp_path / "additions.csv"
    identifiers = tmp_path / "identifiers.csv"
    _write_csv(
        substances,
        curation.SUBSTANCE_FIELDNAMES,
        (
            {
                "name": "Triethylamine",
                "cas": "121-44-8",
                "role_1": "base",
            },
            {
                "name": "dppf-PdCl2",
                "cas": "72287-26-4",
                "role_1": "metal_catalyst",
            },
        ),
    )
    _write_csv(additions, curation.ADDITION_FIELDNAMES)
    _write_csv(identifiers, curation.IDENTIFIER_FIELDNAMES)
    return substances, additions, identifiers


def test_split_component_value_preserves_systematic_name_commas() -> None:
    assert reconciliation._split_component_value(
        "1,10-Phenanthroline, NEt3, Missing"
    ) == ("1,10-Phenanthroline", "NEt3")


def test_extract_weak_label_reagents_deduplicates_row_mentions(
    tmp_path,
) -> None:
    source = tmp_path / "weak.csv"
    _write_csv(
        source,
        SOURCE_FIELDS,
        (
            {
                "Base": "K2CO3, NEt3",
                "Ligand": "Xantphos, XantPhos",
                "Solvent": "1,10-Phenanthroline",
            },
            {"Base": "NEt3", "Ligand": "No ligand"},
        ),
    )

    reagents = reconciliation.extract_weak_label_reagents(source)
    by_name = {item.normalized_name: item for item in reagents}

    assert len(reagents) == 4
    assert by_name["net3"].role_counts == {"base": 2}
    assert by_name["xantphos"].role_counts == {"ligand": 1}
    assert set(by_name["xantphos"].aliases) == {"Xantphos", "XantPhos"}
    assert "1 10 phenanthroline" in by_name


def test_run_reconciliation_applies_compact_and_curated_aliases(
    tmp_path,
) -> None:
    source = tmp_path / "weak.csv"
    output = tmp_path / "output"
    substances, additions, identifiers = _registry_paths(tmp_path)
    _write_csv(
        source,
        SOURCE_FIELDS,
        (
            {
                "Base": "NEt3",
                "Catalyst": "dppfPdCl2",
                "Additive": "Unknown reagent",
            },
        ),
    )

    summary = reconciliation.run_weak_label_reconciliation(
        source,
        output,
        apply_aliases=True,
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )

    assert summary["unique_reagents"] == 3
    assert summary["matched_reagents"] == 2
    assert summary["unresolved_reagents"] == 1
    assert summary["aliases_added"] == 2
    registry = ConditionRegistry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    assert registry.resolve(name="NEt3").substance.substance_id == "cas:121-44-8"
    assert (
        registry.resolve(name="dppfPdCl2").substance.substance_id
        == "cas:72287-26-4"
    )
    with (output / "weak_label_unresolved_reagents.csv").open(
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        unresolved = list(csv.DictReader(handle))
    assert unresolved[0]["name"] == "Unknown reagent"

