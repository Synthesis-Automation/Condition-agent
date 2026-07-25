"""Clean source reaction-label CSV files using reactive-taxonomy labels."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable

from reactive_taxonomy import resolve_source_label, validate_source_label_mappings


SOURCE_TO_OUTPUT = {
    "yield%": "yield_pct",
    "Base": "base",
    "Catalyst": "catalyst",
    "Solvent": "primary_solvent",
    "Ligand": "ligand",
    "Additive": "additive",
    "Coupling Reagent": "coupling_reagent",
    "Secondary Solvent": "secondary_solvent",
    "Tertiary Solvent": "tertiary_solvent",
    "Reaction Type": "source_reaction_type",
    "z-Score": "z_score",
    "conditions": "procedure_text",
}

OUTPUT_SITE_SUFFIXES = (
    "display_label",
    "signature",
    "center_class",
    "attachment_class",
    "alpha_branched",
)

OUTPUT_FIELDNAMES = (
    "source_reaction_type",
    *(f"reactive_site_1_{suffix}" for suffix in OUTPUT_SITE_SUFFIXES),
    *(f"reactive_site_2_{suffix}" for suffix in OUTPUT_SITE_SUFFIXES),
    "yield_pct",
    "z_score",
    "catalyst",
    "ligand",
    "base",
    "primary_solvent",
    "secondary_solvent",
    "tertiary_solvent",
    "additive",
    "coupling_reagent",
    "procedure_text",
)


def clean_rows(
    rows: Iterable[dict[str, str]],
) -> tuple[list[dict[str, str]], Counter[str]]:
    """Filter requested rows and attach normalized FG fields."""
    cleaned: list[dict[str, str]] = []
    stats: Counter[str] = Counter()

    for row in rows:
        stats["input_rows"] += 1
        fg_a = (row.get("FG A") or "").strip()
        fg_b = (row.get("FG B") or "").strip()

        both_blank = not fg_a and not fg_b
        identical = fg_a == fg_b
        protecting_group = "Protecting Group" in {fg_a, fg_b}

        if both_blank:
            stats["matched_both_blank"] += 1
        if identical:
            stats["matched_identical"] += 1
        if protecting_group:
            stats["matched_protecting_group"] += 1
        if both_blank or identical or protecting_group:
            stats["removed_union"] += 1
            continue

        output = {
            destination: (row.get(source) or "").strip()
            for source, destination in SOURCE_TO_OUTPUT.items()
        }
        for source_column, prefix in (
            ("FG A", "reactive_site_1"),
            ("FG B", "reactive_site_2"),
        ):
            source = (row.get(source_column) or "").strip()
            mapping = resolve_source_label(source)
            mapping_columns = mapping.to_columns(prefix)
            output.update(
                {
                    f"{prefix}_{suffix}": mapping_columns[f"{prefix}_{suffix}"]
                    for suffix in OUTPUT_SITE_SUFFIXES
                }
            )
            stats[f"mapping_status:{mapping.mapping_status}"] += 1
            if mapping.base_label != source:
                stats[f"replaced:{source}->{mapping.base_label}"] += 1

        cleaned.append(output)
        stats["output_rows"] += 1

    return cleaned, stats


def clean_csv(source: Path, destination: Path) -> Counter[str]:
    """Clean ``source`` and write a new CSV to ``destination``."""
    validation_errors = validate_source_label_mappings()
    if validation_errors:
        raise ValueError(f"Invalid source-label mappings: {validation_errors}")

    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {source}")
        rows, stats = clean_rows(dict(row) for row in reader)
        required = set(SOURCE_TO_OUTPUT) | {"FG A", "FG B"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing source columns: {missing}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDNAMES,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()

    stats = clean_csv(args.source, args.destination)
    for key in sorted(stats):
        print(f"{key}: {stats[key]}")


if __name__ == "__main__":
    main()
