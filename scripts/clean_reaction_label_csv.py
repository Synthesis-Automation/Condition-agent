"""Clean source reaction-label CSV files using reactive-taxonomy labels."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable

from reactive_taxonomy import resolve_source_label, validate_source_label_mappings


FG_COLUMNS = (
    "Source",
    "Display",
    "Signature",
    "Center Class",
    "Attachment Class",
    "Alpha Branched",
    "Qualifier Scope",
    "Mapping Status",
)


def _output_fieldnames(source_fieldnames: list[str]) -> list[str]:
    fieldnames: list[str] = []
    for fieldname in source_fieldnames:
        if fieldname in {"FG A", "FG B"}:
            fieldnames.append(f"{fieldname} Source")
            fieldnames.append(fieldname)
            fieldnames.extend(
                f"{fieldname} {suffix}" for suffix in FG_COLUMNS if suffix != "Source"
            )
        else:
            fieldnames.append(fieldname)
    return fieldnames


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

        for column in ("FG A", "FG B"):
            source = (row.get(column) or "").strip()
            mapping = resolve_source_label(source)
            row.update(mapping.to_columns(column))
            stats[f"mapping_status:{mapping.mapping_status}"] += 1
            if mapping.base_label != source:
                stats[f"replaced:{source}->{mapping.base_label}"] += 1

        cleaned.append(row)
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
        fieldnames = _output_fieldnames(list(reader.fieldnames))

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
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
