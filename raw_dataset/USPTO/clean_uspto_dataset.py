"""Create a chemistry-first, condition-bearing subset of the USPTO CSV.

The source CSV is preserved.  Rows are retained only when their reaction
structures and supplied atom mapping are internally consistent, their mapper
confidence meets the configured threshold, and they contain at least one
parseable catalyst or reagent.  Source split and reaction-class labels are
excluded from the cleaned artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles


SOURCE_COLUMNS = (
    "source",
    "canonical_rxn",
    "catalyst1",
    "solvent1",
    "solvent2",
    "reagent1",
    "reagent2",
    "dataset",
    "rxn_category",
    "rxn_class_name",
    "remapped_rxn",
    "confidence",
)
OUTPUT_COLUMNS = (
    "source",
    "canonical_rxn",
    "catalyst1",
    "solvent1",
    "solvent2",
    "reagent1",
    "reagent2",
    "remapped_rxn",
    "confidence",
)
DROPPED_COLUMNS = ("dataset", "rxn_category", "rxn_class_name")
CONDITION_COLUMNS = (
    "catalyst1",
    "solvent1",
    "solvent2",
    "reagent1",
    "reagent2",
)
ACTIVE_CONDITION_COLUMNS = ("catalyst1", "reagent1", "reagent2")
DEFAULT_MINIMUM_CONFIDENCE = 0.5
DEFAULT_BATCH_SIZE = 512


@dataclass(frozen=True)
class RowValidation:
    """Validation result for one source row."""

    accepted: bool
    reasons: tuple[str, ...] = ()


@dataclass(frozen=True)
class _MappedSide:
    atoms: dict[int, str]
    bonds: dict[tuple[int, int], str]
    hydrogen_counts: dict[int, int]
    formal_charges: dict[int, int]
    heavy_atom_count: int
    mapped_heavy_atom_count: int
    duplicate_map: bool


def _split_reaction(reaction_smiles: str) -> tuple[str, str] | None:
    text = str(reaction_smiles or "").strip()
    if text.count(">>") == 1:
        left, right = text.split(">>", 1)
    else:
        parts = text.split(">")
        if len(parts) != 3:
            return None
        left, _, right = parts
    if not left.strip() or not right.strip():
        return None
    return left, right


def _parse_components(side: str) -> list[object] | None:
    components = [token.strip() for token in side.split(".") if token.strip()]
    if not components:
        return None
    molecules = [parse_smiles(component) for component in components]
    if any(molecule is None for molecule in molecules):
        return None
    return molecules


def _mapped_side(molecules: Sequence[object]) -> _MappedSide:
    atoms: dict[int, str] = {}
    bonds: dict[tuple[int, int], str] = {}
    hydrogen_counts: dict[int, int] = {}
    formal_charges: dict[int, int] = {}
    heavy_atom_count = 0
    mapped_heavy_atom_count = 0
    duplicate_map = False

    for molecule in molecules:
        for atom in molecule.GetAtoms():
            is_heavy = int(atom.GetAtomicNum()) > 1
            if is_heavy:
                heavy_atom_count += 1
            map_number = int(atom.GetAtomMapNum())
            if not map_number:
                continue
            if is_heavy:
                mapped_heavy_atom_count += 1
            if map_number in atoms:
                duplicate_map = True
                continue
            atoms[map_number] = str(atom.GetSymbol())
            hydrogen_counts[map_number] = int(
                atom.GetTotalNumHs(includeNeighbors=True)
            )
            formal_charges[map_number] = int(atom.GetFormalCharge())

        for bond in molecule.GetBonds():
            left_map = int(bond.GetBeginAtom().GetAtomMapNum())
            right_map = int(bond.GetEndAtom().GetAtomMapNum())
            if not left_map or not right_map:
                continue
            key = tuple(sorted((left_map, right_map)))
            order = str(bond.GetBondType()).upper()
            if key in bonds and bonds[key] != order:
                duplicate_map = True
            bonds[key] = order

    return _MappedSide(
        atoms=atoms,
        bonds=bonds,
        hydrogen_counts=hydrogen_counts,
        formal_charges=formal_charges,
        heavy_atom_count=heavy_atom_count,
        mapped_heavy_atom_count=mapped_heavy_atom_count,
        duplicate_map=duplicate_map,
    )


def _has_mapped_change(left: _MappedSide, right: _MappedSide) -> bool:
    shared_maps = set(left.atoms).intersection(right.atoms)
    if any(
        left.bonds.get(key) != right.bonds.get(key)
        for key in set(left.bonds).union(right.bonds)
    ):
        return True
    return any(
        left.hydrogen_counts[map_number] != right.hydrogen_counts[map_number]
        or left.formal_charges[map_number] != right.formal_charges[map_number]
        for map_number in shared_maps
    )


def _validate_reaction(reaction_smiles: str, *, mapped: bool) -> tuple[str, ...]:
    split = _split_reaction(reaction_smiles)
    if split is None:
        return ("invalid_reaction_format",)
    left_molecules = _parse_components(split[0])
    right_molecules = _parse_components(split[1])
    if left_molecules is None or right_molecules is None:
        return ("unparseable_reaction_structure",)
    if not mapped:
        return ()

    left = _mapped_side(left_molecules)
    right = _mapped_side(right_molecules)
    reasons: list[str] = []
    if left.duplicate_map or right.duplicate_map:
        reasons.append("duplicate_or_conflicting_atom_map")
    if not left.atoms or not right.atoms:
        reasons.append("missing_atom_mapping")
    if right.mapped_heavy_atom_count != right.heavy_atom_count:
        reasons.append("incomplete_product_atom_mapping")
    if set(right.atoms) - set(left.atoms):
        reasons.append("product_map_without_reactant_source")
    if any(
        left.atoms[map_number] != right.atoms[map_number]
        for map_number in set(left.atoms).intersection(right.atoms)
    ):
        reasons.append("atom_map_element_mismatch")
    if not reasons and not _has_mapped_change(left, right):
        reasons.append("no_mapped_transformation")
    return tuple(reasons)


def validate_row(
    row: dict[str, str],
    *,
    minimum_confidence: float = DEFAULT_MINIMUM_CONFIDENCE,
) -> RowValidation:
    """Return deterministic rejection reasons for one USPTO source row."""
    reasons: list[str] = []
    if not any((row.get(column) or "").strip() for column in ACTIVE_CONDITION_COLUMNS):
        reasons.append("missing_active_condition")
    if any(
        value and parse_smiles(value) is None
        for column in CONDITION_COLUMNS
        for value in ((row.get(column) or "").strip(),)
    ):
        reasons.append("unparseable_condition")

    try:
        confidence = float((row.get("confidence") or "").strip())
    except ValueError:
        confidence = math.nan
    if not math.isfinite(confidence) or not 0.0 <= confidence <= 1.0:
        reasons.append("invalid_mapping_confidence")
    elif confidence < minimum_confidence:
        reasons.append("low_mapping_confidence")

    reasons.extend(
        f"canonical_{reason}"
        for reason in _validate_reaction(row.get("canonical_rxn", ""), mapped=False)
    )
    reasons.extend(
        f"remapped_{reason}"
        for reason in _validate_reaction(row.get("remapped_rxn", ""), mapped=True)
    )
    unique_reasons = tuple(dict.fromkeys(reasons))
    return RowValidation(not unique_reasons, unique_reasons)


def _validate_batch(
    payload: tuple[list[dict[str, str]], float],
) -> tuple[list[dict[str, str]], list[RowValidation]]:
    rows, minimum_confidence = payload
    return rows, [
        validate_row(row, minimum_confidence=minimum_confidence)
        for row in rows
    ]


def _batches(
    rows: Iterable[dict[str, str]],
    *,
    batch_size: int,
    minimum_confidence: float,
) -> Iterator[tuple[list[dict[str, str]], float]]:
    batch: list[dict[str, str]] = []
    for row in rows:
        batch.append(row)
        if len(batch) == batch_size:
            yield batch, minimum_confidence
            batch = []
    if batch:
        yield batch, minimum_confidence


def _cleaned_row(row: dict[str, str]) -> dict[str, str]:
    return {column: (row.get(column) or "").strip() for column in OUTPUT_COLUMNS}


def _row_digest(row: dict[str, str]) -> bytes:
    encoded = json.dumps(
        [row[column] for column in OUTPUT_COLUMNS],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).digest()


def clean_dataset(
    source_path: Path,
    output_path: Path,
    report_path: Path,
    *,
    minimum_confidence: float = DEFAULT_MINIMUM_CONFIDENCE,
    processes: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> dict[str, object]:
    """Filter the USPTO CSV and write a reproducible cleanup report."""
    if not 0.0 <= minimum_confidence <= 1.0:
        raise ValueError("minimum_confidence must be in [0, 1]")
    if source_path.resolve() == output_path.resolve():
        raise ValueError("output_path must not overwrite the source dataset")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")
    reason_counts: Counter[str] = Counter()
    seen_rows: set[bytes] = set()
    source_rows = 0
    accepted_rows = 0
    duplicate_rows = 0

    with source_path.open("r", encoding="utf-8-sig", newline="") as source_handle:
        reader = csv.DictReader(source_handle)
        if tuple(reader.fieldnames or ()) != SOURCE_COLUMNS:
            raise ValueError(
                "Unexpected source schema: " + repr(tuple(reader.fieldnames or ()))
            )
        with temporary_output.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=list(OUTPUT_COLUMNS),
                lineterminator="\n",
            )
            writer.writeheader()
            payloads = _batches(
                reader,
                batch_size=batch_size,
                minimum_confidence=minimum_confidence,
            )
            if processes == 1:
                validated_batches = map(_validate_batch, payloads)
                pool = None
            else:
                pool = multiprocessing.Pool(processes=processes)
                validated_batches = pool.imap(_validate_batch, payloads, chunksize=1)
            try:
                for rows, validations in validated_batches:
                    for row, validation in zip(rows, validations):
                        source_rows += 1
                        if not validation.accepted:
                            reason_counts.update(validation.reasons)
                            continue
                        cleaned = _cleaned_row(row)
                        digest = _row_digest(cleaned)
                        if digest in seen_rows:
                            duplicate_rows += 1
                            reason_counts["duplicate_cleaned_row"] += 1
                            continue
                        seen_rows.add(digest)
                        writer.writerow(cleaned)
                        accepted_rows += 1
            finally:
                if pool is not None:
                    pool.close()
                    pool.join()

    temporary_output.replace(output_path)
    report: dict[str, object] = {
        "source_path": source_path.as_posix(),
        "output_path": output_path.as_posix(),
        "source_rows": source_rows,
        "accepted_rows": accepted_rows,
        "removed_rows": source_rows - accepted_rows,
        "duplicate_rows_removed": duplicate_rows,
        "minimum_mapping_confidence": minimum_confidence,
        "required_active_condition_columns": list(ACTIVE_CONDITION_COLUMNS),
        "dropped_columns": list(DROPPED_COLUMNS),
        "output_columns": list(OUTPUT_COLUMNS),
        "rejection_reason_counts": dict(sorted(reason_counts.items())),
    }
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    source_default = Path(__file__).with_name(
        "USPTO_condition_pred_category_maped_rxn.csv"
    )
    parser.add_argument("--source", type=Path, default=source_default)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("USPTO_condition_reactions_cleaned.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path(__file__).with_name("USPTO_condition_reactions_cleanup.json"),
    )
    parser.add_argument(
        "--minimum-confidence",
        type=float,
        default=DEFAULT_MINIMUM_CONFIDENCE,
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    return parser


def main() -> None:
    """Run the command-line cleanup."""
    args = _parser().parse_args()
    report = clean_dataset(
        args.source,
        args.output,
        args.report,
        minimum_confidence=args.minimum_confidence,
        processes=max(1, args.processes),
        batch_size=max(1, args.batch_size),
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
