"""Extract deterministic CAS-number/SMILES pairs from literature CSV files."""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .cas_number_extractor import find_cas_numbers_in_text


OUTPUT_COLUMNS = (
    "cas_no",
    "compound_smiles",
    "reaction_id",
    "citation",
    "source_role",
)
CAS_KEYS = ("casrn", "casno", "casnumber", "casregistrynumber", "casregistry", "cas")
SMILES_KEYS = ("smiles", "canonicalsmiles", "isomericsmiles", "compoundsmiles")
JSON_ROLE_COLUMNS = {
    "reactant": "reactants_json",
    "product": "products_json",
    "reagent": "reagents_json",
    "catalyst": "catalysts_json",
    "solvent": "solvents_json",
}


@dataclass(frozen=True, order=True)
class CASSmilesPair:
    """One source-observed CAS/structure pair with reaction provenance."""

    cas_no: str
    compound_smiles: str
    reaction_id: str = ""
    citation: str = ""
    source_role: str = ""

    def to_row(self) -> dict[str, str]:
        """Return the role-aware output representation."""

        return {
            "cas_no": self.cas_no,
            "compound_smiles": self.compound_smiles,
            "reaction_id": self.reaction_id,
            "citation": self.citation,
            "source_role": self.source_role,
        }


@dataclass(frozen=True)
class CSVExtractionResult:
    """Extraction result and audit counts for one CSV file."""

    source_file: Path
    rows_read: int
    pair_occurrences: int
    pairs: tuple[CASSmilesPair, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class FolderExtractionResult:
    """Aggregated extraction result for a folder of CSV files."""

    files_read: int
    rows_read: int
    pair_occurrences: int
    pairs: tuple[CASSmilesPair, ...]
    warnings: tuple[str, ...]

    @property
    def conflicting_cas_count(self) -> int:
        """Return how many CAS numbers are associated with multiple SMILES."""

        structures: dict[str, set[str]] = {}
        for pair in self.pairs:
            structures.setdefault(pair.cas_no, set()).add(pair.compound_smiles)
        return sum(len(smiles) > 1 for smiles in structures.values())


def discover_csv_files(
    folder_path: str | Path,
    *,
    excluded_paths: Sequence[str | Path] = (),
) -> list[Path]:
    """Return recursively discovered CSV files in deterministic path order."""

    folder = Path(folder_path)
    if not folder.is_dir():
        return []
    excluded = {Path(item).resolve() for item in excluded_paths}
    return [
        path
        for path in sorted(folder.rglob("*.csv"), key=lambda item: str(item).casefold())
        if path.is_file() and path.resolve() not in excluded
    ]


def extract_cas_smiles_pairs_from_csv(path: str | Path) -> CSVExtractionResult:
    """Extract all distinct, explicitly associated CAS/SMILES pairs from a CSV.

    JSON objects are inspected recursively, not only in the currently known
    ``reactants_json`` and ``products_json`` columns. Flat ``*_cas`` and
    ``*_smiles`` column pairs are supported as a fallback when the corresponding
    structured role column is absent. Ambiguous flat many-to-many values are
    skipped instead of inventing a CAS-to-structure correspondence.
    """

    source = Path(path)
    warnings: list[str] = []
    pairs: set[CASSmilesPair] = set()
    pair_occurrences = 0
    rows_read = 0

    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return CSVExtractionResult(source, 0, 0, (), ("CSV has no header.",))

        flat_columns = _matching_flat_columns(reader.fieldnames)
        for row_number, row in enumerate(reader, start=2):
            rows_read += 1
            structured_roles: set[str] = set()
            reaction_id = str(row.get("reaction_id") or "").strip()
            citation = str(row.get("citation") or "").strip()

            for column_name, raw_value in row.items():
                if column_name is None or not _looks_like_json(raw_value):
                    continue
                try:
                    value = json.loads(str(raw_value))
                except (TypeError, ValueError) as exc:
                    if _normalized_key(column_name).endswith("json"):
                        warnings.append(
                            f"Row {row_number}, {column_name}: invalid JSON ({exc})."
                        )
                    continue

                role = _role_for_json_column(column_name)
                if role is not None:
                    structured_roles.add(role)
                extracted = tuple(
                    _pairs_from_json_value(
                        value,
                        reaction_id=reaction_id,
                        citation=citation,
                        source_role=role or "",
                    )
                )
                pair_occurrences += len(extracted)
                pairs.update(extracted)

            for role, cas_column, smiles_column in flat_columns:
                if role in structured_roles:
                    continue
                extracted, warning = _pairs_from_flat_values(
                    row.get(cas_column),
                    row.get(smiles_column),
                    reaction_id=reaction_id,
                    citation=citation,
                    source_role=role,
                )
                if warning:
                    warnings.append(
                        f"Row {row_number}, {cas_column}/{smiles_column}: {warning}"
                    )
                pair_occurrences += len(extracted)
                pairs.update(extracted)

    return CSVExtractionResult(
        source_file=source,
        rows_read=rows_read,
        pair_occurrences=pair_occurrences,
        pairs=tuple(sorted(pairs, key=_pair_sort_key)),
        warnings=tuple(warnings),
    )


def extract_cas_smiles_pairs_from_folder(
    folder_path: str | Path,
    *,
    excluded_paths: Sequence[str | Path] = (),
) -> FolderExtractionResult:
    """Extract and deduplicate pairs from every CSV below ``folder_path``."""

    files = discover_csv_files(folder_path, excluded_paths=excluded_paths)
    pairs: set[CASSmilesPair] = set()
    warnings: list[str] = []
    rows_read = 0
    pair_occurrences = 0

    for path in files:
        try:
            result = extract_cas_smiles_pairs_from_csv(path)
        except (OSError, UnicodeError, csv.Error) as exc:
            warnings.append(f"{path}: failed to read CSV ({exc}).")
            continue
        rows_read += result.rows_read
        pair_occurrences += result.pair_occurrences
        pairs.update(result.pairs)
        warnings.extend(f"{path.name}: {warning}" for warning in result.warnings)

    return FolderExtractionResult(
        files_read=len(files),
        rows_read=rows_read,
        pair_occurrences=pair_occurrences,
        pairs=tuple(sorted(pairs, key=_pair_sort_key)),
        warnings=tuple(warnings),
    )


def write_cas_smiles_pairs(
    pairs: Iterable[CASSmilesPair],
    output_path: str | Path,
) -> int:
    """Write exactly one deterministic record per CAS number.

    If a CAS number has conflicting structures, the structure supported by the
    most distinct reaction records is selected. Ties use lexical SMILES order.
    A reactant-role observation is preferred for retained provenance, followed
    by reaction ID and citation.
    """

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    selected_pairs = _select_one_record_per_cas(pairs)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS))
        writer.writeheader()
        writer.writerows(pair.to_row() for pair in selected_pairs)
    return len(selected_pairs)


def _select_one_record_per_cas(
    pairs: Iterable[CASSmilesPair],
) -> list[CASSmilesPair]:
    records_by_cas: dict[str, set[CASSmilesPair]] = {}
    for pair in pairs:
        records_by_cas.setdefault(pair.cas_no, set()).add(pair)

    selected: list[CASSmilesPair] = []
    for records in records_by_cas.values():
        smiles_counts = {
            smiles: len(
                {
                    (record.reaction_id, record.citation)
                    for record in records
                    if record.compound_smiles == smiles
                }
            )
            for smiles in {record.compound_smiles for record in records}
        }
        selected_smiles = min(
            smiles_counts,
            key=lambda smiles: (-smiles_counts[smiles], smiles),
        )
        selected.append(
            min(
                (
                    record
                    for record in records
                    if record.compound_smiles == selected_smiles
                ),
                key=lambda record: (
                    _source_role_rank(record.source_role),
                    record.reaction_id,
                    record.citation,
                ),
            )
        )
    return sorted(selected, key=_pair_sort_key)


def _pairs_from_json_value(
    value: Any,
    *,
    reaction_id: str,
    citation: str,
    source_role: str,
) -> Iterator[CASSmilesPair]:
    if isinstance(value, Mapping):
        normalized = {_normalized_key(str(key)): item for key, item in value.items()}
        cas_value = next(
            (normalized[key] for key in CAS_KEYS if key in normalized), None
        )
        smiles_value = next(
            (normalized[key] for key in SMILES_KEYS if key in normalized), None
        )
        if cas_value is not None and smiles_value is not None:
            yield from _associate_values(
                cas_value,
                smiles_value,
                reaction_id=reaction_id,
                citation=citation,
                source_role=source_role,
            )
        for child in value.values():
            yield from _pairs_from_json_value(
                child,
                reaction_id=reaction_id,
                citation=citation,
                source_role=source_role,
            )
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from _pairs_from_json_value(
                child,
                reaction_id=reaction_id,
                citation=citation,
                source_role=source_role,
            )


def _associate_values(
    cas_value: Any,
    smiles_value: Any,
    *,
    reaction_id: str,
    citation: str,
    source_role: str,
) -> Iterator[CASSmilesPair]:
    cas_numbers = _cas_numbers(cas_value)
    smiles_values = _smiles_values(smiles_value)
    if not cas_numbers or not smiles_values:
        return
    if len(cas_numbers) == 1:
        for smiles in smiles_values:
            yield CASSmilesPair(
                cas_numbers[0],
                smiles,
                reaction_id,
                citation,
                source_role,
            )
        return
    if len(smiles_values) == 1:
        for cas_no in cas_numbers:
            yield CASSmilesPair(
                cas_no,
                smiles_values[0],
                reaction_id,
                citation,
                source_role,
            )
        return
    if len(cas_numbers) == len(smiles_values):
        for cas_no, smiles in zip(cas_numbers, smiles_values):
            yield CASSmilesPair(
                cas_no,
                smiles,
                reaction_id,
                citation,
                source_role,
            )


def _pairs_from_flat_values(
    cas_value: Any,
    smiles_value: Any,
    *,
    reaction_id: str,
    citation: str,
    source_role: str,
) -> tuple[tuple[CASSmilesPair, ...], str | None]:
    cas_numbers = _cas_numbers(cas_value)
    smiles_text = str(smiles_value or "").strip()
    if not cas_numbers or not smiles_text:
        return (), None
    if len(cas_numbers) == 1:
        return (
            CASSmilesPair(
                cas_numbers[0],
                smiles_text,
                reaction_id,
                citation,
                source_role,
            ),
        ), None

    smiles_values = [item.strip() for item in smiles_text.split(".") if item.strip()]
    if len(smiles_values) != len(cas_numbers):
        return (), "ambiguous multi-value fields were skipped."
    return (
        tuple(
            CASSmilesPair(
                cas_no,
                smiles,
                reaction_id,
                citation,
                source_role,
            )
            for cas_no, smiles in zip(cas_numbers, smiles_values)
        ),
        None,
    )


def _cas_numbers(value: Any) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    scalar_values = _flatten_scalars(value)
    for scalar in scalar_values:
        for cas_no in find_cas_numbers_in_text(str(scalar)):
            if cas_no not in seen:
                seen.add(cas_no)
                ordered.append(cas_no)
    return ordered


def _smiles_values(value: Any) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for scalar in _flatten_scalars(value):
        smiles = str(scalar or "").strip()
        if smiles and smiles not in seen:
            seen.add(smiles)
            ordered.append(smiles)
    return ordered


def _flatten_scalars(value: Any) -> Iterator[Any]:
    if isinstance(value, (list, tuple, set)):
        for child in value:
            yield from _flatten_scalars(child)
    elif value is not None:
        yield value


def _matching_flat_columns(fieldnames: Sequence[str]) -> list[tuple[str, str, str]]:
    by_normalized = {_normalized_key(name): name for name in fieldnames if name}
    matched: list[tuple[str, str, str]] = []
    for normalized, cas_column in by_normalized.items():
        split = _split_cas_column_name(normalized)
        if split is None:
            continue
        role = split
        smiles_column = by_normalized.get(f"{role}smiles")
        if smiles_column is not None:
            matched.append((role.rstrip("s"), cas_column, smiles_column))
    return matched


def _split_cas_column_name(normalized: str) -> str | None:
    for suffix in ("casregistrynumber", "casnumber", "casrn", "casno", "cas"):
        if normalized.endswith(suffix) and len(normalized) > len(suffix):
            return normalized[: -len(suffix)]
    return None


def _role_for_json_column(column_name: str) -> str | None:
    normalized = _normalized_key(column_name)
    for role, known_column in JSON_ROLE_COLUMNS.items():
        if normalized == _normalized_key(known_column):
            return role
    return None


def _looks_like_json(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    stripped = value.lstrip()
    return stripped.startswith("{") or stripped.startswith("[")


def _normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.casefold())


def _source_role_rank(source_role: str) -> tuple[int, str]:
    normalized = str(source_role or "").strip().casefold()
    order = {
        "reactant": 0,
        "starting_material": 1,
        "substrate": 2,
        "reagent": 3,
        "catalyst": 4,
        "solvent": 5,
        "product": 6,
        "": 7,
    }
    return order.get(normalized, 7), normalized


def _pair_sort_key(
    pair: CASSmilesPair,
) -> tuple[tuple[int, ...], str, str, str, str]:
    return (
        tuple(int(part) for part in pair.cas_no.split("-")),
        pair.compound_smiles,
        pair.reaction_id,
        pair.citation,
        pair.source_role,
    )
