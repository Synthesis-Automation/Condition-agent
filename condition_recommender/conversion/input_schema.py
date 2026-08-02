"""Common streaming input schema for heterogeneous reaction datasets."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

from condition_registry import ConditionComponentInput, ConditionProcessStage

from ..condition_normalization import optional_float, split_identifiers

_COLUMN_ALIASES: Dict[str, Tuple[str, ...]] = {
    "source_dataset": ("source_dataset",),
    "source_path": ("source_path",),
    "source_row_number": ("source_row_number",),
    "reaction_id": ("reaction_id", "reactionid", "id"),
    "reaction_type": ("reaction_type", "reactiontype", "family"),
    "reaction_smiles": ("reaction_smiles", "reactionsmiles", "rxn_smiles"),
    "yield_pct": ("yield_pct", "yield", "yield_percent"),
    "temperature_c": ("temperature_c", "temperature", "temp_c"),
    "time_h": ("time_h", "time", "duration_h"),
    "reference": ("reference", "citation", "source_reference"),
    "reactant_cas": ("reactant_cas",),
    "product_cas": ("product_cas",),
    "reagent_cas": ("reagent_cas", "reagents_cas"),
    "catalyst_cas": ("catalyst_cas", "catalysts_cas"),
    "solvent_cas": ("solvent_cas", "solvents_cas"),
    "experimental_procedure": ("experimental_procedure", "procedure"),
    "stages": ("stages",),
    "steps": ("steps",),
    "notes": ("notes",),
}


def _normalized_headers(row: Mapping[str, Any]) -> Dict[str, str]:
    return {str(key).strip().lower(): str(key) for key in row if key is not None}


def _value(row: Mapping[str, Any], field_name: str, headers: Mapping[str, str]) -> str:
    for alias in _COLUMN_ALIASES[field_name]:
        source_key = headers.get(alias)
        if source_key is not None:
            return str(row.get(source_key) or "").strip()
    return ""


@dataclass(frozen=True)
class RawReactionRecord:
    """Source-faithful row normalized only at the dataset contract boundary."""

    source_dataset: str
    source_path: str
    source_row_number: int
    reaction_id: str
    source_declared_family: str
    reaction_smiles: str
    yield_pct: Optional[float]
    temperature_c: Optional[float]
    time_h: Optional[float]
    reference: str
    reactant_cas: Tuple[str, ...]
    product_cas: Tuple[str, ...]
    reagent_cas: Tuple[str, ...]
    catalyst_cas: Tuple[str, ...]
    solvent_cas: Tuple[str, ...]
    experimental_procedure: str
    stages: str
    steps: str
    notes: str
    condition_component_inputs: Tuple[ConditionComponentInput, ...] = ()
    condition_process_stages: Tuple[ConditionProcessStage, ...] = ()
    condition_declared_absences: Tuple[str, ...] = ()
    primary_outcome_type: str = ""
    upstream_observation_id: str = ""
    raw_fields: Dict[str, Any] = field(default_factory=dict)
    schema_version: str = "1.0"


def adapt_row(
    row: Mapping[str, Any],
    *,
    source_dataset: str,
    source_path: str,
    source_row_number: int,
) -> RawReactionRecord:
    """Adapt one source row without performing chemistry classification."""
    headers = _normalized_headers(row)
    preserved_row_number = _value(row, "source_row_number", headers)
    try:
        effective_row_number = (
            int(preserved_row_number) if preserved_row_number else source_row_number
        )
    except ValueError:
        effective_row_number = source_row_number
    effective_dataset = _value(row, "source_dataset", headers) or source_dataset
    effective_source_path = _value(row, "source_path", headers) or source_path
    reaction_id = (
        _value(row, "reaction_id", headers)
        or f"{effective_dataset}:row-{effective_row_number}"
    )
    return RawReactionRecord(
        source_dataset=effective_dataset,
        source_path=effective_source_path,
        source_row_number=effective_row_number,
        reaction_id=reaction_id,
        source_declared_family=_value(row, "reaction_type", headers),
        reaction_smiles=_value(row, "reaction_smiles", headers),
        yield_pct=optional_float(_value(row, "yield_pct", headers)),
        temperature_c=optional_float(_value(row, "temperature_c", headers)),
        time_h=optional_float(_value(row, "time_h", headers)),
        reference=_value(row, "reference", headers),
        reactant_cas=split_identifiers(_value(row, "reactant_cas", headers)),
        product_cas=split_identifiers(_value(row, "product_cas", headers)),
        reagent_cas=split_identifiers(_value(row, "reagent_cas", headers)),
        catalyst_cas=split_identifiers(_value(row, "catalyst_cas", headers)),
        solvent_cas=split_identifiers(_value(row, "solvent_cas", headers)),
        experimental_procedure=_value(row, "experimental_procedure", headers),
        stages=_value(row, "stages", headers),
        steps=_value(row, "steps", headers),
        notes=_value(row, "notes", headers),
        raw_fields={str(key): value for key, value in row.items() if key is not None},
    )


def iter_csv_records(path: str | Path) -> Iterator[RawReactionRecord]:
    """Stream adapted records from one CSV file."""
    csv_path = Path(path)
    dataset = csv_path.stem
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            yield adapt_row(
                row,
                source_dataset=dataset,
                source_path=str(csv_path),
                source_row_number=row_number,
            )


def discover_csv_datasets(path: str | Path) -> Tuple[Path, ...]:
    """Return deterministic CSV paths for a file or recursive directory input."""
    source = Path(path)
    if source.is_file():
        return (source,) if source.suffix.lower() == ".csv" else ()
    if not source.is_dir():
        return ()
    return tuple(
        sorted(
            (
                item
                for item in source.rglob("*")
                if item.is_file() and item.suffix.lower() == ".csv"
            ),
            key=lambda item: item.relative_to(source).as_posix().casefold(),
        )
    )


def iter_conversion_records(path: str | Path) -> Iterator[RawReactionRecord]:
    """Stream either a raw CSV or a chemistry-free intermediate artifact."""
    source = Path(path)
    if source.name.casefold().endswith(".observations.jsonl.gz"):
        from .intermediate import iter_intermediate_records

        yield from iter_intermediate_records(source)
        return
    yield from iter_csv_records(source)


def discover_conversion_datasets(path: str | Path) -> Tuple[Path, ...]:
    """Discover raw CSVs or preprocessed observation files deterministically."""
    source = Path(path)

    def supported(item: Path) -> bool:
        return item.suffix.casefold() == ".csv" or item.name.casefold().endswith(
            ".observations.jsonl.gz"
        )

    if source.is_file():
        return (source,) if supported(source) else ()
    if not source.is_dir():
        return ()
    return tuple(
        sorted(
            (item for item in source.rglob("*") if item.is_file() and supported(item)),
            key=lambda item: item.relative_to(source).as_posix().casefold(),
        )
    )


__all__ = [
    "RawReactionRecord",
    "adapt_row",
    "discover_conversion_datasets",
    "discover_csv_datasets",
    "iter_conversion_records",
    "iter_csv_records",
]
