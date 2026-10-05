"""Explicit registry and deterministic detection for source adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

from .adapters import (
    HiTeaCsvAdapter,
    LiteratureCsvAdapter,
    UsptoConditionCsvAdapter,
    WeakLabelCsvAdapter,
)
from .adapters.base import SourceAdapter, read_headers
from .adapters.reaction_smiles import (
    HigherLevelAbstractionCsvAdapter,
    ReactionSmilesCsvAdapter,
)

_ADAPTERS: Dict[str, SourceAdapter] = {
    adapter.adapter_id: adapter
    for adapter in (
        LiteratureCsvAdapter(),
        HiTeaCsvAdapter(),
        UsptoConditionCsvAdapter(),
        WeakLabelCsvAdapter(),
        ReactionSmilesCsvAdapter(),
        HigherLevelAbstractionCsvAdapter(),
    )
}


def adapter_ids() -> Tuple[str, ...]:
    """Return stable adapter IDs available to the preprocessor."""
    return tuple(sorted(_ADAPTERS))


def get_adapter(adapter_id: str) -> SourceAdapter:
    """Return one explicitly registered adapter."""
    try:
        return _ADAPTERS[adapter_id]
    except KeyError as exc:
        raise ValueError(f"Unknown source adapter: {adapter_id}") from exc


def detect_adapter(path: str | Path) -> SourceAdapter:
    """Select the unique adapter whose required columns are all present."""
    source = Path(path)
    headers = set(read_headers(source))
    if {"id", "rxn_smiles"}.issubset(headers):
        # These exports share a header but have explicitly different source
        # semantics. Unknown exports require an explicit adapter selection.
        if source.name == "uspto_higher-level.csv":
            return get_adapter("higher_level_abstraction_csv.v1")
        if source.name == "uspto_original.csv":
            return get_adapter("reaction_smiles_csv.v1")
        raise ValueError("Select an explicit reaction CSV adapter for this export")
    matches = tuple(
        adapter
        for adapter in _ADAPTERS.values()
        if set(adapter.required_columns).issubset(headers)
    )
    if not matches:
        raise ValueError(
            f"No registered source adapter matches the columns in {source.name}"
        )
    if len(matches) > 1:
        raise ValueError(
            "Ambiguous source format; matching adapters: "
            + ", ".join(sorted(adapter.adapter_id for adapter in matches))
        )
    return matches[0]


__all__ = ["adapter_ids", "detect_adapter", "get_adapter"]
