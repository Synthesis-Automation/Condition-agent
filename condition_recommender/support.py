"""Independent evidence-unit accounting for generic precedents."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .generic_indexing import GenericIndexedReaction


@dataclass(frozen=True)
class EvidenceSupport:
    """Support counts that distinguish correlated rows from replication."""

    observation_count: int
    reaction_count: int
    condition_series_count: int
    reference_count: int
    dataset_count: int
    independent_count: int


def evidence_unit(row: GenericIndexedReaction) -> str:
    """Use a publication when known, otherwise a canonical reaction."""
    if row.reference_id:
        return f"reference:{row.reference_id}"
    return "reaction:" + (
        row.canonical_reaction_id or row.observation_id or row.reaction_id
    )


def summarize_evidence_support(
    rows: Iterable[GenericIndexedReaction],
) -> EvidenceSupport:
    """Count raw and independent support without treating scope rows as papers."""
    values = tuple(rows)
    reactions = {
        row.canonical_reaction_id or row.observation_id or row.reaction_id
        for row in values
    }
    series = {
        row.reference_condition_series_id
        for row in values
        if row.reference_condition_series_id
    }
    references = {row.reference_id for row in values if row.reference_id}
    datasets = {row.source_dataset for row in values if row.source_dataset}
    independent = {evidence_unit(row) for row in values}
    return EvidenceSupport(
        observation_count=len(values),
        reaction_count=len(reactions),
        condition_series_count=len(series),
        reference_count=len(references),
        dataset_count=len(datasets),
        independent_count=len(independent),
    )


__all__ = ["EvidenceSupport", "evidence_unit", "summarize_evidence_support"]
