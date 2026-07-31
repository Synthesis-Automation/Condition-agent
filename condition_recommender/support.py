"""Independent evidence-unit accounting for generic precedents."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
from typing import Any, Iterable

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
    mapping_equivalence_count: int
    mapping_equivalence_row_count: int
    mapping_deduplicated_independent_count: int


_RULES_PATH = (
    Path(__file__).with_name("definitions") / "evidence_support.v1.json"
)


@lru_cache(maxsize=1)
def load_evidence_support_rules() -> dict[str, Any]:
    """Load the deterministic independent-support identity policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported evidence-support schema")
    if str(rules.get("definition_id") or "") != "evidence_support.v1":
        raise ValueError("unexpected evidence-support definition ID")
    if tuple(rules.get("independence_priority") or ()) != (
        "reference_id",
        "reaction_core.mapping_equivalence_key",
        "canonical_reaction_id",
        "observation_id",
        "reaction_id",
    ):
        raise ValueError("evidence-support independence priority is invalid")
    if rules.get("preserve_all_observations") is not True:
        raise ValueError("evidence-support policy must preserve observations")
    if rules.get("mapping_equivalence_applies_without_reference") is not True:
        raise ValueError("mapping equivalence must deduplicate unreferenced rows")
    return rules


def mapping_equivalence_key(row: GenericIndexedReaction) -> str:
    """Return a validated map-number-independent core key when available."""
    value = str(row.reaction_core.get("mapping_equivalence_key") or "")
    return value if value.startswith("RME1:") else ""


def _legacy_evidence_unit(row: GenericIndexedReaction) -> str:
    if row.reference_id:
        return f"reference:{row.reference_id}"
    return "reaction:" + (
        row.canonical_reaction_id or row.observation_id or row.reaction_id
    )


def evidence_unit(row: GenericIndexedReaction) -> str:
    """Use publication identity, then mapping-equivalent reaction identity."""
    load_evidence_support_rules()
    if row.reference_id:
        return f"reference:{row.reference_id}"
    equivalence = mapping_equivalence_key(row)
    if equivalence:
        return f"mapping_equivalence:{equivalence}"
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
    legacy_independent = {_legacy_evidence_unit(row) for row in values}
    mapping_rows = tuple(
        row for row in values if mapping_equivalence_key(row)
    )
    mapping_equivalences = {
        mapping_equivalence_key(row) for row in mapping_rows
    }
    return EvidenceSupport(
        observation_count=len(values),
        reaction_count=len(reactions),
        condition_series_count=len(series),
        reference_count=len(references),
        dataset_count=len(datasets),
        independent_count=len(independent),
        mapping_equivalence_count=len(mapping_equivalences),
        mapping_equivalence_row_count=len(mapping_rows),
        mapping_deduplicated_independent_count=max(
            0,
            len(legacy_independent) - len(independent),
        ),
    )


__all__ = [
    "EvidenceSupport",
    "evidence_unit",
    "load_evidence_support_rules",
    "mapping_equivalence_key",
    "summarize_evidence_support",
]
