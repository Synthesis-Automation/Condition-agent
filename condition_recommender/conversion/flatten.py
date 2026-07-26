"""Stable CSV review view for generic converted records."""

from __future__ import annotations

import json
from typing import Any, Dict, Iterable

from ..models import RecommendationRecord
from .signature_serialization import flattened_signature_fields

GENERIC_REVIEW_FIELDS = (
    "schema_version",
    "observation_id",
    "reaction_id",
    "canonical_reaction_id",
    "source_dataset",
    "source_row_number",
    "source_declared_family",
    "reference_id",
    "reference_resolution_status",
    "admission_tier",
    "admission_reasons",
    "chemistry_status",
    "condition_status",
    "outcome_status",
    "index_eligibility",
    "reaction_smiles",
    "canonical_reaction_smiles",
    "evidence_quality",
    "transformation_class",
    "transformation_confidence",
    "named_family",
    "family_confidence",
    "reaction_label",
    "reaction_label_status",
    "reaction_display_label_status",
    "reaction_display_label_detailed",
    "reaction_display_label_json",
    "yield_pct",
    "temperature_c",
    "time_h",
    "condition_recipe_id",
    "condition_recipe_core_id",
    "reference_condition_series_id",
    "legacy_condition_recipe_id",
    "raw_recipe_id",
    "catalyst_cas",
    "reagent_cas",
    "solvent_cas",
    "condition_resolved_count",
    "condition_unresolved_count",
    "condition_invalid_count",
    "condition_identity_uncertainty",
    "reaction_scope",
    "formed_ring_sizes",
    "ring_count_delta",
    "reaction_event_count",
    "reaction_event_scope",
    "reaction_signature_id",
    "signature_l0_exact",
    "signature_l1_handle",
    "signature_l2_transformation",
    "signature_l3_bond_edit",
    "signature_l4_environment",
    "reaction_signature_json",
    "reference",
)


def _joined(values: Iterable[Any]) -> str:
    return "|".join(str(value) for value in values)


def flatten_generic_record(record: RecommendationRecord) -> Dict[str, Any]:
    """Flatten high-value review fields while retaining nested JSONL as canonical."""
    status_counts = record.condition_resolution.get("status_counts") or {}
    topology = (record.reaction_signature or {}).get("topology") or {}
    signature = record.reaction_signature or {}
    return {
        "schema_version": record.schema_version,
        "observation_id": record.observation_id,
        "reaction_id": record.reaction_id,
        "canonical_reaction_id": record.canonical_reaction_id or "",
        "source_dataset": record.source_dataset,
        "source_row_number": record.source_row_number,
        "source_declared_family": record.source_declared_family,
        "reference_id": record.reference_id,
        "reference_resolution_status": (
            (record.reference_identity or {}).get("resolution_status", "")
        ),
        "admission_tier": record.admission_tier.value,
        "admission_reasons": _joined(record.admission_reasons),
        "chemistry_status": record.chemistry_status.value,
        "condition_status": record.condition_status.value,
        "outcome_status": record.outcome_status.value,
        "index_eligibility": record.index_eligibility.value,
        "reaction_smiles": record.reaction_smiles,
        "canonical_reaction_smiles": record.canonical_reaction_smiles or "",
        "evidence_quality": record.evidence_quality,
        "transformation_class": record.transformation_class or "",
        "transformation_confidence": record.transformation_confidence,
        "named_family": record.named_family or "",
        "family_confidence": record.family_confidence,
        "reaction_label": record.reaction_label or "",
        "reaction_label_status": record.reaction_label_status,
        "reaction_display_label_status": (
            record.reaction_display_label.get("status", "")
            if record.reaction_display_label
            else ""
        ),
        "reaction_display_label_detailed": (
            record.reaction_display_label.get("detailed", "")
            if record.reaction_display_label
            else ""
        ),
        "reaction_display_label_json": (
            json.dumps(
                record.reaction_display_label,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            if record.reaction_display_label
            else ""
        ),
        "yield_pct": record.yield_pct if record.yield_pct is not None else "",
        "temperature_c": (
            record.temperature_c if record.temperature_c is not None else ""
        ),
        "time_h": record.time_h if record.time_h is not None else "",
        "condition_recipe_id": record.resolved_recipe_id or record.conditions.recipe_id,
        "condition_recipe_core_id": record.resolved_recipe_core_id,
        "reference_condition_series_id": record.reference_condition_series_id,
        "legacy_condition_recipe_id": record.conditions.recipe_id,
        "raw_recipe_id": record.raw_recipe_id,
        "catalyst_cas": _joined(record.conditions.catalyst_cas),
        "reagent_cas": _joined(record.conditions.reagent_cas),
        "solvent_cas": _joined(record.conditions.solvent_cas),
        "condition_resolved_count": status_counts.get("resolved", 0),
        "condition_unresolved_count": status_counts.get("unresolved", 0),
        "condition_invalid_count": status_counts.get("invalid_identifier", 0),
        "condition_identity_uncertainty": record.condition_resolution.get(
            "has_uncertainty", False
        ),
        "reaction_scope": topology.get("reaction_scope", ""),
        "formed_ring_sizes": _joined(topology.get("formed_ring_sizes") or ()),
        "ring_count_delta": (
            topology["ring_count_delta"]
            if topology.get("ring_count_delta") is not None
            else ""
        ),
        "reaction_event_count": signature.get("event_count", ""),
        "reaction_event_scope": signature.get("event_scope", ""),
        **flattened_signature_fields(record.reaction_signature),
        "reference": record.source.get("reference", ""),
    }


__all__ = ["GENERIC_REVIEW_FIELDS", "flatten_generic_record"]
