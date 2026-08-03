"""Stable CSV review view for generic converted records."""

from __future__ import annotations

import json
from typing import Any, Dict, Iterable, Mapping

from ..models import RecommendationRecord
from .signature_serialization import (
    flattened_fallback_fields,
    flattened_reaction_core_fields,
    flattened_ring_change_fields,
    flattened_signature_fields,
)
from .pattern_serialization import (
    REACTION_PATTERN_REVIEW_FIELDS,
    flattened_reaction_pattern_fields,
)

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
    "condition_stage_status",
    "outcome_status",
    "index_eligibility",
    "reaction_smiles",
    "reaction_label",
    "canonical_reaction_smiles",
    *REACTION_PATTERN_REVIEW_FIELDS,
    "evidence_quality",
    "reaction_completeness_status",
    "product_heavy_atom_coverage",
    "reaction_completeness_warnings",
    "reaction_completeness_json",
    "reaction_evidence_providers",
    "external_mapping_status",
    "external_mapping_provider",
    "external_mapping_confidence",
    "external_mapping_matched_hypothesis_ids",
    "external_atom_mapping_json",
    "reaction_edit_hypothesis_count",
    "reaction_edit_hypothesis_ids",
    "reaction_edit_hypotheses_json",
    "transformation_class",
    "transformation_confidence",
    "named_family",
    "family_confidence",
    "reaction_label_status",
    "reaction_label_basis",
    "reaction_label_confidence",
    "reaction_label_warnings",
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
    "ring_change_count",
    "ring_change_summaries",
    "reaction_event_count",
    "reaction_event_scope",
    "reaction_signature_id",
    "signature_l0_exact",
    "signature_l1_handle",
    "signature_l2_transformation",
    "signature_l3_bond_edit",
    "signature_l4_environment",
    "reaction_signature_json",
    "reaction_core_id",
    "reaction_core_exact",
    "reaction_core_typed",
    "reaction_core_shape",
    "reaction_core_center_transition",
    "reaction_core_mapping_equivalence",
    "reaction_core_quality_status",
    "reaction_core_quality_reasons",
    "reaction_core_state_changes",
    "reaction_core_substituent_profiles",
    "reaction_core_evidence_status",
    "reaction_core_json",
    "fallback_descriptor_id",
    "fallback_descriptor_evidence_mode",
    "fallback_descriptor_retrieval_eligible",
    "partial_transformation_key",
    "fragment_source_support_status",
    "fragment_source_support_components",
    "fragment_source_support_json",
    "fallback_descriptor_json",
    "reference",
)


def review_reaction_label_text(label: Mapping[str, Any] | None) -> str:
    """Return display text only when a reaction label is available."""
    if not label or str(label.get("status") or "").casefold() == "unavailable":
        return ""
    return str(label.get("text") or "")


def _joined(values: Iterable[Any]) -> str:
    return "|".join(str(value) for value in values)


def flatten_generic_record(record: RecommendationRecord) -> Dict[str, Any]:
    """Flatten high-value review fields while retaining nested JSONL as canonical."""
    status_counts = record.condition_resolution.get("status_counts") or {}
    topology = (record.reaction_signature or {}).get("topology") or {}
    signature = record.reaction_signature or {}
    completeness = record.reaction_completeness or {}
    external_mapping = record.external_atom_mapping or {}
    source_support = tuple(record.fragment_source_support)
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
        "condition_stage_status": record.condition_stage_status.value,
        "outcome_status": record.outcome_status.value,
        "index_eligibility": record.index_eligibility.value,
        "reaction_smiles": record.reaction_smiles,
        "reaction_label": review_reaction_label_text(record.reaction_label),
        "canonical_reaction_smiles": record.canonical_reaction_smiles or "",
        **flattened_reaction_pattern_fields(
            record.reaction_interpretation,
            named_family=record.named_family,
        ),
        "evidence_quality": record.evidence_quality,
        "reaction_completeness_status": completeness.get("status", ""),
        "product_heavy_atom_coverage": (
            completeness.get("product_heavy_atom_coverage")
            if completeness.get("product_heavy_atom_coverage") is not None
            else ""
        ),
        "reaction_completeness_warnings": _joined(completeness.get("warnings") or ()),
        "reaction_completeness_json": (
            json.dumps(
                completeness,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            if completeness
            else ""
        ),
        "reaction_evidence_providers": _joined(
            candidate.get("provider", "")
            for candidate in record.reaction_evidence_candidates
            if candidate.get("provider")
        ),
        "external_mapping_status": external_mapping.get("status", ""),
        "external_mapping_provider": (
            (external_mapping.get("provider") or {}).get("provider_id", "")
        ),
        "external_mapping_confidence": (
            external_mapping.get("mapper_confidence")
            if external_mapping.get("mapper_confidence") is not None
            else ""
        ),
        "external_mapping_matched_hypothesis_ids": _joined(
            external_mapping.get("matched_hypothesis_ids") or ()
        ),
        "external_atom_mapping_json": (
            json.dumps(
                external_mapping,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            if external_mapping
            else ""
        ),
        "reaction_edit_hypothesis_count": len(
            record.reaction_edit_hypotheses
        ),
        "reaction_edit_hypothesis_ids": _joined(
            hypothesis.get("hypothesis_id", "")
            for hypothesis in record.reaction_edit_hypotheses
            if hypothesis.get("hypothesis_id")
        ),
        "reaction_edit_hypotheses_json": (
            json.dumps(
                record.reaction_edit_hypotheses,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            if record.reaction_edit_hypotheses
            else ""
        ),
        "transformation_class": record.transformation_class or "",
        "transformation_confidence": record.transformation_confidence,
        "named_family": record.named_family or "",
        "family_confidence": record.family_confidence,
        "reaction_label_status": (
            record.reaction_label.get("status", "")
            if record.reaction_label
            else ""
        ),
        "reaction_label_basis": (
            record.reaction_label.get("basis", "")
            if record.reaction_label
            else ""
        ),
        "reaction_label_confidence": (
            record.reaction_label.get("confidence", "")
            if record.reaction_label
            else ""
        ),
        "reaction_label_warnings": _joined(
            record.reaction_label.get("warnings", ())
            if record.reaction_label
            else ()
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
        **flattened_ring_change_fields(record.reaction_signature),
        "reaction_event_count": signature.get("event_count", ""),
        "reaction_event_scope": signature.get("event_scope", ""),
        **flattened_signature_fields(record.reaction_signature),
        **flattened_reaction_core_fields(record.reaction_core),
        "partial_transformation_key": (
            (record.fallback_descriptor or {}).get(
                "partial_transformation_key",
                "",
            )
        ),
        "fragment_source_support_status": _joined(
            value.get("status", "") for value in source_support
        ),
        "fragment_source_support_components": _joined(
            sorted(
                {
                    str(identifier)
                    for value in source_support
                    for identifier in (
                        value.get("component_raw_identifiers") or ()
                    )
                }
            )
        ),
        "fragment_source_support_json": (
            json.dumps(
                source_support,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            if source_support
            else ""
        ),
        **flattened_fallback_fields(record.fallback_descriptor),
        "reference": record.source.get("reference", ""),
    }


__all__ = ["GENERIC_REVIEW_FIELDS", "flatten_generic_record"]
