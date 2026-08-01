"""Canonical serialization helpers for reaction observations and signatures."""

from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any, Dict


def ring_change_summary(change: Dict[str, Any]) -> str:
    """Render compact graph facts for one serialized ring change."""
    elements = tuple(str(value) for value in change.get("element_sequence") or ())
    counts: Dict[str, int] = {}
    for element in elements:
        counts[element] = counts.get(element, 0) + 1
    formula = "".join(
        element + (str(count) if count > 1 else "")
        for element, count in sorted(
            counts.items(), key=lambda value: (value[0] != "C", value[0])
        )
    )
    aromatic = "aromatic " if change.get("aromatic_after") else ""
    formed = ",".join(str(value) for value in change.get("formed_bond_types") or ())
    return (
        f"{aromatic}{len(elements)}-membered {formula} ring"
        + (f"; formed={formed}" if formed else "")
    )


def flattened_ring_change_fields(
    signature: Dict[str, Any] | None,
) -> Dict[str, Any]:
    """Expose lower-level ring observations in a compact review form."""
    topology = (signature or {}).get("topology") or {}
    changes = tuple(
        value
        for value in topology.get("ring_changes") or ()
        if isinstance(value, dict)
    )
    return {
        "ring_change_count": len(changes),
        "ring_change_summaries": " | ".join(
            ring_change_summary(change) for change in changes
        ),
    }


def signature_record_fields(analysis: Any) -> Dict[str, Any]:
    """Return RecommendationRecord keyword fields from a reaction analysis."""
    signature = analysis.reaction_signature
    evidence_fields = {
        "reaction_core": (
            asdict(analysis.reaction_core)
            if analysis.reaction_core is not None
            else None
        ),
        "reaction_evidence_candidates": tuple(
            asdict(candidate)
            for candidate in analysis.evidence_candidates
        ),
        "reaction_edit_hypotheses": tuple(
            asdict(hypothesis)
            for hypothesis in analysis.edit_hypotheses
        ),
    }
    if signature is None:
        partial = analysis.partial_product_transformation
        return {
            "reaction_signature": None,
            "fallback_descriptor": (
                asdict(analysis.fallback_descriptor)
                if analysis.fallback_descriptor is not None
                else None
            ),
            "transformation_class": analysis.transformation_class,
            "transformation_confidence": (
                partial.confidence if partial is not None else 0.0
            ),
            "family_confidence": 0.0,
            "taxonomy_definition_versions": {},
            **evidence_fields,
        }
    return {
        "reaction_signature": asdict(signature),
        "fallback_descriptor": (
            asdict(analysis.fallback_descriptor)
            if analysis.fallback_descriptor is not None
            else None
        ),
        "transformation_class": analysis.transformation_class,
        "transformation_confidence": signature.transformation_confidence,
        "family_confidence": 1.0 if analysis.named_family else 0.0,
        "taxonomy_definition_versions": dict(signature.definition_versions),
        **evidence_fields,
    }


def flattened_signature_fields(signature: Dict[str, Any] | None) -> Dict[str, Any]:
    """Expose retrieval keys and canonical nested JSON in a CSV review view."""
    value = signature or {}
    return {
        "reaction_signature_id": value.get("signature_id", ""),
        "signature_l0_exact": value.get("exact_signature_key", ""),
        "signature_l1_handle": value.get("handle_signature_key", ""),
        "signature_l2_transformation": value.get("transformation_signature_key", ""),
        "signature_l3_bond_edit": value.get("bond_edit_signature_key", ""),
        "signature_l4_environment": value.get("environment_signature_key", ""),
        "reaction_signature_json": json.dumps(
            signature, ensure_ascii=False, sort_keys=True
        )
        if signature
        else "",
    }


def flattened_fallback_fields(
    descriptor: Dict[str, Any] | None,
) -> Dict[str, Any]:
    """Expose the fallback descriptor identity and canonical review JSON."""
    value = descriptor or {}
    return {
        "fallback_descriptor_id": value.get("descriptor_id", ""),
        "fallback_descriptor_evidence_mode": value.get("evidence_mode", ""),
        "fallback_descriptor_retrieval_eligible": value.get(
            "retrieval_eligible", False
        ),
        "fallback_descriptor_json": json.dumps(
            descriptor, ensure_ascii=False, sort_keys=True
        )
        if descriptor
        else "",
    }


def flattened_reaction_core_fields(
    core: Dict[str, Any] | None,
) -> Dict[str, Any]:
    """Expose minimized-core identity and canonical review JSON."""
    value = core or {}
    quality = value.get("quality") or {}
    presentation = value.get("presentation") or {}
    summary_sections = (
        presentation.get("equation"),
        *tuple(presentation.get("bond_changes") or ()),
        *tuple(presentation.get("atom_state_changes") or ()),
        *tuple(presentation.get("retained_context") or ()),
        *tuple(presentation.get("departing_context") or ()),
        *tuple(presentation.get("appearing_context") or ()),
    )
    return {
        "reaction_core_id": value.get("core_id", ""),
        "reaction_core_exact": value.get("exact_core_key", ""),
        "reaction_core_typed": value.get("typed_core_key", ""),
        "reaction_core_shape": value.get("shape_core_key", ""),
        "reaction_core_center_transition": value.get(
            "center_transition_key",
            "",
        ),
        "reaction_core_mapping_equivalence": value.get(
            "mapping_equivalence_key", ""
        ),
        "reaction_core_quality_status": quality.get("status", ""),
        "reaction_core_quality_reasons": "; ".join(
            tuple(quality.get("review_reasons") or ())
            + tuple(quality.get("blocking_reasons") or ())
        ),
        "reaction_core_state_changes": json.dumps(
            value.get("state_changes") or (),
            ensure_ascii=False,
            sort_keys=True,
        ),
        "reaction_core_chemist_summary": " | ".join(
            str(section) for section in summary_sections if section
        ),
        "reaction_core_label": value.get("generic_label", ""),
        "reaction_core_evidence_status": value.get("evidence_status", ""),
        "reaction_core_json": json.dumps(
            core,
            ensure_ascii=False,
            sort_keys=True,
        )
        if core
        else "",
    }


__all__ = [
    "flattened_fallback_fields",
    "flattened_reaction_core_fields",
    "flattened_ring_change_fields",
    "flattened_signature_fields",
    "signature_record_fields",
    "ring_change_summary",
]
