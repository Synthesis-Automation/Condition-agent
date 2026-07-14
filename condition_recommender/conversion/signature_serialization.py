"""Temporary shared serialization helpers for Phase A reaction signatures."""

from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any, Dict


def signature_record_fields(analysis: Any) -> Dict[str, Any]:
    """Return RecommendationRecord keyword fields from a reaction analysis."""
    signature = analysis.reaction_signature
    if signature is None:
        return {
            "reaction_signature": None,
            "transformation_class": analysis.transformation_class,
            "transformation_confidence": 0.0,
            "family_confidence": 0.0,
            "taxonomy_definition_versions": {},
        }
    return {
        "reaction_signature": asdict(signature),
        "transformation_class": signature.transformation_class,
        "transformation_confidence": signature.transformation_confidence,
        "family_confidence": signature.family_confidence,
        "taxonomy_definition_versions": dict(signature.definition_versions),
    }


def flattened_signature_fields(signature: Dict[str, Any] | None) -> Dict[str, Any]:
    """Expose retrieval keys and canonical nested JSON in a CSV review view."""
    value = signature or {}
    return {
        "reaction_signature_id": value.get("signature_id", ""),
        "signature_l0_exact": value.get("exact_signature_key", ""),
        "signature_l1_handle": value.get("handle_signature_key", ""),
        "signature_l2_transformation": value.get(
            "transformation_signature_key", ""
        ),
        "signature_l3_bond_edit": value.get("bond_edit_signature_key", ""),
        "signature_l4_environment": value.get("environment_signature_key", ""),
        "reaction_signature_json": json.dumps(
            signature, ensure_ascii=False, sort_keys=True
        )
        if signature
        else "",
    }


__all__ = ["flattened_signature_fields", "signature_record_fields"]
