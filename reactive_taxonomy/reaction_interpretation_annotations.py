"""Optional interpretation annotations over structural reconstructions."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Tuple


REACTION_INTERPRETATION_ANNOTATION_SCHEMA_VERSION = "1.0"
REACTION_INTERPRETATION_ANNOTATION_DEFINITION_VERSION = (
    "reaction_interpretation_annotations.v1"
)
_PATH = (
    Path(__file__).with_name("definitions")
    / "reaction_interpretation_annotations.v1.json"
)


@lru_cache(maxsize=1)
def load_reaction_interpretation_annotations() -> Tuple[Dict[str, Any], ...]:
    """Load optional interpretation and rendering annotations."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if (
        payload.get("schema_version")
        != REACTION_INTERPRETATION_ANNOTATION_SCHEMA_VERSION
    ):
        raise ValueError("Unsupported reaction-interpretation annotation schema")
    annotations = payload.get("annotations")
    if not isinstance(annotations, list) or not annotations:
        raise ValueError(
            "Reaction-interpretation annotations must be a nonempty list"
        )
    identifiers = [str(item.get("id") or "") for item in annotations]
    if any(not identifier for identifier in identifiers) or len(identifiers) != len(
        set(identifiers)
    ):
        raise ValueError(
            "Invalid or duplicate reaction-interpretation annotation ID"
        )
    structural_keys = {
        "operator_id",
        "operator_slot_bindings",
        "slot_relationships",
        "slots",
    }
    for annotation in annotations:
        annotation_id = str(annotation["id"])
        if structural_keys.intersection(annotation):
            raise ValueError(
                f"Interpretation annotation contains structural execution data: "
                f"{annotation_id}"
            )
        if not str(annotation.get("reconstruction_rule_id") or ""):
            raise ValueError(
                "Interpretation annotation has no reconstruction rule: "
                f"{annotation_id}"
            )
        if not str(annotation.get("transformation_class") or ""):
            raise ValueError(
                "Interpretation annotation has no transformation class: "
                f"{annotation_id}"
            )
        if not isinstance(annotation.get("role_bindings"), dict) or not annotation[
            "role_bindings"
        ]:
            raise ValueError(
                f"Interpretation annotation has no semantic roles: {annotation_id}"
            )
    return tuple(
        sorted(
            (dict(item) for item in annotations),
            key=lambda item: (-int(item.get("priority", 0)), str(item["id"])),
        )
    )


__all__ = [
    "REACTION_INTERPRETATION_ANNOTATION_DEFINITION_VERSION",
    "REACTION_INTERPRETATION_ANNOTATION_SCHEMA_VERSION",
    "load_reaction_interpretation_annotations",
]
