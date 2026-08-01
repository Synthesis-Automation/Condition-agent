"""One public renderer for reaction observations and interpretations."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .reaction_contextual_labels import build_contextual_transformation_label
from .reaction_display_labels import (
    _build_reaction_label,
    load_reaction_label_rendering,
)
from .reaction_models import (
    ReactionInterpretation,
    ReactionObservation,
    RenderedReactionLabel,
)
from .reaction_signatures import build_observation_signature

_REACTION_LABEL_DEFINITION_FILES = (
    "reaction_label_rendering.v1.json",
    "reaction_label_patterns.v1.json",
)


def reaction_label_definition_versions() -> dict[str, str]:
    """Return content-addressed versions of display-label definitions.

    Display labels are excluded from reaction identity, but converted records
    serialize them. Sharded conversion therefore uses this separate provenance
    contract to invalidate cached records when presentation rules change.
    """
    root = Path(__file__).with_name("definitions")
    versions = {}
    for filename in _REACTION_LABEL_DEFINITION_FILES:
        raw = (root / filename).read_bytes()
        payload = json.loads(raw)
        version = str(
            payload.get("label_schema_version")
            or payload.get("schema_version")
            or "unknown"
        )
        versions[filename] = (
            f"{version}@sha256:{hashlib.sha256(raw).hexdigest()[:16]}"
        )
    return dict(sorted(versions.items()))


def render_reaction(
    observation: ReactionObservation,
    interpretation: ReactionInterpretation | None = None,
    *,
    style: str = "unicode",
    fallback_label: str | None = None,
    fallback_status: str = "unavailable",
    fallback_detailed_label: str | None = None,
    fallback_evidence: str | None = None,
    fallback_confidence: float | None = None,
) -> RenderedReactionLabel:
    """Render exactly one evidence-ranked label for one observation."""
    contextual = (
        build_contextual_transformation_label(
            observation.reactants,
            observation.edits,
            style=style,
        )
        if observation.edits
        and observation.evidence_quality
        not in {
            "conflicting_edit_evidence",
            "conflicting_stereochemical_evidence",
        }
        else None
    )
    signature = (
        build_observation_signature(
            observation,
            contextual_product_label=contextual.after if contextual else None,
        )
        if observation.topology is not None
        and observation.completeness is not None
        and observation.completeness.status != "incomplete"
        else None
    )
    selected = interpretation.selected_candidate if interpretation else None
    interpretation_conflict = bool(
        interpretation
        and "INTERPRETATION_OBSERVATION_CONFLICT" in interpretation.warnings
    )
    selected_exact = bool(
        selected
        and selected.verification == "exact_product_reconstruction"
        and not interpretation_conflict
    )
    rendered = _build_reaction_label(
        reactants=observation.reactants,
        edits=observation.edits,
        selected_label=selected.grammar_label if selected_exact else None,
        selected_exact=selected_exact,
        grammar_id=selected.grammar_id if selected_exact else None,
        contextual_label=contextual,
        named_family=(
            interpretation.named_family if interpretation and selected_exact else None
        ),
        fallback_label=fallback_label,
        fallback_status=fallback_status,
        evidence=fallback_evidence or observation.evidence_quality,
        confidence=(
            fallback_confidence
            if fallback_confidence is not None
            else observation.evidence_confidence
        ),
        events=signature.events if signature is not None else (),
        topology=observation.topology,
        reaction_core=observation.core,
        warnings=(
            tuple(observation.warnings)
            + (tuple(interpretation.warnings) if interpretation else ())
        ),
        style=style,
        fallback_detailed_label=fallback_detailed_label,
    )
    if rendered is not None:
        return rendered
    definition = load_reaction_label_rendering()
    return RenderedReactionLabel(
        concise="Unavailable",
        detailed="No normalized edits or supported structural interpretation",
        status="unavailable",
        source="unavailable",
        clauses=(),
        evidence=observation.evidence_quality,
        confidence=observation.evidence_confidence,
        warnings=tuple(sorted(set(observation.warnings))),
        style=style,
        definition_version=str(definition["label_schema_version"]),
    )


__all__ = ["reaction_label_definition_versions", "render_reaction"]
