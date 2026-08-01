"""Optional interpretation annotations over structural observations."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Mapping, Tuple

from .reaction_edits import normalize_predicted_edits
from .reaction_environments import build_reaction_family_environment
from .reaction_interpretation_annotations import (
    load_reaction_interpretation_annotations,
)
from .reaction_labels import render_reaction_label
from .reaction_models import (
    ReactionInterpretationCandidate,
    ReactionInterpretation,
    ReactionObservation,
    ReactionPartner,
    ReactionReconstructionCandidate,
    ReactionSiteReference,
)
from .reaction_products import build_product_connection


ReactionInterpretationCandidateSource = Tuple[
    Mapping[str, Any],
    Mapping[str, ReactionSiteReference],
]


@dataclass(frozen=True)
class ReactionInterpretationBuild:
    """Interpretation annotations derived from a structural observation."""

    candidate_sources: Tuple[ReactionInterpretationCandidateSource, ...]
    candidates: Tuple[ReactionInterpretationCandidate, ...]
    selected_candidate: ReactionInterpretationCandidate | None
    selected_events: Tuple[ReactionInterpretationCandidate, ...]
    evidence_quality: str
    compatible_named_families: Tuple[str, ...]
    named_family: str | None
    product_contradicted_candidates: bool
    warnings: Tuple[str, ...]


def build_interpreted_partners(
    observation: ReactionObservation,
    selected: ReactionInterpretationCandidate | None,
    selected_events: Tuple[ReactionInterpretationCandidate, ...] = (),
) -> Tuple[ReactionPartner, ...]:
    """Overlay optional semantic roles on graph-derived site environments."""
    candidates = selected_events or ((selected,) if selected is not None else ())
    partners = []
    seen = set()
    for candidate in candidates:
        confidence = (
            1.0
            if candidate.verification
            in {
                "exact_product_reconstruction",
                "exact_multi_event_reconstruction",
            }
            else 0.7
        )
        for role, reference in sorted(candidate.role_assignments.items()):
            key = (str(role), int(reference.component_index), str(reference.site_id))
            if key in seen:
                continue
            seen.add(key)
            component = next(
                (
                    item
                    for item in observation.reactants
                    if item.component_index == reference.component_index
                ),
                None,
            )
            environment = next(
                (
                    item
                    for item in (
                        component.compound_analysis.site_environments
                        if component is not None
                        else ()
                    )
                    if item.site_id == reference.site_id
                ),
                None,
            )
            handles = tuple(
                sorted(
                    {
                        str(value)
                        for value in (
                            reference.details.get("handle_token"),
                            reference.details.get("center_token"),
                        )
                        if value
                    }
                )
            )
            contexts = tuple(
                sorted(
                    {
                        str(value)
                        for value in (
                            [reference.details.get("anchor_context")]
                            + list(reference.details.get("contexts") or ())
                        )
                        if value
                    }
                )
            )
            partners.append(
                ReactionPartner(
                    partner_id=(
                        f"RPI1:{reference.component_index}:"
                        f"{reference.site_id}:{role}"
                    ),
                    component_index=reference.component_index,
                    role=str(role),
                    role_confidence=confidence,
                    reactive_site_ids=(reference.site_id,),
                    handle_tokens=handles,
                    anchor_contexts=contexts,
                    chemist_label=reference.chemist_label,
                    nearby_groups=(environment.nearby_groups if environment else ()),
                    spectator_group_ids=tuple(
                        sorted(
                            group.group_id
                            for group in observation.spectator_groups
                            if group.component_index == reference.component_index
                        )
                    ),
                    reactivity_profile=(
                        environment.reactivity_profile if environment else None
                    ),
                )
            )
    return tuple(
        sorted(
            partners,
            key=lambda partner: (
                partner.component_index,
                partner.role or "",
                partner.partner_id,
            ),
        )
    )


def _semantic_assignment(
    annotation: Mapping[str, Any],
    reconstruction: ReactionReconstructionCandidate,
) -> dict[str, ReactionSiteReference]:
    return {
        str(role): reconstruction.slot_assignments[str(slot)]
        for role, slot in (annotation.get("role_bindings") or {}).items()
    }


def _annotate_candidate(
    annotation: Mapping[str, Any],
    reconstruction: ReactionReconstructionCandidate,
    *,
    label_style: str,
) -> tuple[
    ReactionInterpretationCandidate,
    ReactionInterpretationCandidateSource,
]:
    assignment = _semantic_assignment(annotation, reconstruction)
    semantic_by_slot = {
        str(slot): str(role)
        for role, slot in (annotation.get("role_bindings") or {}).items()
    }

    def semantic_path(path: str | None) -> str | None:
        if path is None:
            return None
        slot, separator, atom_role = path.partition(".")
        semantic_role = semantic_by_slot.get(slot, slot)
        return (
            f"{semantic_role}.{atom_role}"
            if separator
            else semantic_role
        )

    predicted_bond_changes = tuple(
        replace(
            change,
            atom_1_role=str(semantic_path(change.atom_1_role)),
            atom_2_role=semantic_path(change.atom_2_role),
        )
        for change in reconstruction.predicted_bond_changes
    )
    predicted_stereo_changes = tuple(
        replace(
            change,
            atom_1_role=str(semantic_path(change.atom_1_role)),
            atom_2_role=semantic_path(change.atom_2_role),
        )
        for change in reconstruction.predicted_stereo_changes
    )
    candidate = ReactionInterpretationCandidate(
        annotation_id=str(annotation["id"]),
        rewrite_outcome_id=reconstruction.rewrite_outcome_id,
        edit_archetype=reconstruction.edit_archetype,
        transformation_class=str(annotation["transformation_class"]),
        role_assignments=assignment,
        predicted_bond_changes=predicted_bond_changes,
        predicted_product_smiles=reconstruction.predicted_product_smiles,
        verification=reconstruction.verification,
        interpretation_label=render_reaction_label(
            dict(annotation), assignment, style=label_style
        ),
        predicted_stereo_changes=predicted_stereo_changes,
        compatible_named_families=tuple(
            annotation.get("compatible_named_families") or ()
        ),
        warnings=reconstruction.warnings,
    )
    return candidate, (annotation, assignment)


def build_reaction_interpretation_candidates(
    observation: ReactionObservation,
    *,
    label_style: str,
) -> ReactionInterpretationBuild:
    """Annotate lower-level reconstruction evidence; never reconstruct graphs."""
    annotations_by_rule: dict[str, list[Mapping[str, Any]]] = {}
    for annotation in load_reaction_interpretation_annotations():
        annotations_by_rule.setdefault(
            str(annotation["reconstruction_rule_id"]), []
        ).append(annotation)

    candidate_pairs = [
        _annotate_candidate(annotation, reconstruction, label_style=label_style)
        for reconstruction in observation.reconstruction_candidates
        for annotation in annotations_by_rule.get(reconstruction.rule_id, ())
    ]
    candidates = tuple(pair[0] for pair in candidate_pairs)
    sources = tuple(pair[1] for pair in candidate_pairs)

    def annotate_selected(
        reconstruction: ReactionReconstructionCandidate | None,
    ) -> ReactionInterpretationCandidate | None:
        if reconstruction is None:
            return None
        annotations = annotations_by_rule.get(reconstruction.rule_id, ())
        if not annotations:
            return None
        return _annotate_candidate(
            annotations[0], reconstruction, label_style=label_style
        )[0]

    selected = annotate_selected(observation.selected_reconstruction)
    selected_events = tuple(
        candidate
        for reconstruction in observation.selected_reconstruction_events
        for candidate in (annotate_selected(reconstruction),)
        if candidate is not None
    )
    compatible_families = (
        selected.compatible_named_families
        if selected is not None
        else tuple(
            sorted(
                {
                    family
                    for event in selected_events
                    for family in event.compatible_named_families
                }
            )
        )
    )
    named_family = (
        compatible_families[0]
        if selected is not None and len(compatible_families) == 1
        else None
    )
    product_contradicted = (
        selected is None
        and not selected_events
        and bool(observation.products)
        and bool(candidates)
        and all(
            candidate.verification in {"construction_failed", "product_mismatch"}
            for candidate in candidates
        )
    )
    evidence = (
        observation.selected_reconstruction.verification
        if selected is not None and observation.selected_reconstruction is not None
        else "exact_multi_event_reconstruction"
        if selected_events
        else "reactant_interpretation_only"
        if candidates
        else "unresolved"
    )
    warnings = ()
    if observation.selected_reconstruction is not None and selected is None:
        warnings = ("UNANNOTATED_SELECTED_RECONSTRUCTION",)
    return ReactionInterpretationBuild(
        candidate_sources=sources,
        candidates=candidates,
        selected_candidate=selected,
        selected_events=selected_events,
        evidence_quality=evidence,
        compatible_named_families=compatible_families,
        named_family=named_family,
        product_contradicted_candidates=product_contradicted,
        warnings=warnings,
    )


def _edit_key(edit: Any) -> tuple[Any, ...]:
    endpoints = []
    for atom in (edit.atom_1, edit.atom_2):
        if atom is None:
            endpoints.append(("H",))
        elif atom.atom_map_number is not None:
            endpoints.append(("map", atom.atom_map_number))
        else:
            endpoints.append(
                ("atom", atom.component_index, atom.atom_index, atom.element)
            )
    return (
        edit.edit_type,
        tuple(sorted(endpoints)),
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def interpret_reaction(
    observation: ReactionObservation,
    *,
    label_style: str = "unicode",
) -> ReactionInterpretation:
    """Add optional interpretation semantics to an immutable observation."""
    if not observation.valid:
        return ReactionInterpretation(
            evidence_quality="unresolved",
            warnings=("INVALID_REACTION_OBSERVATION",),
        )
    build = build_reaction_interpretation_candidates(
        observation, label_style=label_style
    )
    conflict = observation.evidence_quality in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }
    if build.selected_candidate is not None and observation.edits:
        predicted = normalize_predicted_edits(
            build.selected_candidate, observation.reactants
        )
        if predicted.valid:
            conflict = {_edit_key(edit) for edit in predicted.edits} != {
                _edit_key(edit) for edit in observation.edits
            }
    warnings = set(build.warnings)
    if conflict:
        warnings.add("INTERPRETATION_OBSERVATION_CONFLICT")
    selected = build.selected_candidate
    family_environment = (
        build_reaction_family_environment(
            observation.reactants,
            selected,
            observation.spectator_groups,
            build.evidence_quality,
        )
        if not conflict
        else None
    )
    product_connection = (
        build_product_connection(selected, build.evidence_quality, style=label_style)
        if not conflict
        else None
    )
    return ReactionInterpretation(
        candidates=build.candidates,
        selected_candidate=selected,
        selected_events=build.selected_events,
        partners=build_interpreted_partners(
            observation, selected, build.selected_events
        ),
        compatible_named_families=(
            () if conflict else build.compatible_named_families
        ),
        named_family=None if conflict else build.named_family,
        family_environment=family_environment,
        product_connection=product_connection,
        evidence_quality=(
            "conflicting_interpretation_evidence"
            if conflict
            else build.evidence_quality
        ),
        warnings=tuple(sorted(warnings)),
    )


__all__ = [
    "ReactionInterpretationCandidateSource",
    "ReactionInterpretationBuild",
    "build_interpreted_partners",
    "build_reaction_interpretation_candidates",
    "interpret_reaction",
]
