"""Optional grammar and family interpretation of reaction observations.

This module owns grammar enumeration and operator-backed product
reconstruction.  It returns interpretations and edit-evidence proposals; it
does not build generic topology, reaction cores, or signature identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Mapping, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .connectivity_rewrite import apply_connectivity_rewrite
from .reaction_archetypes import infer_bond_change_archetype
from .reaction_candidates import enumerate_reaction_candidates
from .reaction_labels import render_reaction_label
from .reaction_edits import normalize_predicted_edits
from .reaction_environments import build_reaction_family_environment
from .reaction_models import (
    ReactionCandidate,
    ReactionComponent,
    ReactionInterpretation,
    ReactionObservation,
    ReactionPartner,
    ReactionSiteReference,
)
from .reaction_products import build_product_connection
from .reaction_multi_events import (
    apply_rewrite_sequence,
    equivalent_multi_event_interpretations,
    exact_multi_event_reconstructions,
)


RawReactionCandidate = Tuple[
    Mapping[str, Any],
    Mapping[str, ReactionSiteReference],
]


@dataclass(frozen=True)
class ReactionInterpretationBuild:
    """Internal build result retaining inputs needed by later reconciliation."""

    raw_candidates: Tuple[RawReactionCandidate, ...]
    candidate_sources: Tuple[RawReactionCandidate, ...]
    candidates: Tuple[ReactionCandidate, ...]
    selected_candidate: ReactionCandidate | None
    selected_events: Tuple[ReactionCandidate, ...]
    evidence_quality: str
    compatible_named_families: Tuple[str, ...]
    named_family: str | None
    product_contradicted_candidates: bool
    warnings: Tuple[str, ...]


def build_interpreted_partners(
    observation: ReactionObservation,
    selected: ReactionCandidate | None,
    selected_events: Tuple[ReactionCandidate, ...] = (),
) -> Tuple[ReactionPartner, ...]:
    """Overlay optional grammar roles on graph-derived site environments."""
    candidates = selected_events or ((selected,) if selected is not None else ())
    partners = []
    seen = set()
    for candidate in candidates:
        confidence = 1.0 if candidate.verification == "exact_product_reconstruction" else 0.7
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
                    nearby_groups=(
                        environment.nearby_groups if environment else ()
                    ),
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


def canonical_without_maps(smiles: str) -> str | None:
    """Return canonical isomeric SMILES after removing map numbers."""
    from rdkit import Chem

    molecule = parse_smiles(smiles)
    if molecule is None:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        )
    except Exception:
        return None


def _candidate(
    *,
    grammar: Mapping[str, Any],
    assignment: Mapping[str, ReactionSiteReference],
    outcome: Any,
    observed_products: set[str],
    label_style: str,
) -> ReactionCandidate:
    predicted_canonical = (
        canonical_without_maps(outcome.predicted_product_smiles)
        if outcome.predicted_product_smiles
        else None
    )
    verification = (
        "exact_product_reconstruction"
        if predicted_canonical in observed_products
        else "product_mismatch"
        if predicted_canonical
        else "construction_failed"
    )
    return ReactionCandidate(
        grammar_id=str(grammar["id"]),
        rewrite_outcome_id=outcome.outcome_id,
        edit_archetype=infer_bond_change_archetype(
            outcome.predicted_bond_changes
        ),
        transformation_class=str(grammar["transformation_class"]),
        role_assignments=dict(assignment),
        predicted_bond_changes=outcome.predicted_bond_changes,
        predicted_product_smiles=predicted_canonical,
        verification=verification,
        grammar_label=render_reaction_label(
            dict(grammar),
            dict(assignment),
            style=label_style,
        ),
        predicted_stereo_changes=outcome.predicted_stereo_changes,
        compatible_named_families=tuple(
            grammar.get("compatible_named_families") or ()
        ),
        warnings=outcome.warnings,
    )


def _equivalent_exact_candidates(
    candidates: Tuple[ReactionCandidate, ...],
) -> bool:
    signatures = {
        (
            candidate.grammar_id,
            tuple(
                sorted(
                    site.canonical_signature
                    for site in candidate.role_assignments.values()
                )
            ),
        )
        for candidate in candidates
    }
    return len(signatures) == 1


def build_reaction_interpretation_candidates(
    *,
    reactants: Tuple[ReactionComponent, ...],
    observed_products: set[str],
    invalid_supplied_mapping: bool,
    label_style: str,
    max_candidates: int,
) -> ReactionInterpretationBuild:
    """Build optional grammar candidates without constructing generic facts."""
    raw_values = list(enumerate_reaction_candidates(reactants))
    warnings: List[str] = []
    if len(raw_values) > max_candidates:
        raw_values = raw_values[:max_candidates]
        warnings.append("CANDIDATE_LIMIT_REACHED")
    raw = tuple((grammar, assignment) for grammar, assignment in raw_values)
    candidates: List[ReactionCandidate] = []
    candidate_sources: List[RawReactionCandidate] = []
    exact: List[ReactionCandidate] = []
    for grammar, assignment in raw:
        outcomes = apply_connectivity_rewrite(grammar, assignment, reactants)
        for outcome in outcomes:
            candidate = _candidate(
                grammar=grammar,
                assignment=assignment,
                outcome=outcome,
                observed_products=observed_products,
                label_style=label_style,
            )
            candidates.append(candidate)
            candidate_sources.append((grammar, assignment))
            if candidate.verification == "exact_product_reconstruction":
                exact.append(candidate)
            if len(candidates) >= max_candidates:
                warnings.append("CANDIDATE_LIMIT_REACHED")
                break
        if len(candidates) >= max_candidates:
            break

    selected = None
    selected_events: Tuple[ReactionCandidate, ...] = ()
    evidence = "reactant_grammar_only" if candidates else "unresolved"
    if invalid_supplied_mapping:
        evidence = "unresolved"
    elif len(exact) == 1:
        selected, evidence = exact[0], "exact_product_reconstruction"
    elif len(exact) > 1:
        exact_tuple = tuple(exact)
        if _equivalent_exact_candidates(exact_tuple):
            selected, evidence = exact[0], "exact_product_reconstruction"
            warnings.append("SYMMETRY_EQUIVALENT_ASSIGNMENTS")
        else:
            evidence = "ambiguous"
            warnings.append("AMBIGUOUS_PARTICIPATING_SITES")
    elif candidates and not invalid_supplied_mapping:
        multi_exact = exact_multi_event_reconstructions(
            raw,
            reactants,
            observed_products,
        )
        if multi_exact and equivalent_multi_event_interpretations(multi_exact):
            chosen = multi_exact[0]
            composite_product = apply_rewrite_sequence(chosen, reactants)
            event_candidates = []
            for grammar, assignment in chosen:
                outcomes = apply_connectivity_rewrite(
                    grammar,
                    assignment,
                    reactants,
                )
                if not outcomes:
                    continue
                outcome = outcomes[0]
                event_candidates.append(
                    ReactionCandidate(
                        grammar_id=str(grammar["id"]),
                        rewrite_outcome_id=outcome.outcome_id,
                        edit_archetype=infer_bond_change_archetype(
                            outcome.predicted_bond_changes
                        ),
                        transformation_class=str(
                            grammar["transformation_class"]
                        ),
                        role_assignments=dict(assignment),
                        predicted_bond_changes=outcome.predicted_bond_changes,
                        predicted_product_smiles=composite_product,
                        verification="exact_multi_event_reconstruction",
                        grammar_label=render_reaction_label(
                            dict(grammar),
                            dict(assignment),
                            style=label_style,
                        ),
                        predicted_stereo_changes=(
                            outcome.predicted_stereo_changes
                        ),
                        compatible_named_families=tuple(
                            grammar.get("compatible_named_families") or ()
                        ),
                    )
                )
            selected_events = tuple(event_candidates)
            evidence = "exact_multi_event_reconstruction"
            if len(multi_exact) > 1:
                warnings.append(
                    "SYMMETRY_EQUIVALENT_MULTI_EVENT_ASSIGNMENTS"
                )
        elif multi_exact:
            evidence = "ambiguous"
            warnings.append("AMBIGUOUS_MULTI_EVENT_ASSIGNMENTS")

    candidate_tuple = tuple(candidates)
    product_contradicted = (
        selected is None
        and not selected_events
        and bool(observed_products)
        and bool(candidate_tuple)
        and all(
            candidate.verification in {
                "construction_failed",
                "product_mismatch",
            }
            for candidate in candidate_tuple
        )
    )
    named_family = (
        selected.compatible_named_families[0]
        if selected and len(selected.compatible_named_families) == 1
        else None
    )
    compatible_families = (
        selected.compatible_named_families
        if selected
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
    return ReactionInterpretationBuild(
        raw_candidates=raw,
        candidate_sources=tuple(candidate_sources),
        candidates=candidate_tuple,
        selected_candidate=selected,
        selected_events=selected_events,
        evidence_quality=evidence,
        compatible_named_families=compatible_families,
        named_family=named_family,
        product_contradicted_candidates=product_contradicted,
        warnings=tuple(sorted(set(warnings))),
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
                (
                    "atom",
                    atom.component_index,
                    atom.atom_index,
                    atom.element,
                )
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
    max_candidates: int = 500,
) -> ReactionInterpretation:
    """Interpret an immutable observation with optional grammar semantics."""
    if not observation.valid:
        return ReactionInterpretation(
            evidence_quality="unresolved",
            warnings=("INVALID_REACTION_OBSERVATION",),
        )
    observed_products = {
        canonical
        for component in observation.products
        for canonical in (canonical_without_maps(component.input_smiles),)
        if canonical is not None
    }
    invalid_mapping = any(
        candidate.evidence == "invalid_atom_mapping"
        for candidate in observation.evidence_candidates
    )
    build = build_reaction_interpretation_candidates(
        reactants=observation.reactants,
        observed_products=observed_products,
        invalid_supplied_mapping=invalid_mapping,
        label_style=label_style,
        max_candidates=max_candidates,
    )
    conflict = observation.evidence_quality in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }
    if build.selected_candidate is not None and observation.edits:
        predicted = normalize_predicted_edits(
            build.selected_candidate,
            observation.reactants,
        )
        if predicted.valid:
            conflict = {
                _edit_key(edit) for edit in predicted.edits
            } != {_edit_key(edit) for edit in observation.edits}
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
        build_product_connection(
            selected,
            build.evidence_quality,
            style=label_style,
        )
        if not conflict
        else None
    )
    return ReactionInterpretation(
        candidates=build.candidates,
        selected_candidate=selected,
        selected_events=build.selected_events,
        partners=build_interpreted_partners(
            observation,
            selected,
            build.selected_events,
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
    "RawReactionCandidate",
    "ReactionInterpretationBuild",
    "build_interpreted_partners",
    "build_reaction_interpretation_candidates",
    "canonical_without_maps",
    "interpret_reaction",
]
