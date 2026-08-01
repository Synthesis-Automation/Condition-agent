"""Public reaction featurization API."""

from __future__ import annotations

from .labels import available_styles
from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_contextual_labels import build_contextual_transformation_label
from .reaction_edits import (
    EditNormalizationResult,
    normalize_mapped_edits,
    resolve_reaction_evidence,
)
from .reaction_labels import render_reactant_label
from .reaction_environments import build_reaction_family_environment
from .reaction_fallback_descriptors import build_reaction_fallback_descriptor
from .reaction_models import ReactionAnalysis, ReactionInterpretation
from .reaction_interpretation import (
    build_interpreted_partners,
    build_reaction_interpretation_candidates,
    canonical_without_maps,
)
from .reaction_observation import build_reaction_observation
from .reaction_rendering import render_reaction
from .reaction_parser import parse_reaction_smiles
from .partial_product_correspondence import (
    infer_partial_product_transformation,
    render_partial_product_transformation,
)
from .reaction_products import build_product_connection
from .reaction_signatures import build_observation_signature


def featurize_reaction(
    reaction_smiles: str,
    *,
    label_style: str = "unicode",
    max_candidates: int = 500,
    _mapped_edit_override: EditNormalizationResult | None = None,
    _mapped_provider: str = "supplied_atom_mapping",
) -> ReactionAnalysis:
    """Analyze a reaction with site grammars and product reconstruction."""
    parsed = parse_reaction_smiles(reaction_smiles)
    if label_style not in available_styles():
        return ReactionAnalysis(
            reaction_smiles, False, error=f"UNKNOWN_LABEL_STYLE:{label_style}"
        )
    if not parsed.valid:
        return ReactionAnalysis(
            reaction_smiles,
            False,
            parsed.reactants,
            parsed.agents,
            parsed.products,
            warnings=parsed.warnings,
            error=parsed.error,
        )
    observed_products = {
        canonical_without_maps(component.input_smiles)
        for component in parsed.products
    }
    observed_products.discard(None)
    supplied_mapping = (
        _mapped_edit_override
        if _mapped_edit_override is not None
        else normalize_mapped_edits(parsed.reactants, parsed.products)
    )
    invalid_supplied_mapping = supplied_mapping.evidence == "invalid_atom_mapping"
    interpretation_build = build_reaction_interpretation_candidates(
        reactants=parsed.reactants,
        observed_products={str(product) for product in observed_products},
        invalid_supplied_mapping=invalid_supplied_mapping,
        label_style=label_style,
        max_candidates=max_candidates,
    )
    raw = interpretation_build.raw_candidates
    candidates = interpretation_build.candidates
    candidate_sources = interpretation_build.candidate_sources
    selected = interpretation_build.selected_candidate
    selected_events = interpretation_build.selected_events
    evidence = interpretation_build.evidence_quality
    warnings = list(parsed.warnings) + list(interpretation_build.warnings)
    mapped_changes = tuple(supplied_map_bond_changes(reaction_smiles))
    named_family = interpretation_build.named_family
    compatible_named_families = (
        interpretation_build.compatible_named_families
    )
    edit_result = resolve_reaction_evidence(
        parsed.reactants,
        parsed.products,
        selected,
        selected_events,
        tuple(candidates),
        mapped_override=_mapped_edit_override,
        mapped_provider=_mapped_provider,
    )
    interpretation_conflict = edit_result.evidence in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }
    if interpretation_conflict:
        named_family = None
        compatible_named_families = ()
        warnings.append("INTERPRETATION_OBSERVATION_CONFLICT")
    observation = build_reaction_observation(
        input_reaction_smiles=reaction_smiles,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        edit_result=edit_result,
        mapped_bond_changes=mapped_changes,
        operator_candidates=raw,
        selected_operator_candidate=selected,
        selected_operator_events=selected_events,
        warnings=warnings,
    )
    spectators = observation.spectator_groups
    reaction_topology = observation.topology
    reaction_core = observation.core
    reaction_completeness = observation.completeness
    assert reaction_completeness is not None
    warnings = list(observation.warnings)
    family_environment = (
        build_reaction_family_environment(
            parsed.reactants,
            selected,
            spectators,
            evidence,
        )
        if not interpretation_conflict
        else None
    )
    product_connection = (
        build_product_connection(
            selected,
            evidence,
            style=label_style,
        )
        if not interpretation_conflict
        else None
    )
    product_contradicted_candidates = (
        interpretation_build.product_contradicted_candidates
    )
    partial_product_transformation = (
        infer_partial_product_transformation(
            reactants=parsed.reactants,
            agents=parsed.agents,
            products=parsed.products,
            completeness=reaction_completeness,
        )
        if (
            selected is None
            and not selected_events
            and not edit_result.edits
            and not invalid_supplied_mapping
        )
        else None
    )
    if partial_product_transformation is not None:
        warnings.extend(partial_product_transformation.warnings)
    if product_contradicted_candidates:
        warnings.append("PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES")
    contextual_label = (
        None
        if edit_result.evidence
        in {
            "conflicting_edit_evidence",
            "conflicting_stereochemical_evidence",
        }
        else build_contextual_transformation_label(
            parsed.reactants, edit_result.edits, style=label_style
        )
    )
    edit_hypotheses = observation.edit_hypotheses
    evidence_candidates = observation.evidence_candidates
    effective_evidence = evidence
    if edit_result.evidence in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }:
        effective_evidence = edit_result.evidence
    elif edit_result.evidence.startswith("validated_mapping"):
        effective_evidence = edit_result.evidence
    elif edit_result.evidence.startswith("external_mapping"):
        effective_evidence = edit_result.evidence
    elif selected is None and (
        edit_result.valid or edit_result.evidence == "ambiguous_atom_correspondence"
    ):
        effective_evidence = edit_result.evidence
    if partial_product_transformation is not None:
        effective_evidence = partial_product_transformation.evidence
    display_arrow = "→" if label_style == "unicode" else "->"
    selected_product_label = None
    if (
        selected is not None
        and selected.verification == "exact_product_reconstruction"
        and selected.grammar_label
        and display_arrow in selected.grammar_label
    ):
        selected_product_label = selected.grammar_label.split(display_arrow, 1)[
            1
        ].strip()
    reaction_signature = (
        build_observation_signature(
            observation,
            contextual_product_label=(
                contextual_label.after
                if contextual_label is not None
                else selected_product_label
            ),
            warnings=warnings,
        )
        if (
            reaction_topology is not None
            and reaction_completeness.status != "incomplete"
        )
        else None
    )
    fallback_label = selected.grammar_label if selected else None
    fallback_status = "exact_product" if selected else "unavailable"
    fallback_detailed_label = None
    if partial_product_transformation is not None:
        fallback_label, fallback_detailed_label = render_partial_product_transformation(
            partial_product_transformation,
            reactants=parsed.reactants,
            products=parsed.products,
            style=label_style,
        )
        fallback_status = "partial_product_correspondence"
    elif selected is None and candidates:
        exact_candidate_indices = tuple(
            index
            for index, candidate in enumerate(candidates)
            if candidate.verification == "exact_product_reconstruction"
        )
        fallback_raw = (
            [candidate_sources[index] for index in exact_candidate_indices]
            if exact_candidate_indices
            else candidate_sources
        )
        if exact_candidate_indices and len(exact_candidate_indices) < len(candidates):
            warnings.append("PRODUCT_MISMATCH_CANDIDATES_EXCLUDED_FROM_LABEL")
        if any(len(assignment) > 1 for _, assignment in fallback_raw):
            fallback_raw = [
                (grammar, assignment)
                for grammar, assignment in fallback_raw
                if len(assignment) > 1
            ]
        reactant_labels = sorted(
            {
                render_reactant_label(assignment, style=label_style)
                for _, assignment in fallback_raw
            }
        )
        if len(reactant_labels) == 1:
            fallback_label = f"{reactant_labels[0]} {display_arrow}"
            fallback_status = (
                "product_contradicted_reactants"
                if product_contradicted_candidates
                else "reactant_only"
            )
        elif reactant_labels:
            fallback_label = (
                " OR ".join(f"({label})" for label in reactant_labels)
                + f" {display_arrow}"
            )
            fallback_status = (
                "product_contradicted_reactants"
                if product_contradicted_candidates
                else "ambiguous_reactants"
            )
    elif selected is not None and product_connection is not None:
        reactants_label = render_reactant_label(
            selected.role_assignments, style=label_style
        )
        fallback_label = (
            f"{reactants_label} {display_arrow} {product_connection.concise_label}"
        )
    interpretation_warnings = set(interpretation_build.warnings)
    if product_contradicted_candidates:
        interpretation_warnings.add(
            "PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES"
        )
    if interpretation_conflict:
        interpretation_warnings.add("INTERPRETATION_OBSERVATION_CONFLICT")
    interpretation = ReactionInterpretation(
        candidates=candidates,
        selected_candidate=selected,
        selected_events=selected_events,
        partners=build_interpreted_partners(
            observation,
            selected,
            selected_events,
        ),
        compatible_named_families=compatible_named_families,
        named_family=named_family,
        family_environment=family_environment,
        product_connection=product_connection,
        evidence_quality=(
            "conflicting_interpretation_evidence"
            if interpretation_conflict
            else evidence
        ),
        warnings=tuple(sorted(interpretation_warnings)),
    )
    reaction_label = render_reaction(
        observation,
        interpretation,
        style=label_style,
        fallback_label=fallback_label,
        fallback_status=fallback_status,
        fallback_detailed_label=fallback_detailed_label,
        fallback_evidence=(
            partial_product_transformation.evidence
            if partial_product_transformation is not None and not edit_result.edits
            else None
        ),
        fallback_confidence=(
            partial_product_transformation.confidence
            if partial_product_transformation is not None and not edit_result.edits
            else None
        ),
    )
    fallback_descriptor = build_reaction_fallback_descriptor(
        reactants=parsed.reactants,
        products=parsed.products,
        candidates=tuple(candidates),
        signature=reaction_signature,
        partial_transformation=partial_product_transformation,
        completeness=reaction_completeness,
        evidence_quality=effective_evidence,
        warnings=warnings,
    )
    return ReactionAnalysis(
        input_reaction_smiles=reaction_smiles,
        valid=True,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        candidates=tuple(candidates),
        selected_candidate=selected,
        selected_events=selected_events,
        edit_archetype=(
            selected.edit_archetype
            if selected is not None
            else reaction_signature.edit_archetype
            if reaction_signature is not None
            else "unresolved"
        ),
        transformation_class=(
            selected.transformation_class
            if selected
            else reaction_signature.transformation_class
            if reaction_signature is not None
            else partial_product_transformation.transformation_class
            if partial_product_transformation is not None
            else None
        ),
        compatible_named_families=compatible_named_families,
        named_family=named_family,
        reaction_label=reaction_label,
        evidence_quality=effective_evidence,
        mapped_bond_changes=mapped_changes,
        spectator_groups=spectators,
        family_environment=family_environment,
        product_connection=product_connection,
        reaction_topology=reaction_topology,
        evidence_candidates=evidence_candidates,
        edit_hypotheses=edit_hypotheses,
        reaction_signature=reaction_signature,
        fallback_descriptor=fallback_descriptor,
        partial_product_transformation=partial_product_transformation,
        reaction_completeness=reaction_completeness,
        reaction_core=reaction_core,
        observation=observation,
        interpretation=interpretation,
        warnings=tuple(sorted(set(warnings))),
    )


__all__ = ["featurize_reaction"]
