"""Public reaction featurization API."""

from __future__ import annotations

from dataclasses import replace

from .labels import available_styles
from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_edits import (
    EditNormalizationResult,
    normalize_mapped_edits,
    resolve_structural_evidence,
)
from .reaction_fallback_descriptors import build_reaction_fallback_descriptor
from .reaction_models import ReactionAnalysis, ReactionInterpretation
from .reaction_interpretation import interpret_reaction
from .reaction_observation import build_reaction_observation
from .reaction_spectators import derive_observed_spectator_groups
from .reaction_rendering import render_reaction
from .reaction_render_context import build_reaction_render_context
from .reaction_r_group_context import build_r_group_functional_contexts
from .reaction_parser import interpret_parsed_molecules, parse_reaction_smiles
from .reaction_stoichiometry import infer_reactant_multiplicity
from .partial_product_correspondence import (
    infer_partial_product_transformation,
)
from .reaction_signatures import build_observation_signature


def featurize_reaction(
    reaction_smiles: str,
    *,
    label_style: str = "unicode",
    _mapped_edit_override: EditNormalizationResult | None = None,
    _mapped_provider: str = "supplied_atom_mapping",
) -> ReactionAnalysis:
    """Analyze structural evidence before optional chemistry interpretation."""
    parsed = parse_reaction_smiles(
        reaction_smiles,
        include_molecular_interpretation=False,
    )
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
    multiplicity_warnings: tuple[str, ...] = ()
    if _mapped_edit_override is None:
        multiplicity = infer_reactant_multiplicity(
            parsed.reactants,
            parsed.products,
        )
        multiplicity_warnings = multiplicity.warnings
        if multiplicity.inferred_copy_count:
            parsed = replace(
                parsed,
                reactants=multiplicity.reactants,
                warnings=tuple(
                    sorted(set(parsed.warnings + multiplicity.warnings))
                ),
            )
    supplied_mapping = (
        _mapped_edit_override
        if _mapped_edit_override is not None
        else normalize_mapped_edits(parsed.reactants, parsed.products)
    )
    invalid_supplied_mapping = supplied_mapping.evidence == "invalid_atom_mapping"
    warnings = list(parsed.warnings)
    mapped_changes = tuple(supplied_map_bond_changes(reaction_smiles))
    edit_result = resolve_structural_evidence(
        parsed.reactants,
        parsed.products,
        mapped_override=_mapped_edit_override,
        mapped_provider=_mapped_provider,
        additional_warnings=multiplicity_warnings,
    )
    observation = build_reaction_observation(
        input_reaction_smiles=reaction_smiles,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        edit_result=edit_result,
        mapped_bond_changes=mapped_changes,
        warnings=warnings,
    )
    reaction_topology = observation.topology
    reaction_core = observation.core
    reaction_completeness = observation.completeness
    assert reaction_completeness is not None
    warnings = list(observation.warnings)
    reaction_signature = (
        build_observation_signature(
            observation,
            warnings=warnings,
        )
        if (
            reaction_topology is not None
            and reaction_completeness.status != "incomplete"
        )
        else None
    )

    # Optional molecular annotations start only after the observation-only
    # signature has been finalized.
    parsed = interpret_parsed_molecules(parsed, label_style=label_style)
    spectators = derive_observed_spectator_groups(
        parsed.reactants,
        observation.edits,
        observation.evidence_quality,
    )
    r_group_functional_contexts = build_r_group_functional_contexts(
        reactants=parsed.reactants,
        core=reaction_core,
        spectator_groups=spectators,
    )
    partial_product_transformation = (
        infer_partial_product_transformation(
            reactants=parsed.reactants,
            agents=parsed.agents,
            products=parsed.products,
            completeness=reaction_completeness,
        )
        if (
            not edit_result.edits
            and not invalid_supplied_mapping
        )
        else None
    )
    if partial_product_transformation is not None:
        warnings.extend(partial_product_transformation.warnings)

    # Optional interpretation starts only after the generic signature exists.
    interpretation = replace(
        interpret_reaction(observation, label_style=label_style),
        spectator_groups=spectators,
        r_group_functional_contexts=r_group_functional_contexts,
    )
    edit_hypotheses = observation.edit_hypotheses
    evidence_candidates = observation.evidence_candidates
    effective_evidence = (
        partial_product_transformation.evidence
        if partial_product_transformation is not None
        else edit_result.evidence
    )
    render_context = build_reaction_render_context(
        observation=observation,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        signature=reaction_signature,
        partial_transformation=partial_product_transformation,
        style=label_style,
        interpretation=interpretation,
    )
    reaction_label = render_reaction(render_context)
    fallback_descriptor = build_reaction_fallback_descriptor(
        reactants=parsed.reactants,
        products=parsed.products,
        signature=reaction_signature,
        partial_transformation=partial_product_transformation,
        completeness=reaction_completeness,
        evidence_quality=effective_evidence,
        has_edit_hypotheses=bool(edit_hypotheses),
        warnings=warnings,
    )
    return ReactionAnalysis(
        input_reaction_smiles=reaction_smiles,
        valid=True,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        edit_archetype=(
            reaction_signature.edit_archetype
            if reaction_signature is not None
            else "unresolved"
        ),
        transformation_class=(
            reaction_signature.transformation_class
            if reaction_signature is not None
            else partial_product_transformation.transformation_class
            if partial_product_transformation is not None
            else None
        ),
        compatible_named_families=interpretation.compatible_named_families,
        named_family=interpretation.named_family,
        reaction_label=reaction_label,
        evidence_quality=effective_evidence,
        mapped_bond_changes=mapped_changes,
        spectator_groups=spectators,
        family_environment=interpretation.family_environment,
        product_connection=interpretation.product_connection,
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


def identify_reaction_patterns(
    reaction_smiles: str,
    *,
    label_style: str = "unicode",
) -> ReactionInterpretation | None:
    """Identify optional structural patterns directly from reaction SMILES."""
    return featurize_reaction(
        reaction_smiles,
        label_style=label_style,
    ).interpretation


__all__ = ["featurize_reaction", "identify_reaction_patterns"]
