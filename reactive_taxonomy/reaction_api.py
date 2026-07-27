"""Public reaction featurization API."""

from __future__ import annotations

from typing import List

from .chemistry.rdkit_utils import parse_smiles

from .labels import available_styles
from .reaction_bond_changes import supplied_map_bond_changes
from .reaction_candidates import enumerate_reaction_candidates
from .reaction_completeness import build_reaction_completeness
from .reaction_contextual_labels import build_contextual_transformation_label
from .reaction_display_labels import build_reaction_display_label
from .reaction_edits import normalize_mapped_edits, normalize_reaction_edits
from .reaction_labels import render_reactant_label, render_reaction_label
from .reaction_environments import build_reaction_family_environment
from .reaction_models import ReactionAnalysis, ReactionCandidate
from .reaction_multi_events import (
    equivalent_multi_event_interpretations,
    exact_multi_event_reconstructions,
)
from .reaction_operators import apply_operator, apply_operator_sequence
from .reaction_parser import parse_reaction_smiles
from .reaction_products import build_product_connection
from .reaction_spectators import derive_spectator_groups
from .reaction_signatures import build_reaction_signature
from .reaction_topology import build_reaction_topology


def _canonical_without_maps(smiles: str) -> str | None:
    from rdkit import Chem

    mol = parse_smiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def featurize_reaction(
    reaction_smiles: str,
    *,
    label_style: str = "unicode",
    max_candidates: int = 500,
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
        _canonical_without_maps(component.input_smiles) for component in parsed.products
    }
    observed_products.discard(None)
    raw = enumerate_reaction_candidates(parsed.reactants)
    warnings = list(parsed.warnings)
    supplied_mapping = normalize_mapped_edits(parsed.reactants, parsed.products)
    invalid_supplied_mapping = supplied_mapping.evidence == "invalid_atom_mapping"
    if len(raw) > max_candidates:
        raw = raw[:max_candidates]
        warnings.append("CANDIDATE_LIMIT_REACHED")
    candidates: List[ReactionCandidate] = []
    exact: List[ReactionCandidate] = []
    for grammar, assignment in raw:
        predicted, changes = apply_operator(grammar, assignment, parsed.reactants)
        predicted_canonical = _canonical_without_maps(predicted) if predicted else None
        verification = (
            "exact_product_reconstruction"
            if predicted_canonical in observed_products
            else ("product_mismatch" if predicted else "construction_failed")
        )
        label = render_reaction_label(grammar, assignment, style=label_style)
        candidate = ReactionCandidate(
            grammar_id=grammar["id"],
            transformation_class=grammar["transformation_class"],
            role_assignments=assignment,
            predicted_bond_changes=changes,
            predicted_product_smiles=predicted_canonical,
            verification=verification,
            reaction_label=label,
            compatible_named_families=tuple(
                grammar.get("compatible_named_families") or []
            ),
        )
        candidates.append(candidate)
        if verification == "exact_product_reconstruction":
            exact.append(candidate)
    selected = None
    selected_events: tuple[ReactionCandidate, ...] = ()
    evidence = "reactant_grammar_only" if candidates else "unresolved"
    if invalid_supplied_mapping:
        evidence = "unresolved"
    elif len(exact) == 1:
        selected, evidence = exact[0], "exact_product_reconstruction"
    elif len(exact) > 1:
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
            for candidate in exact
        }
        if len(signatures) == 1:
            selected, evidence = exact[0], "exact_product_reconstruction"
            warnings.append("SYMMETRY_EQUIVALENT_ASSIGNMENTS")
        else:
            evidence = "ambiguous"
            warnings.append("AMBIGUOUS_PARTICIPATING_SITES")
    elif candidates and not invalid_supplied_mapping:
        multi_exact = exact_multi_event_reconstructions(
            raw,
            parsed.reactants,
            {str(product) for product in observed_products},
        )
        if multi_exact and equivalent_multi_event_interpretations(multi_exact):
            chosen = multi_exact[0]
            composite_product = apply_operator_sequence(chosen, parsed.reactants)
            selected_events = tuple(
                ReactionCandidate(
                    grammar_id=str(grammar["id"]),
                    transformation_class=str(grammar["transformation_class"]),
                    role_assignments=assignment,
                    predicted_bond_changes=apply_operator(
                        grammar, assignment, parsed.reactants
                    )[1],
                    predicted_product_smiles=composite_product,
                    verification="exact_multi_event_reconstruction",
                    reaction_label=render_reaction_label(
                        grammar, assignment, style=label_style
                    ),
                    compatible_named_families=tuple(
                        grammar.get("compatible_named_families") or []
                    ),
                )
                for grammar, assignment in chosen
            )
            evidence = "exact_multi_event_reconstruction"
            if len(multi_exact) > 1:
                warnings.append("SYMMETRY_EQUIVALENT_MULTI_EVENT_ASSIGNMENTS")
        elif multi_exact:
            evidence = "ambiguous"
            warnings.append("AMBIGUOUS_MULTI_EVENT_ASSIGNMENTS")
    mapped_changes = tuple(supplied_map_bond_changes(reaction_smiles))
    spectators = derive_spectator_groups(
        parsed.reactants,
        selected,
        evidence,
        selected_events,
    )
    family_environment = build_reaction_family_environment(
        parsed.reactants, selected, spectators, evidence
    )
    product_connection = build_product_connection(selected, evidence, style=label_style)
    named_family = (
        selected.compatible_named_families[0]
        if selected and len(selected.compatible_named_families) == 1
        else None
    )
    compatible_named_families = (
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
    edit_result = normalize_reaction_edits(
        parsed.reactants,
        parsed.products,
        selected,
        selected_events,
    )
    reaction_completeness = build_reaction_completeness(
        reactants=parsed.reactants,
        products=parsed.products,
        raw_candidates=raw,
        selected=selected,
        selected_events=selected_events,
        edit_result=edit_result,
    )
    warnings.extend(reaction_completeness.warnings)
    contextual_label = (
        None
        if edit_result.evidence == "conflicting_edit_evidence"
        else build_contextual_transformation_label(
            parsed.reactants, edit_result.edits, style=label_style
        )
    )
    warnings.extend(edit_result.warnings)
    reaction_topology = build_reaction_topology(
        reactants=parsed.reactants,
        products=parsed.products,
        selected=selected,
        edit_result=edit_result,
    )
    effective_evidence = evidence
    if edit_result.evidence == "conflicting_edit_evidence":
        effective_evidence = edit_result.evidence
    elif edit_result.evidence.startswith("validated_mapping"):
        effective_evidence = edit_result.evidence
    elif selected is None and (
        edit_result.valid or edit_result.evidence == "ambiguous_atom_correspondence"
    ):
        effective_evidence = edit_result.evidence
    display_arrow = "→" if label_style == "unicode" else "->"
    selected_product_label = None
    if (
        selected is not None
        and selected.verification == "exact_product_reconstruction"
        and selected.reaction_label
        and display_arrow in selected.reaction_label
    ):
        selected_product_label = selected.reaction_label.split(
            display_arrow, 1
        )[1].strip()
    reaction_signature = (
        build_reaction_signature(
            reactants=parsed.reactants,
            selected=selected,
            selected_events=selected_events,
            edit_result=edit_result,
            family_environment=family_environment,
            product_connection=product_connection,
            spectators=spectators,
            named_family=named_family,
            compatible_named_families=compatible_named_families,
            topology=reaction_topology,
            completeness=reaction_completeness,
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
    reaction_label = selected.reaction_label if selected else None
    reaction_label_status = "exact_product" if selected else "unavailable"
    if selected is None and candidates:
        fallback_raw = raw
        if any(len(assignment) > 1 for _, assignment in raw):
            fallback_raw = [
                (grammar, assignment)
                for grammar, assignment in raw
                if len(assignment) > 1
            ]
        reactant_labels = sorted(
            {
                render_reactant_label(assignment, style=label_style)
                for _, assignment in fallback_raw
            }
        )
        if len(reactant_labels) == 1:
            reaction_label = f"{reactant_labels[0]} {display_arrow}"
            reaction_label_status = "reactant_only"
        elif reactant_labels:
            reaction_label = (
                " OR ".join(f"({label})" for label in reactant_labels)
                + f" {display_arrow}"
            )
            reaction_label_status = "ambiguous_reactants"
    elif selected is not None and product_connection is not None:
        reactants_label = render_reactant_label(
            selected.role_assignments, style=label_style
        )
        reaction_label = (
            f"{reactants_label} {display_arrow} "
            f"{product_connection.concise_label}"
        )
    display_label = build_reaction_display_label(
        edits=edit_result.edits,
        selected_label=reaction_label if selected is not None else None,
        selected_exact=bool(
            selected and selected.verification == "exact_product_reconstruction"
        ),
        grammar_id=selected.grammar_id if selected is not None else None,
        contextual_label=contextual_label,
        named_family=named_family,
        fallback_label=reaction_label,
        fallback_status=reaction_label_status,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
        events=(reaction_signature.events if reaction_signature is not None else ()),
        topology=reaction_topology,
        warnings=warnings,
        style=label_style,
    )
    if display_label is not None:
        reaction_label = display_label.concise
        if display_label.status == "observed_edits":
            reaction_label_status = (
                "mapped_edit_summary"
                if "mapping" in edit_result.evidence
                else "observed_edit_summary"
            )
        elif display_label.status == "generic_pattern":
            reaction_label_status = (
                "mapped_generic_pattern"
                if "mapping" in edit_result.evidence
                else "generic_pattern"
            )
        elif display_label.status == "conflicting_evidence":
            reaction_label_status = "conflicting_edit_summary"
        elif display_label.status == "multi_event":
            reaction_label_status = "multi_event_edit_summary"
    return ReactionAnalysis(
        input_reaction_smiles=reaction_smiles,
        valid=True,
        reactants=parsed.reactants,
        agents=parsed.agents,
        products=parsed.products,
        candidates=tuple(candidates),
        selected_candidate=selected,
        selected_events=selected_events,
        transformation_class=(
            selected.transformation_class
            if selected
            else reaction_signature.transformation_class
            if reaction_signature is not None
            else None
        ),
        compatible_named_families=compatible_named_families,
        named_family=named_family,
        reaction_label=reaction_label,
        reaction_label_status=reaction_label_status,
        display_label=display_label,
        evidence_quality=effective_evidence,
        mapped_bond_changes=mapped_changes,
        spectator_groups=spectators,
        family_environment=family_environment,
        product_connection=product_connection,
        reaction_topology=reaction_topology,
        reaction_signature=reaction_signature,
        reaction_completeness=reaction_completeness,
        warnings=tuple(sorted(set(warnings))),
    )


__all__ = ["featurize_reaction"]
