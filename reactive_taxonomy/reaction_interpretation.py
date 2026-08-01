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
from .reaction_models import (
    ReactionCandidate,
    ReactionComponent,
    ReactionSiteReference,
)
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
        reaction_label=render_reaction_label(
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
                        reaction_label=render_reaction_label(
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


__all__ = [
    "RawReactionCandidate",
    "ReactionInterpretationBuild",
    "build_reaction_interpretation_candidates",
    "canonical_without_maps",
]
