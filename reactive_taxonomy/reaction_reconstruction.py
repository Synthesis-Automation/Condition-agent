"""Interpretation-independent reconstruction evidence from graph operators."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .connectivity_rewrite import apply_reaction_operator
from .reaction_archetypes import infer_bond_change_archetype
from .reaction_models import (
    ReactionComponent,
    ReactionReconstructionCandidate,
    ReactionSiteReference,
)
from .reaction_multi_events import (
    apply_rewrite_sequence,
    equivalent_multi_event_interpretations,
    exact_multi_event_reconstructions,
)
from .reaction_reconstruction_rules import load_reaction_reconstruction_rules


RawReconstructionCandidate = Tuple[
    Mapping[str, Any],
    Mapping[str, ReactionSiteReference],
]


@dataclass(frozen=True)
class ReactionReconstructionBuild:
    """Structural reconstruction proposals before interpretation."""

    raw_candidates: Tuple[RawReconstructionCandidate, ...]
    candidates: Tuple[ReactionReconstructionCandidate, ...]
    selected_candidate: ReactionReconstructionCandidate | None
    selected_events: Tuple[ReactionReconstructionCandidate, ...]
    evidence_quality: str
    product_contradicted_candidates: bool
    warnings: Tuple[str, ...]


def canonical_without_maps(smiles: str) -> str | None:
    """Return canonical isomeric SMILES after removing atom maps."""
    from rdkit import Chem

    molecule = parse_smiles(smiles)
    if molecule is None:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def _site_references(
    components: Iterable[ReactionComponent],
) -> Tuple[ReactionSiteReference, ...]:
    references = []
    for component in components:
        for site in component.compound_analysis.sites:
            references.append(
                ReactionSiteReference(
                    side="reactant",
                    component_index=component.component_index,
                    site_id=site.site_id,
                    site_type=site.site_type,
                    canonical_signature=site.canonical_signature,
                    chemist_label=site.chemist_label,
                    availability=site.availability,
                    atom_roles={
                        role: tuple(int(index) for index in indices)
                        for role, indices in (
                            site.details.get("atom_roles") or {}
                        ).items()
                    },
                    details=dict(site.details),
                )
            )
    return tuple(references)


def _matches(site: ReactionSiteReference, constraint: Mapping[str, Any]) -> bool:
    if site.site_type != constraint.get("site_type"):
        return False
    details = site.details
    scalar_any = (
        ("anchor_context_any", "anchor_context"),
        ("ring_context_any", "ring_context"),
        ("handle_token_any", "handle_token"),
        ("center_any", "center_token"),
        ("availability_any", None),
        ("carbonyl_subtype_any", "carbonyl_subtype"),
        ("activation_state_any", "activation_state"),
        ("donor_class_any", "donor_class"),
        ("source_kind_any", "source_kind"),
        ("departing_element_any", "departing_element"),
    )
    for constraint_name, detail_name in scalar_any:
        allowed = constraint.get(constraint_name)
        observed = site.availability if detail_name is None else details.get(detail_name)
        if allowed and observed not in allowed:
            return False
    if constraint.get("h_count_min") is not None and int(
        details.get("h_count", 0)
    ) < int(constraint["h_count_min"]):
        return False
    if constraint.get("attachment_h_count_min") is not None and int(
        details.get("attachment_h_count") or 0
    ) < int(constraint["attachment_h_count_min"]):
        return False
    if constraint.get("endpoint_h_count_max_min") is not None and max(
        (int(value) for value in details.get("endpoint_h_counts") or (0,)),
        default=0,
    ) < int(constraint["endpoint_h_count_max_min"]):
        return False
    for constraint_name, detail_name in (
        ("activation_relationship", "activation_relationship"),
        ("center_family", "center_family"),
    ):
        expected = constraint.get(constraint_name)
        if expected and details.get(detail_name) != expected:
            return False
    if constraint.get("contexts_any") and not set(
        details.get("contexts") or ()
    ).intersection(constraint["contexts_any"]):
        return False
    return True


def enumerate_reconstruction_assignments(
    components: Tuple[ReactionComponent, ...],
    rule: Mapping[str, Any],
) -> Tuple[Dict[str, ReactionSiteReference], ...]:
    """Enumerate site bindings for one structural reconstruction rule."""
    sites = _site_references(components)
    slot_names = tuple((rule.get("slots") or {}).keys())
    choices = tuple(
        tuple(site for site in sites if _matches(site, rule["slots"][slot]))
        for slot in slot_names
    )
    if any(not values for values in choices):
        return ()
    assignments = []
    for selected in product(*choices):
        assignment = dict(zip(slot_names, selected))
        if len({(site.component_index, site.site_id) for site in selected}) != len(
            selected
        ):
            continue
        valid = True
        for relationship in rule.get("slot_relationships") or ():
            indices = [
                assignment[slot].component_index
                for slot in relationship.get("slots") or ()
            ]
            relation = str(relationship.get("component_relation") or "")
            if relation == "same" and len(set(indices)) != 1:
                valid = False
            elif relation == "different" and len(set(indices)) != len(indices):
                valid = False
            if not valid:
                break
        if valid:
            assignments.append(assignment)
    return tuple(assignments)


def enumerate_reconstruction_candidates(
    components: Tuple[ReactionComponent, ...],
) -> Tuple[RawReconstructionCandidate, ...]:
    """Enumerate structural rule and site bindings."""
    return tuple(
        (rule, assignment)
        for rule in load_reaction_reconstruction_rules()
        for assignment in enumerate_reconstruction_assignments(components, rule)
    )


def _operator_assignment(
    rule: Mapping[str, Any],
    assignment: Mapping[str, ReactionSiteReference],
) -> Dict[str, ReactionSiteReference]:
    bindings = rule.get("operator_slot_bindings") or {}
    return {
        str(operator_slot): assignment[str(rule_slot)]
        for operator_slot, rule_slot in bindings.items()
    }


def _candidate(
    *,
    rule: Mapping[str, Any],
    assignment: Mapping[str, ReactionSiteReference],
    outcome: Any,
    observed_products: set[str],
    verification: str | None = None,
    predicted_product_smiles: str | None = None,
) -> ReactionReconstructionCandidate:
    predicted = predicted_product_smiles or (
        canonical_without_maps(outcome.predicted_product_smiles)
        if outcome.predicted_product_smiles
        else None
    )
    effective_verification = verification or (
        "exact_product_reconstruction"
        if predicted in observed_products
        else "product_mismatch"
        if predicted
        else "construction_failed"
    )
    return ReactionReconstructionCandidate(
        rule_id=str(rule["id"]),
        operator_id=str(rule["operator_id"]),
        rewrite_outcome_id=outcome.outcome_id,
        edit_archetype=infer_bond_change_archetype(
            outcome.predicted_bond_changes
        ),
        slot_assignments=dict(assignment),
        predicted_bond_changes=outcome.predicted_bond_changes,
        predicted_product_smiles=predicted,
        verification=effective_verification,
        predicted_stereo_changes=outcome.predicted_stereo_changes,
        warnings=outcome.warnings,
    )


def _equivalent_exact_candidates(
    candidates: Tuple[ReactionReconstructionCandidate, ...],
) -> bool:
    return len(
        {
            (
                candidate.rule_id,
                tuple(
                    sorted(
                        site.canonical_signature
                        for site in candidate.slot_assignments.values()
                    )
                ),
            )
            for candidate in candidates
        }
    ) == 1


def build_reaction_reconstruction_candidates(
    *,
    reactants: Tuple[ReactionComponent, ...],
    observed_products: set[str],
    invalid_supplied_mapping: bool,
    max_candidates: int,
) -> ReactionReconstructionBuild:
    """Collect exact-product operator evidence without interpretation semantics."""
    raw_values = list(enumerate_reconstruction_candidates(reactants))
    warnings: List[str] = []
    if len(raw_values) > max_candidates:
        raw_values = raw_values[:max_candidates]
        warnings.append("CANDIDATE_LIMIT_REACHED")
    raw = tuple(raw_values)
    candidates = []
    exact = []
    for rule, assignment in raw:
        operator_assignment = _operator_assignment(rule, assignment)
        output_roles = {
            operator_slot: str(rule_slot)
            for operator_slot, rule_slot in (
                rule.get("operator_slot_bindings") or {}
            ).items()
        }
        outcomes = apply_reaction_operator(
            str(rule["operator_id"]),
            operator_assignment,
            reactants,
            output_role_labels=output_roles,
        )
        for outcome in outcomes:
            candidate = _candidate(
                rule=rule,
                assignment=assignment,
                outcome=outcome,
                observed_products=observed_products,
            )
            candidates.append(candidate)
            if candidate.verification == "exact_product_reconstruction":
                exact.append(candidate)
            if len(candidates) >= max_candidates:
                warnings.append("CANDIDATE_LIMIT_REACHED")
                break
        if len(candidates) >= max_candidates:
            break
    selected = None
    selected_events: Tuple[ReactionReconstructionCandidate, ...] = ()
    evidence = "structural_reconstruction_candidates" if candidates else "unresolved"
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
            events = []
            for rule, assignment in chosen:
                operator_assignment = _operator_assignment(rule, assignment)
                output_roles = {
                    operator_slot: str(rule_slot)
                    for operator_slot, rule_slot in (
                        rule.get("operator_slot_bindings") or {}
                    ).items()
                }
                outcomes = apply_reaction_operator(
                    str(rule["operator_id"]),
                    operator_assignment,
                    reactants,
                    output_role_labels=output_roles,
                )
                if outcomes:
                    events.append(
                        _candidate(
                            rule=rule,
                            assignment=assignment,
                            outcome=outcomes[0],
                            observed_products=observed_products,
                            verification="exact_multi_event_reconstruction",
                            predicted_product_smiles=composite_product,
                        )
                    )
            selected_events = tuple(events)
            evidence = "exact_multi_event_reconstruction"
            if len(multi_exact) > 1:
                warnings.append("SYMMETRY_EQUIVALENT_MULTI_EVENT_ASSIGNMENTS")
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
            candidate.verification in {"construction_failed", "product_mismatch"}
            for candidate in candidate_tuple
        )
    )
    return ReactionReconstructionBuild(
        raw_candidates=raw,
        candidates=candidate_tuple,
        selected_candidate=selected,
        selected_events=selected_events,
        evidence_quality=evidence,
        product_contradicted_candidates=product_contradicted,
        warnings=tuple(sorted(set(warnings))),
    )


__all__ = [
    "RawReconstructionCandidate",
    "ReactionReconstructionBuild",
    "build_reaction_reconstruction_candidates",
    "canonical_without_maps",
    "enumerate_reconstruction_assignments",
    "enumerate_reconstruction_candidates",
]
