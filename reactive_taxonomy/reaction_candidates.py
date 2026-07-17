"""Grammar-driven reaction-site role assignment."""

from __future__ import annotations

from itertools import product
from typing import Any, Dict, Iterable, List, Tuple

from .reaction_grammars import load_reaction_grammars
from .reaction_models import ReactionComponent, ReactionSiteReference


def _site_references(
    components: Iterable[ReactionComponent],
) -> List[ReactionSiteReference]:
    references: List[ReactionSiteReference] = []
    for component in components:
        for site in component.compound_analysis.sites:
            raw_roles = site.details.get("atom_roles") or {}
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
                        for role, indices in raw_roles.items()
                    },
                    details=dict(site.details),
                )
            )
    return references


def _matches(site: ReactionSiteReference, constraint: Dict[str, Any]) -> bool:
    if site.site_type != constraint.get("site_type"):
        return False
    details = site.details
    if (
        constraint.get("anchor_context_any")
        and details.get("anchor_context") not in constraint["anchor_context_any"]
    ):
        return False
    if (
        constraint.get("handle_token_any")
        and details.get("handle_token") not in constraint["handle_token_any"]
    ):
        return False
    if (
        constraint.get("center_any")
        and details.get("center_token") not in constraint["center_any"]
    ):
        return False
    if constraint.get("h_count_min") is not None and int(
        details.get("h_count", 0)
    ) < int(constraint["h_count_min"]):
        return False
    if (
        constraint.get("availability_any")
        and site.availability not in constraint["availability_any"]
    ):
        return False
    if (
        constraint.get("center_family")
        and details.get("center_family") != constraint["center_family"]
    ):
        return False
    if (
        constraint.get("carbonyl_subtype_any")
        and details.get("carbonyl_subtype") not in constraint["carbonyl_subtype_any"]
    ):
        return False
    if (
        constraint.get("activation_state_any")
        and details.get("activation_state") not in constraint["activation_state_any"]
    ):
        return False
    if constraint.get("endpoint_h_count_max_min") is not None and max(
        (int(value) for value in details.get("endpoint_h_counts") or (0,)),
        default=0,
    ) < int(constraint["endpoint_h_count_max_min"]):
        return False
    return True


def enumerate_role_assignments(
    components: Tuple[ReactionComponent, ...], grammar: Dict[str, Any]
) -> List[Dict[str, ReactionSiteReference]]:
    sites = _site_references(components)
    role_names = list((grammar.get("roles") or {}).keys())
    choices = [
        [site for site in sites if _matches(site, grammar["roles"][role])]
        for role in role_names
    ]
    if any(not values for values in choices):
        return []
    assignments: List[Dict[str, ReactionSiteReference]] = []
    for selected in product(*choices):
        assignment = dict(zip(role_names, selected))
        if len({(site.component_index, site.site_id) for site in selected}) != len(
            selected
        ):
            continue
        valid = True
        for left_role, right_role in grammar.get("distinct_components") or []:
            if (
                assignment[left_role].component_index
                == assignment[right_role].component_index
            ):
                valid = False
                break
        if valid:
            assignments.append(assignment)
    return assignments


def enumerate_reaction_candidates(
    components: Tuple[ReactionComponent, ...],
) -> List[Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]]:
    return [
        (grammar, assignment)
        for grammar in load_reaction_grammars()
        for assignment in enumerate_role_assignments(components, grammar)
    ]


__all__ = ["enumerate_reaction_candidates", "enumerate_role_assignments"]
