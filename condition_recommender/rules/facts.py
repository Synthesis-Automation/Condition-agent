"""Project verified reactive-taxonomy observations into rule facts."""

from __future__ import annotations

from typing import Any, Optional

from .models import PartnerRuleFacts, RuleQueryFacts


def build_rule_query_facts(
    analysis: Any,
) -> tuple[Optional[RuleQueryFacts], Optional[str]]:
    """Build facts only for one verified, grammar-assigned reaction event."""
    if not analysis.valid:
        return None, analysis.error or "INVALID_REACTION"
    signature = analysis.reaction_signature
    if signature is None:
        return None, "QUERY_HAS_NO_USABLE_REACTION_SIGNATURE"
    if signature.event_scope != "single_event":
        return None, "RULE_QUERY_REQUIRES_SINGLE_EVENT"
    interpretation = analysis.interpretation
    selected = (
        interpretation.selected_candidate
        if interpretation is not None
        else None
    )
    if selected is None:
        return None, "QUERY_HAS_NO_SELECTED_REACTION_GRAMMAR"
    if (
        interpretation.evidence_quality
        == "conflicting_interpretation_evidence"
    ):
        return None, "RULE_QUERY_INTERPRETATION_CONFLICT"
    topology = analysis.reaction_topology
    if topology is None:
        return None, "QUERY_HAS_NO_REACTION_TOPOLOGY"

    environment_by_role = {
        str(partner.role): partner
        for partner in interpretation.partners
        if partner.role is not None
    }

    partners = []
    for role, reference in sorted(selected.role_assignments.items()):
        details = reference.details or {}
        h_count_value = details.get("h_count")
        h_count = int(h_count_value) if h_count_value is not None else None
        environment = environment_by_role.get(str(role))
        profile = (
            environment.reactivity_profile
            if environment is not None
            else None
        )
        context = profile.context if profile is not None else None
        context_kind = profile.context_kind if profile is not None else None
        ortho_value = (
            context.ortho_occupancy_count
            if context_kind == "aromatic"
            else None
        )
        alpha_branched_group_count = (
            int(context.alpha_branched_group_count)
            if context_kind == "heteroatom"
            else int(bool(context.alpha_branched))
            if context_kind == "alkyl"
            else 0
        )
        partners.append(
            PartnerRuleFacts(
                role=str(role),
                component_index=reference.component_index,
                site_id=reference.site_id,
                site_type=reference.site_type,
                availability=reference.availability,
                anchor_context=str(details["anchor_context"])
                if details.get("anchor_context")
                else None,
                handle_token=str(details["handle_token"])
                if details.get("handle_token")
                else None,
                center_token=str(details["center_token"])
                if details.get("center_token")
                else None,
                derived_family=str(details["derived_family"])
                if details.get("derived_family")
                else None,
                h_count=h_count,
                retained_contexts=tuple(
                    sorted(str(value) for value in details.get("contexts") or ())
                ),
                ortho_occupancy=int(ortho_value)
                if ortho_value is not None
                else None,
                alpha_branched_group_count=alpha_branched_group_count,
                environment_flags=tuple(
                    sorted(
                        str(value)
                        for value in (
                            environment.flags
                            if environment is not None
                            else ()
                        )
                    )
                ),
                context_kind=context_kind,
                ring_family=(
                    str(context.ring_family)
                    if context_kind == "aromatic"
                    else None
                ),
                steric_accessibility=(
                    profile.steric.accessibility_class
                    if profile is not None
                    else None
                ),
                ortho_burden_class=(
                    str(context.ortho_burden_class)
                    if context_kind == "aromatic"
                    else None
                ),
                electronic_axis=(
                    profile.electronic.activation_axis
                    if profile is not None
                    else None
                ),
                alkyl_substitution=(
                    str(context.carbon_substitution)
                    if context_kind == "alkyl"
                    else None
                ),
                beta_hydrogen_status=(
                    "present"
                    if context_kind == "alkyl"
                    and int(context.beta_hydrogen_count) > 0
                    else "absent"
                    if context_kind == "alkyl"
                    else None
                ),
                lone_pair_availability=(
                    profile.reactive_center.lone_pair_availability
                    if profile is not None
                    else None
                ),
                reactivity_modifiers=tuple(
                    sorted(
                        f"{modifier.modifier_type}:{modifier.class_name}"
                        for modifier in (
                            profile.modifiers if profile is not None else ()
                        )
                    )
                ),
            )
        )
    return (
        RuleQueryFacts(
            signature_id=signature.signature_id,
            reaction_signature_schema_version=signature.schema_version,
            transformation_class=str(selected.transformation_class or ""),
            event_scope=signature.event_scope,
            evidence_quality=analysis.evidence_quality,
            reaction_scope=str(topology.reaction_scope),
            partners=tuple(partners),
            taxonomy_definition_versions=tuple(
                sorted(
                    (str(key), str(value))
                    for key, value in signature.definition_versions.items()
                )
            ),
        ),
        None,
    )


__all__ = ["build_rule_query_facts"]
