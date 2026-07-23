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
    selected = analysis.selected_candidate
    if selected is None:
        return None, "QUERY_HAS_NO_SELECTED_REACTION_GRAMMAR"
    if selected.transformation_class != signature.transformation_class:
        return None, "RULE_QUERY_TRANSFORMATION_CONFLICT"
    topology = analysis.reaction_topology
    if topology is None:
        return None, "QUERY_HAS_NO_REACTION_TOPOLOGY"

    environment_by_role = {
        str(partner.role): partner
        for partner in (
            analysis.family_environment.partners
            if analysis.family_environment is not None
            else ()
        )
    }

    partners = []
    for role, reference in sorted(selected.role_assignments.items()):
        details = reference.details or {}
        h_count_value = details.get("h_count")
        h_count = int(h_count_value) if h_count_value is not None else None
        environment = environment_by_role.get(str(role))
        steric = dict(environment.steric) if environment is not None else {}
        electronic = (
            dict(environment.electronic) if environment is not None else {}
        )
        ortho_value = steric.get("ortho_substituent_count")
        alpha_branched_group_count = sum(
            bool(group.get("alpha_branched"))
            for group in steric.get("attached_groups") or ()
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
                steric_class=str(steric["class"])
                if steric.get("class")
                else None,
                electronic_class=str(electronic["class"])
                if electronic.get("class")
                else None,
                ortho_substituent_count=int(ortho_value)
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
            )
        )
    return (
        RuleQueryFacts(
            signature_id=signature.signature_id,
            reaction_signature_schema_version=signature.schema_version,
            transformation_class=str(signature.transformation_class or ""),
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
