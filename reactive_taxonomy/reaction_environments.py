"""Reaction-family overlays derived from selected sites and raw environments."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from .reaction_models import (
    ReactionCandidate,
    ReactionComponent,
    ReactionFamilyEnvironment,
    ReactionPartnerEnvironment,
    ReactionSpectatorGroup,
)

_RULES_PATH = Path(__file__).with_name("definitions") / "descriptor_rules.v1.json"


@lru_cache(maxsize=1)
def _family_rules() -> Dict[str, Dict[str, Any]]:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle).get("reaction_families") or {})


def _component(components: Tuple[ReactionComponent, ...], index: int) -> ReactionComponent | None:
    return next((component for component in components if component.component_index == index), None)


def _unique(values: Iterable[str]) -> Tuple[str, ...]:
    return tuple(sorted({str(value) for value in values if value}))


def build_reaction_family_environment(
    components: Tuple[ReactionComponent, ...],
    selected: ReactionCandidate | None,
    spectators: Tuple[ReactionSpectatorGroup, ...],
    evidence_quality: str,
) -> ReactionFamilyEnvironment | None:
    """Build an interpretable family feature overlay for a selected event."""
    if selected is None or "suzuki_miyaura" not in selected.compatible_named_families:
        return None
    family_id = "suzuki_miyaura"
    rules = _family_rules().get(family_id) or {}
    coordination_tags = set(rules.get("coordination_tags") or [])
    xh_tags = set(rules.get("unprotected_xh_tags") or [])
    competing_types = set(rules.get("competing_site_types") or [])
    partners: List[ReactionPartnerEnvironment] = []
    family_flags: List[str] = []
    for role in rules.get("roles") or selected.role_assignments:
        reference = selected.role_assignments.get(str(role))
        if reference is None:
            continue
        component = _component(components, reference.component_index)
        if component is None:
            continue
        environment = next((
            item for item in component.compound_analysis.site_environments
            if item.site_id == reference.site_id
        ), None)
        role_spectators = [
            group for group in spectators
            if group.component_index == reference.component_index
        ]
        competing = [
            site.chemist_label for site in component.compound_analysis.sites
            if site.site_id != reference.site_id and site.site_type in competing_types
        ]
        coordination = _unique(
            group.group_id for group in role_spectators
            if coordination_tags.intersection(group.tags)
        )
        unprotected_xh = _unique(
            group.group_id for group in role_spectators
            if xh_tags.intersection(group.tags)
        )
        flags: List[str] = []
        anchor_context = reference.details.get("anchor_context")
        handle_token = reference.details.get("handle_token")
        steric = dict(environment.steric) if environment else {}
        electronic = dict(environment.electronic) if environment else {}
        if anchor_context == "HeteroAr":
            flags.append("heteroaryl_partner")
        if steric.get("class") in {"ortho_hindered", "tertiary"}:
            flags.append("sterically_hindered")
        if electronic.get("class") in {"electron_poor", "electron_rich"}:
            flags.append(str(electronic["class"]))
        if competing:
            flags.append("competing_reactive_handle")
        if coordination:
            flags.append("coordination_risk")
        if unprotected_xh:
            flags.append("unprotected_xh")
        if role == "electrophile" and handle_token in set(rules.get("challenging_leaving_groups") or []):
            flags.append("challenging_oxidative_addition")
        if (
            role == "transfer_partner"
            and anchor_context in set(rules.get("protodeboronation_contexts") or [])
            and handle_token in set(rules.get("protodeboronation_handles") or [])
        ):
            flags.append("protodeboronation_risk")
        family_flags.extend(f"{role}:{flag}" for flag in flags)
        partners.append(ReactionPartnerEnvironment(
            role=str(role),
            component_index=reference.component_index,
            site_id=reference.site_id,
            handle_token=str(handle_token) if handle_token else None,
            anchor_context=str(anchor_context) if anchor_context else None,
            chemist_label=reference.chemist_label,
            steric=steric,
            electronic=electronic,
            nearby_groups=environment.nearby_groups if environment else (),
            spectator_group_ids=_unique(group.group_id for group in role_spectators),
            competing_site_labels=_unique(competing),
            coordination_group_ids=coordination,
            unprotected_xh_group_ids=unprotected_xh,
            flags=_unique(flags),
        ))
    return ReactionFamilyEnvironment(
        family_id=family_id,
        partners=tuple(partners),
        flags=_unique(family_flags),
        evidence=evidence_quality,
    )


__all__ = ["build_reaction_family_environment"]
