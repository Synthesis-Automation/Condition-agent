"""Stable categorical tokens for typed reactivity profiles."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Mapping, Tuple


def _mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if is_dataclass(value):
        return asdict(value)
    return {}


def reactivity_profile_tokens(profile: Any) -> Tuple[str, ...]:
    """Return identity-safe binned tokens without raw scores or atom indices."""
    value = _mapping(profile)
    if not value:
        return ()
    context = _mapping(value.get("context"))
    steric = _mapping(value.get("steric"))
    electronic = _mapping(value.get("electronic"))
    center = _mapping(value.get("reactive_center"))
    kind = str(value.get("context_kind") or context.get("context_kind") or "")
    tokens = []
    if kind:
        tokens.append(f"context:{kind}")
    access = steric.get("accessibility_class")
    if access:
        tokens.append(f"steric_access:{access}")
    axis = electronic.get("activation_axis")
    electronic_class = electronic.get("activation_class")
    if axis and electronic_class:
        tokens.append(f"electronic:{axis}:{electronic_class}")
    if kind == "aromatic":
        family = context.get("ring_family")
        if family:
            tokens.append(f"aromatic_family:{family}")
        for size in sorted(set(context.get("ring_sizes") or ())):
            tokens.append(f"aromatic_ring_size:{int(size)}")
        tokens.append(
            "ortho_occupancy:"
            f"{int(context.get('ortho_occupancy_count') or 0)}_of_"
            f"{int(context.get('ortho_capacity') or 0)}"
        )
        burden = context.get("ortho_burden_class")
        if burden:
            tokens.append(f"ortho_burden:{burden}")
        for item in context.get("heteroatoms") or ():
            record = _mapping(item)
            if record.get("positional_relation") == "anchor":
                continue
            tokens.append(
                "aromatic_heteroatom:"
                f"{record.get('element')}@{record.get('positional_relation')}:"
                f"{record.get('aromatic_role')}"
            )
    elif kind == "alkyl":
        tokens.extend(
            (
                f"alkyl_center:{context.get('carbon_substitution')}",
                "alpha_branching:"
                + ("branched" if context.get("alpha_branched") else "unbranched"),
                "beta_hydrogen:"
                + (
                    "present"
                    if int(context.get("beta_hydrogen_count") or 0) > 0
                    else "absent"
                ),
            )
        )
        for key in ("benzylic", "allylic", "propargylic", "cyclic"):
            if context.get(key):
                tokens.append(f"alkyl_activation:{key}")
    elif kind in {"alkenyl", "alkynyl"}:
        endpoints = context.get("endpoint_substitution") or ()
        if endpoints:
            tokens.append(
                "pi_substitution:" + "_".join(str(value) for value in endpoints)
            )
        conjugation = context.get("conjugation_class")
        if conjugation:
            tokens.append(f"pi_conjugation:{conjugation}")
        if context.get("stereochemistry"):
            tokens.append(f"pi_stereo:{context.get('stereochemistry')}")
        if kind == "alkynyl":
            tokens.append(
                "alkyne_terminus:"
                + ("terminal" if context.get("terminal") else "internal")
            )
    elif kind in {"acyl", "sulfonyl", "phosphoryl"}:
        tokens.append(f"activated_center:{context.get('center_class')}")
        tokens.append(f"conjugation:{context.get('conjugation_class')}")
    elif kind == "heteroatom":
        tokens.extend(
            (
                f"heteroatom:{context.get('element')}",
                f"heteroatom_substitution:{context.get('substitution_class')}",
                f"lone_pair:{context.get('lone_pair_class')}",
                "attached_alpha_branching:"
                f"{int(context.get('alpha_branched_group_count') or 0)}",
            )
        )
        availability = center.get("lone_pair_availability")
        if availability:
            tokens.append(f"lone_pair_availability:{availability}")
        for group in context.get("attached_groups") or ():
            record = _mapping(group)
            if record.get("attachment_carbon_class"):
                tokens.append(
                    "attached_group:"
                    f"{record.get('context')}:"
                    f"{record.get('attachment_carbon_class')}"
                )
    for modifier in value.get("modifiers") or ():
        record = _mapping(modifier)
        modifier_type = record.get("modifier_type")
        class_name = record.get("class_name")
        if modifier_type and class_name:
            tokens.append(f"modifier:{modifier_type}:{class_name}")
    return tuple(sorted(set(str(token) for token in tokens if token)))


__all__ = ["reactivity_profile_tokens"]
