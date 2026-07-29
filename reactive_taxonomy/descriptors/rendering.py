"""Shared compact and expanded rendering for reactivity profiles."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Mapping

from .registry import rendering_rules


def _mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if is_dataclass(value):
        return asdict(value)
    return {}


def _readable(value: Any) -> str:
    return str(value or "").replace("_", " ").strip()


def _context_identity(context: Mapping[str, Any]) -> str:
    kind = str(context.get("context_kind") or "other")
    if kind == "aromatic":
        family = _readable(context.get("ring_family")) or "aromatic system"
        sizes = tuple(context.get("ring_sizes") or ())
        size_text = "/".join(str(value) for value in sizes)
        identity = f"{family} ({size_text}-membered)" if size_text else family
        hetero = []
        for item in context.get("heteroatoms") or ():
            record = _mapping(item)
            if record.get("positional_relation") == "anchor":
                continue
            hetero.append(
                f"{record.get('element')}@"
                f"{_readable(record.get('positional_relation'))}"
            )
        if hetero:
            identity += f", heteroatoms {','.join(hetero)}"
        if context.get("fused"):
            identity += ", fused"
        return identity
    if kind == "alkyl":
        identity = f"{_readable(context.get('carbon_substitution'))} alkyl"
        activations = [
            label
            for key, label in (
                ("benzylic", "benzylic"),
                ("allylic", "allylic"),
                ("propargylic", "propargylic"),
            )
            if context.get(key)
        ]
        if activations:
            identity += f", {'/'.join(activations)}"
        ring_sizes = tuple(context.get("ring_sizes") or ())
        if ring_sizes:
            identity += f", cyclic ({min(ring_sizes)}-membered)"
        return identity
    if kind == "alkenyl":
        identity = _readable(context.get("alkene_class")) or "alkenyl"
        stereo = _readable(context.get("stereochemistry"))
        if stereo:
            identity += f", {stereo}"
        if context.get("cyclic"):
            identity += f", cyclic ({context.get('ring_size')}-membered)"
        conjugation = _readable(context.get("conjugation_class"))
        if conjugation and conjugation != "isolated":
            identity += f", {conjugation}"
        return identity
    if kind == "alkynyl":
        identity = "terminal alkynyl" if context.get("terminal") else "internal alkynyl"
        conjugation = _readable(context.get("conjugation_class"))
        if conjugation and conjugation != "isolated":
            identity += f", {conjugation}"
        return identity
    if kind in {"acyl", "sulfonyl", "phosphoryl"}:
        identity = _readable(context.get("center_class")) or _readable(kind)
        conjugation = _readable(context.get("conjugation_class"))
        if conjugation and conjugation != "not aryl conjugated":
            identity += f", {conjugation}"
        return identity
    if kind == "heteroatom":
        element = str(context.get("element") or "heteroatom")
        substitution = _readable(context.get("substitution_class"))
        resonance = _readable(context.get("resonance_class"))
        attached = []
        for item in context.get("attached_groups") or ():
            record = _mapping(item)
            attachment_class = _readable(
                record.get("attachment_carbon_class")
            )
            if attachment_class:
                attached.append(
                    f"{record.get('context')}: {attachment_class}"
                )
        identity = ", ".join(
            value
            for value in (f"{substitution} {element}", resonance)
            if value
        )
        if attached:
            identity += f", attached {', '.join(sorted(set(attached)))}"
        return identity
    return _readable(context.get("reason")) or "unclassified context"


def _steric_text(
    context: Mapping[str, Any],
    steric: Mapping[str, Any],
) -> str:
    kind = str(context.get("context_kind") or "")
    access = _readable(steric.get("accessibility_class")) or "unresolved"
    if kind == "aromatic":
        occupied = context.get("ortho_occupancy_count")
        capacity = context.get("ortho_capacity")
        burden = _readable(context.get("ortho_burden_class"))
        if occupied is not None and capacity:
            return f"ortho burden {burden} ({occupied}/{capacity})"
    if kind == "alkyl":
        details = []
        if context.get("alpha_branched"):
            details.append("alpha-branched")
        beta = int(context.get("beta_branch_count") or 0)
        if beta:
            details.append(f"beta branches {beta}")
        return f"access {access}" + (f" ({', '.join(details)})" if details else "")
    return f"access {access}"


def _electronic_text(electronic: Mapping[str, Any]) -> str:
    axis = str(electronic.get("activation_axis") or "")
    label = str(
        (rendering_rules().get("labels") or {}).get(axis)
        or _readable(axis)
        or "electronic influence"
    )
    value = _readable(electronic.get("activation_class")) or "unresolved"
    return f"{label} {value}"


def _center_text(
    context: Mapping[str, Any],
    center: Mapping[str, Any],
) -> str:
    values = []
    if center.get("element") in {"N", "O", "S", "P"}:
        substitution = _readable(center.get("substitution_class"))
        hydrogens = int(center.get("hydrogen_count") or 0)
        if substitution:
            values.append(
                f"{substitution} {center.get('element')}"
                + (f"-H{hydrogens}" if hydrogens > 1 else "-H" if hydrogens else "")
            )
        lone_pair = _readable(center.get("lone_pair_class"))
        availability = _readable(center.get("lone_pair_availability"))
        if lone_pair:
            values.append(
                f"lone pair {lone_pair}"
                + (f"/{availability}" if availability else "")
            )
    if context.get("context_kind") == "alkyl":
        beta_h = context.get("beta_hydrogen_count")
        if beta_h is not None:
            values.append(
                f"beta-H {'present' if int(beta_h) > 0 else 'absent'}"
            )
    return "; ".join(values)


def _liability_text(modifiers: Any) -> str:
    values = []
    for item in modifiers or ():
        modifier = _mapping(item)
        modifier_type = str(modifier.get("modifier_type") or "")
        class_name = _readable(modifier.get("class_name"))
        if modifier_type in {"strain", "coordination", "redox"}:
            values.append(class_name)
    return ", ".join(sorted(set(value for value in values if value)))


def render_reactivity_profile(profile: Any) -> str:
    """Render one concise, fixed-order chemist-facing profile."""
    value = _mapping(profile)
    if not value:
        return "reactivity profile unavailable"
    context = _mapping(value.get("context"))
    steric = _mapping(value.get("steric"))
    electronic = _mapping(value.get("electronic"))
    center = _mapping(value.get("reactive_center"))
    clauses = [
        _context_identity(context),
        _steric_text(context, steric),
        _electronic_text(electronic),
        _center_text(context, center),
        _liability_text(value.get("modifiers")),
    ]
    return "; ".join(clause for clause in clauses if clause)


def render_reactivity_profile_expanded(profile: Any) -> str:
    """Render compact identity plus scores, methods, and contribution evidence."""
    value = _mapping(profile)
    if not value:
        return "reactivity profile unavailable"
    steric = _mapping(value.get("steric"))
    electronic = _mapping(value.get("electronic"))
    steric_evidence = _mapping(steric.get("evidence"))
    electronic_evidence = _mapping(electronic.get("evidence"))
    electronic_contributors = ", ".join(
        f"{_mapping(item).get('source_id')}:"
        f"{float(_mapping(item).get('contribution') or 0.0):+.2f}"
        for item in electronic.get("contributions") or ()
    ) or "none recognized"
    return (
        f"{render_reactivity_profile(value)} | "
        f"steric score={float(steric.get('accessibility_score') or 0.0):.2f} "
        f"[{steric_evidence.get('method') or 'unknown'}]; "
        f"electronic score={float(electronic.get('activation_score') or 0.0):+.2f} "
        f"[{electronic_evidence.get('method') or 'unknown'}]; "
        f"contributors={electronic_contributors}"
    )


__all__ = [
    "render_reactivity_profile",
    "render_reactivity_profile_expanded",
]
