"""Shared chemist-facing review summary for reaction analyses and records.

The summary keeps the structural display label and minimized reaction-core
label side by side.  Neither label is promoted to observed chemistry when its
underlying evidence is unavailable, and neither is used for signature identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Mapping, Optional, Tuple

from .descriptors import render_reactivity_profile
from .reaction_models import ReactionAnalysis


REACTION_REVIEW_SUMMARY_SCHEMA_VERSION = "2.0"
_SUBSCRIPT_TRANSLATION = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


@dataclass(frozen=True)
class ReactionReviewSummary:
    """Priority chemistry fields shared by GUI and dataset review views."""

    detailed_reaction_label: str
    reaction_core_equation: str
    core_limiter: str
    atom_level_core_equation: str
    spectators: str
    electronic_steric_analysis: str
    display_status: str
    display_evidence: str
    core_evidence_status: str
    schema_version: str = REACTION_REVIEW_SUMMARY_SCHEMA_VERSION


def _member(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def _readable_token(value: Any) -> str:
    raw = value.value if isinstance(value, Enum) else value
    return str(raw or "").replace("_", " ").strip()


def _formula_text(value: Any) -> str:
    return str(value or "").translate(_SUBSCRIPT_TRANSLATION)


def _spectator_summary(groups: Any) -> str:
    grouped: dict[tuple[str, str], list[Optional[int]]] = {}
    for group in groups or ():
        group_id = str(_member(group, "group_id", "") or "")
        label = str(_member(group, "chemist_label", "") or "")
        if not group_id and not label:
            continue
        distance = _member(group, "graph_distance")
        grouped.setdefault((label, group_id), []).append(
            int(distance) if distance is not None else None
        )

    values = []
    for (label, group_id), occurrences in sorted(grouped.items()):
        readable_group = _readable_token(group_id)
        display = _formula_text(label) or readable_group
        if readable_group and readable_group.casefold() != display.casefold():
            display += f" [{readable_group}]"
        if len(occurrences) > 1:
            display = f"{len(occurrences)}× {display}"
        distances = sorted(
            {distance for distance in occurrences if distance is not None}
        )
        if distances:
            display += f" (d={'/'.join(str(value) for value in distances)})"
        values.append(display)
    return "; ".join(values)


def _partner_label(partner: Any) -> str:
    role = _readable_token(_member(partner, "role"))
    if role:
        return role
    component = int(_member(partner, "component_index", 0) or 0) + 1
    chemist_label = _formula_text(_member(partner, "chemist_label"))
    return f"P{component} ({chemist_label})" if chemist_label else f"P{component}"


def _partner_sort_key(partner: Any) -> Tuple[int, str, str]:
    return (
        int(_member(partner, "component_index", 0) or 0),
        str(_member(partner, "role", "") or ""),
        str(_member(partner, "partner_id", "") or ""),
    )


def _partner_environment_summary(partners: Any) -> str:
    values = []
    for partner in sorted(tuple(partners or ()), key=_partner_sort_key):
        profile = _member(partner, "reactivity_profile")
        if profile is None:
            continue
        values.append(
            f"{_partner_label(partner)}: {render_reactivity_profile(profile)}"
        )
    return " | ".join(values)


def build_reaction_review_summary(
    source: ReactionAnalysis | Mapping[str, Any],
) -> ReactionReviewSummary:
    """Build one review summary from live analysis or its serialized record.

    Live analyses and recommendation records both use the same nested
    ``reaction_label`` contract.
    """
    display = _member(source, "reaction_label")
    core = _member(source, "reaction_core")
    abstraction = _member(core, "abstraction")
    signature = _member(source, "reaction_signature")

    detailed_label = str(_member(display, "detailed", "") or "")

    spectators = _member(signature, "spectator_groups")
    if spectators is None:
        spectators = _member(source, "spectator_groups", ())

    interpretation = _member(source, "interpretation")
    family_environment = _member(interpretation, "family_environment")
    if family_environment is None:
        family_environment = _member(source, "family_environment")
    partners = _member(interpretation, "partners")
    if not partners:
        partners = _member(family_environment, "partners")
    if not partners:
        partners = _member(signature, "partners", ())

    return ReactionReviewSummary(
        detailed_reaction_label=detailed_label,
        reaction_core_equation=str(
            _member(abstraction, "general_label", "")
            or _member(core, "generic_label", "")
            or ""
        ),
        core_limiter=str(
            _member(abstraction, "limiter_label", "") or ""
        ),
        atom_level_core_equation=(
            str(_member(core, "generic_label", "") or "")
            if abstraction is not None
            else ""
        ),
        spectators=_spectator_summary(spectators),
        electronic_steric_analysis=_partner_environment_summary(partners),
        display_status=str(_member(display, "status", "") or ""),
        display_evidence=str(_member(display, "evidence", "") or ""),
        core_evidence_status=str(
            _member(core, "evidence_status", "") or ""
        ),
    )


def format_reaction_review_summary(summary: ReactionReviewSummary) -> str:
    """Render the four evidence-preserving priority lines for a GUI review."""
    detailed = summary.detailed_reaction_label or "Unavailable"
    graphic = summary.reaction_core_equation or "Unavailable"
    spectators = summary.spectators or "None detected"
    profile = summary.electronic_steric_analysis or "Unavailable"
    lines = [
        f"Detailed reaction label: {detailed}",
        f"Graphic core label: {graphic}",
    ]
    if summary.core_limiter:
        lines.append(f"Graphic core limiter: {summary.core_limiter}")
    if summary.atom_level_core_equation:
        lines.append(f"Atom-level core: {summary.atom_level_core_equation}")
    lines.extend(
        (
            f"Spectators: {spectators}",
            f"Electronic / steric analysis: {profile}",
        )
    )
    return "\n".join(lines)


__all__ = [
    "REACTION_REVIEW_SUMMARY_SCHEMA_VERSION",
    "ReactionReviewSummary",
    "build_reaction_review_summary",
    "format_reaction_review_summary",
]
