"""Deterministic observation-first reaction display-label rendering."""

from __future__ import annotations

import json
import re
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Tuple

from .reaction_contextual_labels import ContextualTransformationLabel
from .reaction_label_patterns import match_reaction_label_pattern
from .reaction_labels import load_fragment_context_symbols
from .reaction_models import (
    ReactionDisplayLabel,
    ReactionEdit,
    ReactionEvent,
    ReactionLabelClause,
    ReactionTopology,
)
from .reaction_topology import topology_label_prefix

_PATH = Path(__file__).with_name("definitions") / "reaction_label_rendering.v1.json"


@lru_cache(maxsize=1)
def load_reaction_label_rendering() -> dict[str, Any]:
    """Load the versioned declarative edit-label rendering rules."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != "1.4":
        raise ValueError("Unsupported reaction-label rendering schema")
    return dict(payload)


def _style(name: str) -> dict[str, str]:
    styles = load_reaction_label_rendering()["styles"]
    if name not in styles:
        raise ValueError(f"Unknown reaction-label style: {name}")
    return dict(styles[name])


def _bond(order: Optional[str], style: dict[str, str]) -> str:
    return {
        "SINGLE": style["single"],
        "DOUBLE": style["double"],
        "TRIPLE": style["triple"],
        "AROMATIC": style["aromatic"],
    }.get(str(order or "SINGLE").upper(), style["single"])


def _ordered_atoms(edit: ReactionEdit) -> Tuple[Any, ...]:
    atoms = tuple(atom for atom in (edit.atom_1, edit.atom_2) if atom is not None)
    precedence = {
        element: index
        for index, element in enumerate(
            load_reaction_label_rendering()["element_precedence"]
        )
    }
    return tuple(
        sorted(
            atoms,
            key=lambda atom: (
                precedence.get(atom.element, 999),
                atom.element,
                atom.atom_map_number if atom.atom_map_number is not None else -1,
                atom.component_index,
                atom.atom_index,
            ),
        )
    )


def _mapped_atom(atom: Any) -> str:
    template = load_reaction_label_rendering()["templates"]["mapped_atom"]
    if atom.atom_map_number is None:
        return str(atom.element)
    return str(template).format(
        element=atom.element, map_number=atom.atom_map_number
    )


def render_reaction_label_clause(
    edit: ReactionEdit,
    *,
    style: str = "unicode",
) -> ReactionLabelClause:
    """Render one normalized edit without inferring a reaction family."""
    styling = _style(style)
    templates = load_reaction_label_rendering()["templates"]
    atoms = _ordered_atoms(edit)
    elements = tuple(atom.element for atom in atoms)
    maps = tuple(
        int(atom.atom_map_number)
        for atom in atoms
        if atom.atom_map_number is not None
    )
    if edit.edit_type == "hydrogen_change":
        center = edit.atom_1.element
        detailed_center = _mapped_atom(edit.atom_1)
        gained = bool(edit.new_order and not edit.old_order)
        key = "hydrogen_gain" if gained else "hydrogen_loss"
        concise = str(templates[key]).format(
            center=center, bond=styling["single"]
        )
        detailed = str(templates[key]).format(
            center=detailed_center, bond=styling["single"]
        )
        elements = (edit.atom_1.element, "H")
    else:
        left, right = atoms
        concise_values = {"left": left.element, "right": right.element}
        detailed_values = {"left": _mapped_atom(left), "right": _mapped_atom(right)}
        if edit.edit_type == "formed":
            values = {"bond": _bond(edit.new_order, styling)}
        elif edit.edit_type == "broken":
            values = {"bond": _bond(edit.old_order, styling)}
        else:
            values = {
                "old_bond": _bond(edit.old_order, styling),
                "new_bond": _bond(edit.new_order, styling),
                "arrow": styling["arrow"],
            }
        concise = str(templates[edit.edit_type]).format(**concise_values, **values)
        detailed = str(templates[edit.edit_type]).format(**detailed_values, **values)
    return ReactionLabelClause(
        edit_type=edit.edit_type,
        concise=concise,
        detailed=detailed,
        elements=elements,
        atom_map_numbers=maps,
        old_order=edit.old_order,
        new_order=edit.new_order,
        evidence=edit.evidence,
        confidence=edit.confidence,
    )


def _ordered_clauses(clauses: Iterable[ReactionLabelClause]) -> tuple[ReactionLabelClause, ...]:
    order = {
        edit_type: index
        for index, edit_type in enumerate(
            load_reaction_label_rendering()["clause_order"]
        )
    }
    return tuple(
        sorted(
            clauses,
            key=lambda clause: (
                order.get(clause.edit_type, 999),
                clause.concise,
                clause.detailed,
            ),
        )
    )


def _compose_concise(
    clauses: Sequence[ReactionLabelClause], *, style: str
) -> str:
    rendering = load_reaction_label_rendering()
    styling = _style(style)
    template = rendering["templates"]["counted_clause"]
    counts = Counter(clause.concise for clause in clauses)
    first_position = {
        clause.concise: index for index, clause in reversed(list(enumerate(clauses)))
    }
    parts = []
    for concise in sorted(counts, key=first_position.get):
        count = counts[concise]
        parts.append(
            str(template).format(
                count=count, times=styling["times"], clause=concise
            )
            if count > 1
            else concise
        )
    return styling["separator"].join(parts)


def _edit_transition_label(
    edits: Sequence[ReactionEdit], *, style: str
) -> str:
    """Render a mechanism-neutral before-to-after view from observed edits."""
    styling = _style(style)
    before = []
    after = []
    endpoint_elements = []
    for edit in edits:
        atoms = _ordered_atoms(edit)
        endpoint_elements.extend(atom.element for atom in atoms)
        if edit.edit_type == "hydrogen_change":
            token = (
                f"{edit.atom_1.element}"
                f"{_bond('SINGLE', styling)}H"
            )
            if edit.old_order and not edit.new_order:
                before.append(token)
            elif edit.new_order and not edit.old_order:
                after.append(token)
            continue
        if len(atoms) != 2:
            continue
        left, right = atoms
        if edit.edit_type in {"broken", "order_changed"}:
            before.append(
                f"{left.element}{_bond(edit.old_order, styling)}{right.element}"
            )
        if edit.edit_type in {"formed", "order_changed"}:
            after.append(
                f"{left.element}{_bond(edit.new_order, styling)}{right.element}"
            )
    if not before:
        before.extend(sorted(set(endpoint_elements)))
    if not after:
        after.extend(sorted(set(endpoint_elements)))

    def counted(tokens: Sequence[str]) -> str:
        counts = Counter(tokens)
        parts = []
        for token in sorted(counts):
            count = counts[token]
            parts.append(
                str(
                    load_reaction_label_rendering()["templates"][
                        "counted_clause"
                    ]
                ).format(
                    count=count,
                    times=styling["times"],
                    clause=token,
                )
                if count > 1
                else token
            )
        return " + ".join(parts)

    return f"{counted(before)} {styling['arrow']} {counted(after)}"


def _event_labels(
    events: Sequence[ReactionEvent],
    *,
    style: str,
    prefer_transition: bool = False,
) -> tuple[tuple[str, ...], str, str]:
    """Render and aggregate event-local labels without changing identity."""
    rendering = load_reaction_label_rendering()
    styling = _style(style)
    labels = []
    details = []
    for number, event in enumerate(events, start=1):
        event_clauses = _ordered_clauses(
            render_reaction_label_clause(edit, style=style) for edit in event.edits
        )
        pattern = match_reaction_label_pattern(event.edits, style=style)
        label = (
            _edit_transition_label(event.edits, style=style)
            if prefer_transition
            else (
                pattern.label
                if pattern
                else _compose_concise(event_clauses, style=style)
            )
        )
        labels.append(label)
        details.append(
            str(rendering["templates"]["event_detail"]).format(
                number=number,
                label=label,
                clauses=styling["separator"].join(
                    clause.detailed for clause in event_clauses
                ),
            )
        )
    counts = Counter(labels)
    first_position = {
        label: index for index, label in reversed(list(enumerate(labels)))
    }
    concise_parts = []
    for label in sorted(counts, key=first_position.get):
        count = counts[label]
        concise_parts.append(
            str(rendering["templates"]["counted_clause"]).format(
                count=count,
                times=styling["times"],
                clause=label,
            )
            if count > 1
            else label
        )
    concise = styling["event_separator"].join(concise_parts)
    event_details = styling["event_detail_separator"].join(details)
    detailed = str(rendering["templates"]["multi_event_detail"]).format(
        label=concise,
        events=event_details,
    )
    return (
        tuple(labels),
        concise,
        detailed,
    )


def _render_fragment_indices(text: str, *, style: str) -> str:
    """Render reaction-local generic fragment indices as Unicode superscripts."""
    if _style(style).get("fragment_index_format") != "superscript":
        return text
    symbols = load_fragment_context_symbols()
    if not symbols:
        return text
    symbol_pattern = "|".join(re.escape(symbol) for symbol in symbols)
    superscript_digits = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

    def replace(match: re.Match[str]) -> str:
        return match.group(1) + match.group(2).translate(superscript_digits)

    return re.sub(rf"({symbol_pattern})([0-9]+)", replace, text)


def _render_formula_counts(text: str, *, style: str) -> str:
    """Render molecular-formula counts without changing other numeric context."""
    if _style(style).get("formula_count_format") != "subscript":
        return text
    subscript_digits = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    def replace(match: re.Match[str]) -> str:
        return match.group(1) + match.group(2).translate(subscript_digits)

    return re.sub(r"([A-Z][a-z]?|\))([0-9]+)", replace, text)


def _topology_context(
    topology: Optional[ReactionTopology],
    *,
    rendering: dict[str, Any],
) -> str:
    """Render topology as secondary context rather than a leading phrase."""
    if topology is None or topology.reaction_scope != "intramolecular":
        return ""
    templates = rendering["templates"]
    if len(topology.formed_ring_sizes) == 1:
        return str(templates["intramolecular_ring_context"]).format(
            ring_size=topology.formed_ring_sizes[0]
        )
    return str(templates["intramolecular_context"])


def _make_detailed_label_readable(
    detailed: str,
    *,
    topology: Optional[ReactionTopology],
    style: str,
    rendering: dict[str, Any],
) -> str:
    """Put the transformation first and apply display-only fragment typography."""
    context = _topology_context(topology, rendering=rendering)
    if context:
        scope_prefix = "intramolecular "
        had_leading_scope = detailed.startswith(scope_prefix)
        core = detailed[len(scope_prefix) :] if had_leading_scope else detailed
        if had_leading_scope:
            separator = _style(style)["separator"]
            head, delimiter, tail = core.partition(separator)
            detailed = (
                f"{head}{separator}{context}{separator}{tail}"
                if delimiter
                else f"{head}{separator}{context}"
            )
        else:
            detailed = f"{core}{_style(style)['separator']}{context}"
    detailed = _render_fragment_indices(detailed, style=style)
    return _render_formula_counts(detailed, style=style)


def build_reaction_display_label(
    *,
    edits: Sequence[ReactionEdit],
    selected_label: Optional[str],
    selected_exact: bool,
    grammar_id: Optional[str],
    contextual_label: Optional[ContextualTransformationLabel],
    named_family: Optional[str],
    fallback_label: Optional[str],
    fallback_status: str,
    evidence: str,
    confidence: float,
    events: Sequence[ReactionEvent] = (),
    topology: Optional[ReactionTopology] = None,
    warnings: Iterable[str] = (),
    style: str = "unicode",
    fallback_detailed_label: Optional[str] = None,
) -> Optional[ReactionDisplayLabel]:
    """Build the best display label while retaining its evidence and clauses."""
    rendering = load_reaction_label_rendering()
    styling = _style(style)
    clauses = _ordered_clauses(
        render_reaction_label_clause(edit, style=style) for edit in edits
    )
    concise_clauses = _compose_concise(clauses, style=style) if clauses else ""
    detailed_clauses = styling["separator"].join(
        clause.detailed for clause in clauses
    )
    warning_tuple = tuple(sorted(set(str(warning) for warning in warnings)))
    pattern = match_reaction_label_pattern(edits, style=style) if clauses else None
    inferred_transition = (
        _edit_transition_label(edits, style=style)
        if evidence == "global_atom_correspondence" and clauses
        else None
    )
    structural_label = concise_clauses or None
    transformation_label = pattern.label if pattern else None
    grammar_label = selected_label if selected_exact and selected_label else None
    exact_display_label = selected_label
    if selected_exact and contextual_label is not None and grammar_id is not None:
        from .reaction_labels import load_reaction_rendering

        rule = load_reaction_rendering().get(grammar_id) or {}
        if bool(rule.get("prefer_contextual_label")):
            exact_display_label = contextual_label.concise
    rendered_event_labels: tuple[str, ...] = ()
    if evidence in {
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    } and clauses:
        concise = str(rendering["templates"]["conflict"]).format(
            clauses=concise_clauses
        )
        detailed = str(rendering["templates"]["conflict_detail"]).format(
            label=concise_clauses,
            clauses=detailed_clauses,
        )
        status = "conflicting_evidence"
    elif len(events) > 1 and clauses:
        rendered_event_labels, concise, detailed = _event_labels(
            events,
            style=style,
            prefer_transition=evidence == "global_atom_correspondence",
        )
        transformation_label = concise
        status = "multi_event"
    elif selected_exact and exact_display_label:
        concise = exact_display_label
        detailed = str(rendering["templates"]["exact_detail"]).format(
            label=exact_display_label,
            clauses=detailed_clauses or "none",
        )
        status = "family_overlay" if named_family else "exact_reconstruction"
    elif pattern is not None:
        concise = (
            contextual_label.concise
            if contextual_label
            else inferred_transition or pattern.label
        )
        detail_template = (
            rendering["templates"]["contextual_detail"]
            if contextual_label
            else rendering["templates"]["exact_detail"]
        )
        detailed = str(detail_template).format(
            label=concise,
            transformation=pattern.label,
            clauses=detailed_clauses,
        )
        status = "generic_pattern"
    elif clauses:
        concise = (
            contextual_label.concise
            if contextual_label
            else inferred_transition or concise_clauses
        )
        detailed = (
            str(rendering["templates"]["exact_detail"]).format(
                label=contextual_label.concise,
                clauses=detailed_clauses,
            )
            if contextual_label
            else str(rendering["templates"]["observed_detail"]).format(
                label=concise,
                clauses=detailed_clauses,
            )
        )
        status = "observed_edits"
    elif fallback_label:
        concise = fallback_label
        detailed = fallback_detailed_label or fallback_label
        status = (
            fallback_status
            if fallback_status
            in {
                "partial_product_correspondence",
                "reactant_only",
                "ambiguous_reactants",
            }
            else "unavailable"
        )
        if status == "partial_product_correspondence":
            structural_label = fallback_label
            transformation_label = fallback_label
    else:
        return None
    prefix = topology_label_prefix(topology)
    if prefix:
        scope_prefix = "intramolecular "
        concise = (
            prefix + concise[len(scope_prefix) :]
            if concise.startswith(scope_prefix)
            else prefix + concise
        )
    detailed = _make_detailed_label_readable(
        detailed,
        topology=topology,
        style=style,
        rendering=rendering,
    )
    return ReactionDisplayLabel(
        concise=concise,
        detailed=detailed,
        status=status,
        clauses=clauses,
        evidence=evidence,
        confidence=confidence,
        warnings=warning_tuple,
        style=style,
        definition_version=str(rendering["label_schema_version"]),
        structural_label=structural_label,
        transformation_label=transformation_label,
        grammar_label=grammar_label,
        family_label=named_family,
        pattern_id=pattern.pattern_id if pattern else None,
        pattern_definition_version=pattern.definition_version if pattern else None,
        grammar_id=grammar_id if grammar_label else None,
        contextual_label=contextual_label.concise if contextual_label else None,
        reactant_context_label=contextual_label.before if contextual_label else None,
        product_context_label=contextual_label.after if contextual_label else None,
        event_labels=rendered_event_labels,
        event_count=len(events),
    )


__all__ = [
    "build_reaction_display_label",
    "load_reaction_label_rendering",
    "render_reaction_label_clause",
]
