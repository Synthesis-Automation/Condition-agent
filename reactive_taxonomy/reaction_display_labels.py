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
from .reaction_models import (
    RenderedReactionLabel,
    ReactionComponent,
    ReactionCoreProjection,
    ReactionEdit,
    ReactionEvent,
    ReactionLabelClause,
    ReactionTopology,
)
from .reaction_ring_rendering import render_ring_change
from .reaction_topology import topology_label_prefix

_PATH = Path(__file__).with_name("definitions") / "reaction_label_rendering.v1.json"


@lru_cache(maxsize=1)
def load_reaction_label_rendering() -> dict[str, Any]:
    """Load the versioned declarative edit-label rendering rules."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != "2.0":
        raise ValueError("Unsupported reaction-label rendering schema")
    return dict(payload)


def _fragment_context_symbols() -> tuple[str, ...]:
    """Return generic fragments eligible for reaction-local indexing."""
    values = (
        load_reaction_label_rendering().get("fragment_context_symbols") or ()
    )
    return tuple(
        sorted(
            {str(value) for value in values},
            key=lambda value: (-len(value), value),
        )
    )


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
    symbols = _fragment_context_symbols()
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


def _render_core_text(text: str, *, style: str) -> str:
    """Apply display style to presentation-only reaction-core text."""
    if style == "unicode":
        return text
    return (
        text.replace("→", "->")
        .replace("–", "-")
        .replace("≡", "#")
        .replace("×", "x")
        .replace("∅", "none")
    )


def _core_context_tokens(label: str) -> frozenset[str]:
    """Return broad substituent classes that make a core label informative."""
    return frozenset(
        match.group(1)
        for match in re.finditer(
            r"(?<![A-Za-z])(Cycloalkyl|HetAr|Alkenyl|Alkynyl|Acyl|Ar|R)"
            r"(?![a-z])",
            label,
        )
    )


def _core_display_parts(
    reaction_core: Optional[ReactionCoreProjection],
    *,
    style: str,
    rendering: dict[str, Any],
) -> Optional[tuple[str, str]]:
    """Render a non-blocked core as a concise equation and typed audit detail."""
    if reaction_core is None or reaction_core.quality.status == "blocked":
        return None
    templates = rendering["templates"]
    styling = _style(style)
    presentation = reaction_core.presentation
    atom_equation = _render_core_text(presentation.equation, style=style)
    concise = _render_core_text(reaction_core.generic_label, style=style)
    sections = [
        str(templates["core_projection"]).format(equation=atom_equation),
    ]
    typed_sections = (
        ("core_bond_changes", "changes", presentation.bond_changes),
        ("core_atom_state_changes", "changes", presentation.atom_state_changes),
        ("core_retained_context", "contexts", presentation.retained_context),
        ("core_departing_context", "contexts", presentation.departing_context),
        ("core_appearing_context", "contexts", presentation.appearing_context),
    )
    for template_key, value_key, values in typed_sections:
        if values:
            sections.append(
                str(templates[template_key]).format(
                    **{
                        value_key: ", ".join(
                            _render_core_text(value, style=style)
                            for value in values
                        )
                    }
                )
            )
    sections.extend(
        (
            str(templates["core_evidence"]).format(
                evidence=presentation.evidence_label
            ),
            str(templates["core_quality"]).format(
                quality=presentation.quality_label
            ),
        )
    )
    return concise, styling["separator"].join(sections)


def _core_primary_label(
    reaction_core: ReactionCoreProjection,
    core_concise: str,
    *,
    rendering: dict[str, Any],
) -> str:
    """Mark review-only or inferred core equations as provisional."""
    if (
        reaction_core.evidence_status == "verified"
        and reaction_core.quality.status == "pass"
    ):
        return core_concise
    return str(rendering["templates"]["provisional_core"]).format(
        evidence=reaction_core.presentation.evidence_label,
        core=core_concise,
    )


def _core_adds_context(base_label: str, core_label: str) -> bool:
    """Whether the core retains substituent classes absent from a base label."""
    return bool(_core_context_tokens(core_label) - _core_context_tokens(base_label))


def _core_should_lead(
    reaction_core: Optional[ReactionCoreProjection],
    *,
    base_label: str,
) -> bool:
    """Prefer a trustworthy, materially richer core over a generic edit label."""
    if reaction_core is None or reaction_core.quality.status == "blocked":
        return False
    return (
        reaction_core.quality.status == "pass"
        and _core_adds_context(base_label, reaction_core.generic_label)
    )


def _build_reaction_label(
    *,
    reactants: Sequence[ReactionComponent],
    edits: Sequence[ReactionEdit],
    contextual_label: Optional[ContextualTransformationLabel],
    fallback_label: Optional[str],
    fallback_status: str,
    evidence: str,
    confidence: float,
    events: Sequence[ReactionEvent] = (),
    topology: Optional[ReactionTopology] = None,
    reaction_core: Optional[ReactionCoreProjection] = None,
    interpretation_pattern_label: Optional[str] = None,
    interpretation_pattern_id: Optional[str] = None,
    warnings: Iterable[str] = (),
    style: str = "unicode",
    fallback_detailed_label: Optional[str] = None,
) -> Optional[RenderedReactionLabel]:
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
    rendered_event_labels: tuple[str, ...] = ()
    ring_change = next(
        (
            change
            for change in (topology.ring_changes if topology is not None else ())
            if len(change.formed_bond_types) >= 2
        ),
        None,
    )
    ring_display = (
        render_ring_change(
            ring_change,
            reactants=reactants,
            style=styling,
            templates=rendering["templates"],
            raw_edit_audit=detailed_clauses or "none",
        )
        if ring_change is not None
        else None
    )
    core_display = _core_display_parts(
        reaction_core,
        style=style,
        rendering=rendering,
    )
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
        source = "literal_edits"
    elif len(events) > 1 and clauses:
        rendered_event_labels, concise, detailed = _event_labels(
            events,
            style=style,
            prefer_transition=evidence == "global_atom_correspondence",
        )
        transformation_label = concise
        status = "multi_event"
        source = "literal_edits"
    elif ring_display is not None:
        concise = ring_display.concise
        detailed = ring_display.detailed
        structural_label = ring_display.concise
        transformation_label = ring_display.concise
        status = "ring_formation"
        source = "generic_topology"
    elif interpretation_pattern_label and clauses:
        concise = interpretation_pattern_label
        detailed = str(rendering["templates"]["exact_detail"]).format(
            label=concise,
            clauses=detailed_clauses,
        )
        transformation_label = concise
        status = "generic_pattern"
        source = "optional_pattern"
    elif pattern is not None:
        pattern_concise = (
            contextual_label.concise
            if contextual_label
            else inferred_transition or pattern.label
        )
        if (
            core_display is not None
            and _core_should_lead(reaction_core, base_label=pattern_concise)
        ):
            assert reaction_core is not None
            core_concise, core_detailed = core_display
            concise = _core_primary_label(
                reaction_core,
                core_concise,
                rendering=rendering,
            )
            detailed = styling["separator"].join(
                (
                    concise,
                    core_detailed,
                    str(rendering["templates"]["pattern_overlay"]).format(
                        label=pattern_concise
                    ),
                    str(rendering["templates"]["literal_edit_audit"]).format(
                        clauses=detailed_clauses
                    ),
                )
            )
            structural_label = reaction_core.generic_label
            transformation_label = core_concise
            status = "core_projection"
            source = "reaction_core"
        else:
            concise = pattern_concise
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
            source = "literal_edits"
    elif clauses:
        observed_concise = (
            contextual_label.concise
            if contextual_label
            else inferred_transition or concise_clauses
        )
        if (
            core_display is not None
            and _core_should_lead(reaction_core, base_label=observed_concise)
        ):
            assert reaction_core is not None
            core_concise, core_detailed = core_display
            concise = _core_primary_label(
                reaction_core,
                core_concise,
                rendering=rendering,
            )
            detailed = styling["separator"].join(
                (
                    concise,
                    core_detailed,
                    str(rendering["templates"]["literal_edit_audit"]).format(
                        clauses=detailed_clauses
                    ),
                )
            )
            structural_label = reaction_core.generic_label
            transformation_label = core_concise
            status = "core_projection"
            source = "reaction_core"
        else:
            concise = observed_concise
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
            source = "literal_edits"
    elif fallback_status == "partial_product_correspondence" and fallback_label:
        concise = fallback_label
        detailed = fallback_detailed_label or fallback_label
        structural_label = fallback_label
        transformation_label = fallback_label
        status = "partial_product_correspondence"
        source = "partial_product_correspondence"
    elif core_display is not None:
        assert reaction_core is not None
        core_concise, core_detailed = core_display
        concise = _core_primary_label(
            reaction_core,
            core_concise,
            rendering=rendering,
        )
        detailed = styling["separator"].join((concise, core_detailed))
        structural_label = reaction_core.generic_label
        transformation_label = core_concise
        status = "core_projection"
        source = "reaction_core"
    elif fallback_label:
        concise = fallback_label
        detailed = fallback_detailed_label or fallback_label
        status = "unavailable"
        source = "unavailable"
    else:
        return None
    if (
        core_display is not None
        and source != "reaction_core"
        and status != "conflicting_evidence"
        and reaction_core is not None
        and _core_adds_context(concise, core_display[0])
    ):
        detailed = styling["separator"].join((detailed, core_display[1]))
    label_evidence = evidence
    label_confidence = confidence
    if source == "reaction_core" and reaction_core is not None:
        label_evidence = reaction_core.evidence
        label_confidence = min(confidence, reaction_core.confidence)
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
    return RenderedReactionLabel(
        concise=concise,
        detailed=detailed,
        status=status,
        source=source,
        clauses=clauses,
        evidence=label_evidence,
        confidence=label_confidence,
        warnings=warning_tuple,
        style=style,
        definition_version=str(rendering["label_schema_version"]),
        structural_label=structural_label,
        transformation_label=transformation_label,
        pattern_id=(
            interpretation_pattern_id
            if interpretation_pattern_label
            else pattern.pattern_id if pattern else None
        ),
        pattern_definition_version=(
            "2.0"
            if interpretation_pattern_label
            else pattern.definition_version if pattern else None
        ),
        contextual_label=contextual_label.concise if contextual_label else None,
        reactant_context_label=contextual_label.before if contextual_label else None,
        product_context_label=contextual_label.after if contextual_label else None,
        event_labels=rendered_event_labels,
        event_count=len(events),
    )


__all__ = [
    "load_reaction_label_rendering",
    "render_reaction_label_clause",
]
