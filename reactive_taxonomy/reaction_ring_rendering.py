"""Display-only rendering of graph-derived ring-change observations."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    ReactionComponent,
    ReactionRingChange,
)


@dataclass(frozen=True)
class RingChangeDisplay:
    """Concise and auditable text derived from one ring observation."""

    concise: str
    detailed: str


def _formula(elements: Sequence[str], *, subscript: bool) -> str:
    counts = Counter(elements)
    order = sorted(counts, key=lambda value: (value != "C", value))
    translate = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    parts = []
    for element in order:
        count = counts[element]
        suffix = str(count) if count > 1 else ""
        parts.append(element + (suffix.translate(translate) if subscript else suffix))
    return "".join(parts)


def _cycle_segments(change: ReactionRingChange) -> tuple[tuple[Any, ...], ...]:
    segments: list[list[Any]] = []
    for reference in change.atom_references:
        if not segments or segments[-1][-1].component_index != reference.component_index:
            segments.append([reference])
        else:
            segments[-1].append(reference)
    if (
        len(segments) > 1
        and segments[0][0].component_index == segments[-1][0].component_index
    ):
        segments[0] = [*segments[-1], *segments[0]]
        segments.pop()
    return tuple(tuple(segment) for segment in segments)


def _fragment_text(
    segment: Sequence[Any],
    components: Mapping[int, ReactionComponent],
    style: Mapping[str, str],
) -> str:
    if not segment:
        return ""
    component = components.get(int(segment[0].component_index))
    molecule = parse_smiles(component.input_smiles) if component else None
    parts = [str(segment[0].element)]
    for left, right in zip(segment, segment[1:]):
        order = None
        if molecule is not None:
            bond = molecule.GetBondBetweenAtoms(
                int(left.atom_index), int(right.atom_index)
            )
            order = str(bond.GetBondType()).upper() if bond is not None else None
        parts.append(
            {
                "SINGLE": style["single"],
                "DOUBLE": style["double"],
                "TRIPLE": style["triple"],
                "AROMATIC": style["aromatic"],
            }.get(str(order or "").upper(), "…")
        )
        parts.append(str(right.element))
    return "".join(parts)


def _formed_bond_text(
    formed_bond_types: Sequence[str],
    *,
    style: Mapping[str, str],
) -> str:
    counts = Counter(formed_bond_types)
    values = []
    for bond_type, count in sorted(counts.items()):
        rendered = bond_type.replace("-", style["single"])
        values.append(
            f"{count} {style['times']} {rendered} bonds formed"
            if count > 1
            else f"{rendered} bond formed"
        )
    return " + ".join(values)


def render_ring_change(
    change: ReactionRingChange,
    *,
    reactants: Sequence[ReactionComponent],
    style: Mapping[str, str],
    templates: Mapping[str, str],
    raw_edit_audit: str,
) -> RingChangeDisplay:
    """Render a ring observation without changing its graph-derived facts."""
    components = {component.component_index: component for component in reactants}
    fragments = sorted(
        filter(
            None,
            (
                _fragment_text(segment, components, style)
                for segment in _cycle_segments(change)
            ),
        )
    )
    concise = str(templates["ring_formation_concise"]).format(
        fragments=str(style["event_separator"]).join(fragments),
        arrow=style["arrow"],
        aromatic="aromatic " if change.aromatic_after else "",
        ring_size=change.ring_size,
        formula=_formula(
            change.element_sequence,
            subscript=style.get("formula_count_format") == "subscript",
        ),
    )
    electronic = (
        str(templates["ring_aromatization"]).format(arrow=style["arrow"])
        if change.aromatic_after
        else ""
    )
    detailed = str(templates["ring_formation_detail"]).format(
        component_count=len(change.source_component_indices),
        concise=concise,
        formed_bonds=_formed_bond_text(change.formed_bond_types, style=style),
        electronic=electronic,
        evidence=change.evidence.replace("_", " "),
        confidence=f"{change.confidence:.2f}",
        raw_edits=raw_edit_audit,
    )
    return RingChangeDisplay(concise=concise, detailed=detailed)


__all__ = ["RingChangeDisplay", "render_ring_change"]
