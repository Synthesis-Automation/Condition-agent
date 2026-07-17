"""Definition-driven interpretation of normalized reaction-edit patterns."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

from .reaction_models import ReactionEdit


_PATH = Path(__file__).with_name("definitions") / "reaction_label_patterns.v1.json"


@dataclass(frozen=True)
class ReactionLabelPatternMatch:
    """One deterministic generic transformation-label interpretation."""

    pattern_id: str
    label: str
    definition_version: str


@lru_cache(maxsize=1)
def load_reaction_label_patterns() -> dict[str, Any]:
    """Load versioned generic reaction-label pattern definitions."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != "1.0":
        raise ValueError("Unsupported reaction-label pattern schema")
    return dict(payload)


def _atom_key(atom: Any) -> tuple[Any, ...]:
    if atom.atom_map_number is not None:
        return ("map", int(atom.atom_map_number))
    return ("atom", int(atom.component_index), int(atom.atom_index), atom.element)


def _edit_atom_keys(edit: ReactionEdit) -> frozenset[tuple[Any, ...]]:
    return frozenset(
        _atom_key(atom) for atom in (edit.atom_1, edit.atom_2) if atom is not None
    )


def _edits_of_type(edits: Sequence[ReactionEdit], edit_type: str) -> list[ReactionEdit]:
    return [edit for edit in edits if edit.edit_type == edit_type]


def _substitution(edits: Sequence[ReactionEdit]) -> Optional[dict[str, ReactionEdit]]:
    formed = _edits_of_type(edits, "formed")
    broken = _edits_of_type(edits, "broken")
    changed = _edits_of_type(edits, "order_changed")
    if len(formed) != 1 or len(broken) != 1 or changed:
        return None
    if not (_edit_atom_keys(formed[0]) & _edit_atom_keys(broken[0])):
        return None
    return {"formed": formed[0]}


def _hydrogenation(edits: Sequence[ReactionEdit]) -> Optional[dict[str, ReactionEdit]]:
    changed = _edits_of_type(edits, "order_changed")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    if len(changed) != 1 or changed[0].old_order != "DOUBLE" or changed[0].new_order != "SINGLE":
        return None
    if {changed[0].atom_1.element, changed[0].atom_2.element} != {"C"}:
        return None
    gains = [edit for edit in hydrogen if edit.old_order is None and edit.new_order == "SINGLE"]
    if len(gains) != 2 or any(edit.edit_type in {"formed", "broken"} for edit in edits):
        return None
    endpoints = _edit_atom_keys(changed[0])
    if any(_atom_key(edit.atom_1) not in endpoints for edit in gains):
        return None
    return {"changed": changed[0]}


def _alkyne_hydrogenation(
    edits: Sequence[ReactionEdit], *, partial: bool
) -> Optional[dict[str, ReactionEdit]]:
    changed = _edits_of_type(edits, "order_changed")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    expected_order = "DOUBLE" if partial else "SINGLE"
    expected_gains = 2 if partial else 4
    if (
        len(changed) != 1
        or changed[0].old_order != "TRIPLE"
        or changed[0].new_order != expected_order
        or {changed[0].atom_1.element, changed[0].atom_2.element} != {"C"}
    ):
        return None
    gains = [
        edit
        for edit in hydrogen
        if edit.old_order is None and edit.new_order == "SINGLE"
    ]
    endpoints = _edit_atom_keys(changed[0])
    if (
        len(gains) != expected_gains
        or any(edit.edit_type in {"formed", "broken"} for edit in edits)
        or any(_atom_key(edit.atom_1) not in endpoints for edit in gains)
    ):
        return None
    return {"changed": changed[0]}


def _complete_alkyne_hydrogenation(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    return _alkyne_hydrogenation(edits, partial=False)


def _partial_alkyne_hydrogenation(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    return _alkyne_hydrogenation(edits, partial=True)


def _heteroatom_bond_reduction(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    changed = _edits_of_type(edits, "order_changed")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    if (
        len(changed) != 1
        or changed[0].old_order != "DOUBLE"
        or changed[0].new_order != "SINGLE"
        or changed[0].atom_1.element == changed[0].atom_2.element == "C"
    ):
        return None
    gains = [
        edit
        for edit in hydrogen
        if edit.old_order is None and edit.new_order == "SINGLE"
    ]
    endpoints = _edit_atom_keys(changed[0])
    if (
        len(gains) != 2
        or any(edit.edit_type in {"formed", "broken"} for edit in edits)
        or any(_atom_key(edit.atom_1) not in endpoints for edit in gains)
    ):
        return None
    return {"changed": changed[0]}


def _dehydrogenation(edits: Sequence[ReactionEdit]) -> Optional[dict[str, ReactionEdit]]:
    changed = _edits_of_type(edits, "order_changed")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    if len(changed) != 1 or changed[0].old_order != "SINGLE" or changed[0].new_order != "DOUBLE":
        return None
    if {changed[0].atom_1.element, changed[0].atom_2.element} != {"C"}:
        return None
    losses = [edit for edit in hydrogen if edit.old_order == "SINGLE" and edit.new_order is None]
    if len(losses) != 2 or any(edit.edit_type in {"formed", "broken"} for edit in edits):
        return None
    endpoints = _edit_atom_keys(changed[0])
    if any(_atom_key(edit.atom_1) not in endpoints for edit in losses):
        return None
    return {"changed": changed[0]}


def _heteroatom_bond_oxidation(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    changed = _edits_of_type(edits, "order_changed")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    if (
        len(changed) != 1
        or changed[0].old_order != "SINGLE"
        or changed[0].new_order != "DOUBLE"
        or changed[0].atom_1.element == changed[0].atom_2.element == "C"
    ):
        return None
    losses = [
        edit
        for edit in hydrogen
        if edit.old_order == "SINGLE" and edit.new_order is None
    ]
    endpoints = _edit_atom_keys(changed[0])
    if (
        len(losses) != 2
        or any(edit.edit_type in {"formed", "broken"} for edit in edits)
        or any(_atom_key(edit.atom_1) not in endpoints for edit in losses)
    ):
        return None
    return {"changed": changed[0]}


def _reductive_bond_cleavage(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    broken = _edits_of_type(edits, "broken")
    hydrogen = _edits_of_type(edits, "hydrogen_change")
    if len(broken) != 1 or any(
        edit.edit_type in {"formed", "order_changed"} for edit in edits
    ):
        return None
    gains = [
        edit
        for edit in hydrogen
        if edit.old_order is None and edit.new_order == "SINGLE"
    ]
    endpoints = _edit_atom_keys(broken[0])
    if len(gains) < 2 or any(
        _atom_key(edit.atom_1) not in endpoints for edit in gains
    ):
        return None
    return {"broken": broken[0]}


def _intramolecular_bond_formation(
    edits: Sequence[ReactionEdit],
) -> Optional[dict[str, ReactionEdit]]:
    formed = _edits_of_type(edits, "formed")
    non_hydrogen = [edit for edit in edits if edit.edit_type != "hydrogen_change"]
    if len(formed) != 1 or len(non_hydrogen) != 1 or formed[0].atom_2 is None:
        return None
    if formed[0].atom_1.component_index != formed[0].atom_2.component_index:
        return None
    return {"formed": formed[0]}


_MATCHERS: dict[
    str, Callable[[Sequence[ReactionEdit]], Optional[dict[str, ReactionEdit]]]
] = {
    "substitution": _substitution,
    "hydrogenation": _hydrogenation,
    "complete_alkyne_hydrogenation": _complete_alkyne_hydrogenation,
    "partial_alkyne_hydrogenation": _partial_alkyne_hydrogenation,
    "heteroatom_bond_reduction": _heteroatom_bond_reduction,
    "dehydrogenation": _dehydrogenation,
    "heteroatom_bond_oxidation": _heteroatom_bond_oxidation,
    "reductive_bond_cleavage": _reductive_bond_cleavage,
    "intramolecular_bond_formation": _intramolecular_bond_formation,
}


def _ordered_elements(edit: ReactionEdit) -> tuple[str, str]:
    elements = sorted(
        (edit.atom_1.element, edit.atom_2.element if edit.atom_2 else "H"),
        key=lambda element: (element != "C", element),
    )
    return elements[0], elements[1]


def _bond_symbol(order: Optional[str], style: str) -> str:
    symbols = {
        "unicode": {"SINGLE": "–", "DOUBLE": "=", "TRIPLE": "≡", "AROMATIC": ":"},
        "ascii": {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": ":"},
        "hte_legacy": {"SINGLE": "-", "DOUBLE": "=", "TRIPLE": "#", "AROMATIC": ":"},
    }
    return symbols[style].get(str(order or "SINGLE").upper(), symbols[style]["SINGLE"])


def _pair(edit: ReactionEdit, order: Optional[str], style: str) -> str:
    left, right = _ordered_elements(edit)
    return f"{left}{_bond_symbol(order, style)}{right}"


def match_reaction_label_pattern(
    edits: Sequence[ReactionEdit], *, style: str = "unicode"
) -> Optional[ReactionLabelPatternMatch]:
    """Return the highest-priority generic interpretation supported by edits."""
    definitions = load_reaction_label_patterns()
    patterns = sorted(
        definitions.get("patterns") or [],
        key=lambda item: (-int(item.get("priority", 0)), str(item.get("id", ""))),
    )
    for pattern in patterns:
        matcher = _MATCHERS.get(str(pattern.get("matcher") or ""))
        matched = matcher(edits) if matcher else None
        if not matched:
            continue
        values: dict[str, str] = {}
        if "formed" in matched:
            values["formed_pair"] = _pair(
                matched["formed"], matched["formed"].new_order, style
            )
        if "changed" in matched:
            values["old_bond_pair"] = _pair(
                matched["changed"], matched["changed"].old_order, style
            )
            values["new_bond_pair"] = _pair(
                matched["changed"], matched["changed"].new_order, style
            )
        if "broken" in matched:
            values["broken_pair"] = _pair(
                matched["broken"], matched["broken"].old_order, style
            )
        templates = pattern.get("templates") or {}
        template = str(templates.get(style) or "")
        if not template:
            continue
        return ReactionLabelPatternMatch(
            pattern_id=str(pattern["id"]),
            label=template.format(**values),
            definition_version=str(definitions["schema_version"]),
        )
    return None


__all__ = [
    "ReactionLabelPatternMatch",
    "load_reaction_label_patterns",
    "match_reaction_label_pattern",
]
