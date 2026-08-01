"""Deterministic net-edit archetype classification."""

from __future__ import annotations

from typing import Iterable, Optional, Tuple, cast

from .reaction_models import EditArchetype, ReactionEdit


_ORDER_RANK = {
    "SINGLE": 1,
    "DOUBLE": 2,
    "TRIPLE": 3,
    "AROMATIC": 4,
}


def _atom_key(edit: ReactionEdit, position: int) -> tuple[object, ...] | None:
    atom = edit.atom_1 if position == 1 else edit.atom_2
    if atom is None:
        return None
    if atom.atom_map_number is not None:
        return ("map", int(atom.atom_map_number))
    return (
        atom.side,
        int(atom.component_index),
        int(atom.atom_index),
    )


def _order_direction(edit: ReactionEdit) -> int:
    if edit.edit_type != "order_changed":
        return 0
    old = _ORDER_RANK.get(str(edit.old_order or "").upper())
    new = _ORDER_RANK.get(str(edit.new_order or "").upper())
    if old is None or new is None:
        return 0
    return (new > old) - (new < old)


def infer_edit_archetype(edits: Iterable[ReactionEdit]) -> EditArchetype:
    """Classify normalized edits by their net molecular-graph topology."""
    normalized = tuple(edits)
    if not normalized:
        return "unresolved"

    decreased = tuple(edit for edit in normalized if _order_direction(edit) < 0)
    increased = tuple(edit for edit in normalized if _order_direction(edit) > 0)
    formed = tuple(edit for edit in normalized if edit.edit_type == "formed")
    broken = tuple(edit for edit in normalized if edit.edit_type == "broken")
    hydrogen_added = tuple(
        edit
        for edit in normalized
        if edit.edit_type == "hydrogen_change"
        and edit.old_order is None
        and edit.new_order is not None
    )
    hydrogen_removed = tuple(
        edit
        for edit in normalized
        if edit.edit_type == "hydrogen_change"
        and edit.old_order is not None
        and edit.new_order is None
    )

    if decreased and not increased and (formed or hydrogen_added):
        return "addition"
    if increased and not decreased and (broken or hydrogen_removed):
        return "elimination"

    broken_atoms = {
        key
        for edit in broken
        for key in (_atom_key(edit, 1), _atom_key(edit, 2))
        if key is not None
    }
    formed_atoms = {
        key
        for edit in formed
        for key in (_atom_key(edit, 1), _atom_key(edit, 2))
        if key is not None
    }
    if broken and formed and broken_atoms & formed_atoms:
        return "substitution"
    if (decreased or increased) and not (formed or broken):
        return "bond_order_change"
    if len({edit.edit_type for edit in normalized}) > 1:
        return "composite"
    return "unresolved"


def reconcile_edit_archetype(
    edits: Iterable[ReactionEdit],
    declared: Optional[str],
) -> Tuple[EditArchetype, Tuple[str, ...]]:
    """Prefer observed edits and expose interpretation disagreement."""
    inferred = infer_edit_archetype(edits)
    declared_value = (
        cast(EditArchetype, declared)
        if declared in {
            "substitution",
            "addition",
            "elimination",
            "bond_order_change",
            "composite",
            "unresolved",
        }
        else "unresolved"
    )
    if inferred == "unresolved":
        return declared_value, ()
    if declared_value in {"unresolved", inferred}:
        return inferred, ()
    return inferred, (
        f"EDIT_ARCHETYPE_CONFLICT:declared={declared_value}:observed={inferred}",
    )


__all__ = [
    "infer_edit_archetype",
    "reconcile_edit_archetype",
]
