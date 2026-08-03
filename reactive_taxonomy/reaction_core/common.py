"""Private shared types and atom identities for reaction-core modules."""

from __future__ import annotations

from typing import Any, Mapping, Optional, Tuple

from ..reaction_models import ReactionAtomReference, ReactionComponent


AtomIdentity = Tuple[object, ...]
Coordinate = Tuple[int, int]
Location = Tuple[ReactionComponent, Any, int]
EditRecord = Tuple[int, str, Tuple[AtomIdentity, ...]]


def atom_identity(reference: ReactionAtomReference) -> AtomIdentity:
    """Return a stable mapped or coordinate identity for an edit atom."""
    if reference.atom_map_number is not None:
        return ("map", int(reference.atom_map_number))
    return (
        "coordinate",
        str(reference.side),
        int(reference.component_index),
        int(reference.atom_index),
        str(reference.element),
    )


def atom_map_number(
    atom: Any,
    *,
    side: str,
    component_index: int,
    atom_index: int,
    atom_map_overrides: Optional[Mapping[tuple[str, int, int], int]],
) -> int:
    """Return a map number without changing parsed atom coordinates."""
    return int(
        (atom_map_overrides or {}).get(
            (side, component_index, atom_index),
            int(atom.GetAtomMapNum()),
        )
    )


__all__ = [
    "AtomIdentity",
    "Coordinate",
    "EditRecord",
    "Location",
    "atom_identity",
    "atom_map_number",
]
