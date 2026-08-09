"""Normalize hydrogen-coupled C-X edits into one executable event."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence


ALLOWED_HYDROGEN_CARRIERS = frozenset({"N", "O", "S"})


@dataclass(frozen=True)
class NormalizedCxEvent:
    """One heavy-atom C-X event plus directly coupled hydrogen losses."""

    heavy_edit_indices: tuple[int, ...]
    hydrogen_edit_indices: tuple[int, ...]
    source_event_ids: tuple[str, ...]


def _event_edits(
    event: Mapping[str, Any],
    edits: Sequence[Mapping[str, Any]],
) -> tuple[tuple[int, Mapping[str, Any]], ...] | None:
    indices = tuple(int(value) for value in event.get("edit_indices") or ())
    if not indices or any(index < 0 or index >= len(edits) for index in indices):
        return None
    return tuple((index, edits[index]) for index in indices)


def _is_supported_hydrogen_loss(edit: Mapping[str, Any]) -> bool:
    atom = edit.get("atom_1") or {}
    return (
        edit.get("edit_type") == "hydrogen_change"
        and str(atom.get("element") or "") in ALLOWED_HYDROGEN_CARRIERS
        and str(edit.get("old_order") or "").upper() == "SINGLE"
        and edit.get("new_order") is None
    )


def _locally_related(
    heavy_event_id: str,
    hydrogen_event_id: str,
    relations: Sequence[Mapping[str, Any]],
) -> bool:
    pair = {heavy_event_id, hydrogen_event_id}
    for relation in relations:
        if {
            str(relation.get("event_id_1") or ""),
            str(relation.get("event_id_2") or ""),
        } != pair:
            continue
        paths = relation.get("shortest_paths") or ()
        return any(
            path.get("bond_count") is not None
            and 0 <= int(path["bond_count"]) <= 1
            for path in paths
        )
    return False


def normalize_single_cx_event(
    core: Mapping[str, Any],
    observation: Mapping[str, Any],
) -> NormalizedCxEvent | None:
    """Coalesce local N/O/S-H loss with one heavy-atom transformation event.

    Reaction observation remains unchanged. This interpretation accepts the
    common case where inferred tautomer correspondence places an X-H loss in a
    separate event at the same heteroatom or an adjacent aromatic heteroatom.
    Unrelated hydrogen chemistry and multiple heavy-atom events remain blocked.
    """

    events = tuple(core.get("events") or ())
    edits = tuple(observation.get("edits") or ())
    if int(core.get("event_count") or 0) != len(events) or not events or not edits:
        return None

    heavy_events: list[tuple[Mapping[str, Any], tuple[Any, ...]]] = []
    hydrogen_events: list[tuple[Mapping[str, Any], tuple[Any, ...]]] = []
    for event in events:
        indexed = _event_edits(event, edits)
        if indexed is None:
            return None
        if any(edit.get("edit_type") != "hydrogen_change" for _, edit in indexed):
            heavy_events.append((event, indexed))
        else:
            hydrogen_events.append((event, indexed))
    if len(heavy_events) != 1:
        return None

    heavy_event, heavy_indexed = heavy_events[0]
    heavy_id = str(heavy_event.get("event_id") or "")
    if not heavy_id:
        return None
    if not all(
        _is_supported_hydrogen_loss(edit)
        for _, edit in heavy_indexed
        if edit.get("edit_type") == "hydrogen_change"
    ):
        return None
    relations = tuple(core.get("event_relations") or ())
    hydrogen_indices: list[int] = [
        index
        for index, edit in heavy_indexed
        if edit.get("edit_type") == "hydrogen_change"
    ]
    source_ids = [heavy_id]
    for event, indexed in hydrogen_events:
        event_id = str(event.get("event_id") or "")
        if (
            not event_id
            or not all(_is_supported_hydrogen_loss(edit) for _, edit in indexed)
            or not _locally_related(heavy_id, event_id, relations)
        ):
            return None
        source_ids.append(event_id)
        hydrogen_indices.extend(index for index, _ in indexed)

    return NormalizedCxEvent(
        heavy_edit_indices=tuple(
            index
            for index, edit in heavy_indexed
            if edit.get("edit_type") != "hydrogen_change"
        ),
        hydrogen_edit_indices=tuple(sorted(hydrogen_indices)),
        source_event_ids=tuple(source_ids),
    )


__all__ = ["NormalizedCxEvent", "normalize_single_cx_event"]
