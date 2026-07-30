"""Anonymous, grammar-independent reaction edit prototypes."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from typing import Any, Mapping, Tuple


ANONYMOUS_EDIT_PROTOTYPE_VERSION = "1.0"


def _bond_pair(token: Any) -> str:
    value = str(token or "")
    return value.split(":", 1)[0]


def _order_pair(token: Any) -> str:
    value = str(token or "")
    return value.split(":", 1)[0]


def _counter_tokens(values: Any, parser: Any) -> Tuple[str, ...]:
    return tuple(sorted(parser(value) for value in values or () if str(value)))


def _build_prototype(
    *,
    formed: Tuple[str, ...],
    broken: Tuple[str, ...],
    changed: Tuple[str, ...],
    topology: Mapping[str, Any],
    event_count: int,
) -> "AnonymousEditPrototype | None":
    if not formed and not broken and not changed:
        return None
    payload = {
        "formed_element_pairs": formed,
        "broken_element_pairs": broken,
        "order_changed_element_pairs": changed,
        "ring_count_delta": int(topology.get("ring_count_delta") or 0),
        "formed_ring_sizes": tuple(
            sorted(int(value) for value in topology.get("formed_ring_sizes") or ())
        ),
        "event_count": int(event_count),
        "version": ANONYMOUS_EDIT_PROTOTYPE_VERSION,
    }
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    )
    prototype_id = "EP1:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[
        :24
    ]
    return AnonymousEditPrototype(prototype_id=prototype_id, **payload)


@dataclass(frozen=True)
class AnonymousEditPrototype:
    """Coarse edit graph used only after deterministic chemistry gates."""

    prototype_id: str
    formed_element_pairs: Tuple[str, ...]
    broken_element_pairs: Tuple[str, ...]
    order_changed_element_pairs: Tuple[str, ...]
    ring_count_delta: int
    formed_ring_sizes: Tuple[int, ...]
    event_count: int
    version: str = ANONYMOUS_EDIT_PROTOTYPE_VERSION


def anonymous_edit_prototype(
    signature: Mapping[str, Any],
) -> AnonymousEditPrototype | None:
    """Derive a name-free edit prototype from a serialized signature."""
    formed = _counter_tokens(signature.get("formed_bond_types"), _bond_pair)
    broken = _counter_tokens(signature.get("broken_bond_types"), _bond_pair)
    changed = _counter_tokens(signature.get("order_changes"), _order_pair)
    topology = signature.get("topology") or {}
    return _build_prototype(
        formed=formed,
        broken=broken,
        changed=changed,
        topology=topology,
        event_count=int(signature.get("event_count") or 0),
    )


def _hypothesis_atom_key(atom: Mapping[str, Any]) -> tuple[Any, ...]:
    map_number = atom.get("atom_map_number")
    if map_number is not None:
        return ("map", int(map_number))
    return (
        str(atom.get("side") or ""),
        int(atom.get("component_index") or 0),
        int(atom.get("atom_index") or 0),
    )


def _hypothesis_element_pair(edit: Mapping[str, Any]) -> str | None:
    atom_1 = edit.get("atom_1")
    atom_2 = edit.get("atom_2")
    if not isinstance(atom_1, Mapping) or not isinstance(atom_2, Mapping):
        return None
    elements = sorted(
        (str(atom_1.get("element") or ""), str(atom_2.get("element") or ""))
    )
    if not all(elements):
        return None
    return "-".join(elements)


def _hypothesis_event_count(edits: Tuple[Mapping[str, Any], ...]) -> int:
    adjacency: dict[tuple[Any, ...], set[tuple[Any, ...]]] = {}
    for edit in edits:
        atom_1 = edit.get("atom_1")
        atom_2 = edit.get("atom_2")
        if not isinstance(atom_1, Mapping):
            continue
        left = _hypothesis_atom_key(atom_1)
        adjacency.setdefault(left, set())
        if not isinstance(atom_2, Mapping):
            continue
        right = _hypothesis_atom_key(atom_2)
        adjacency.setdefault(right, set())
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(adjacency)
    count = 0
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            current = stack.pop()
            neighbors = adjacency.get(current, set()) & remaining
            remaining.difference_update(neighbors)
            stack.extend(neighbors)
    return count


def anonymous_edit_prototype_from_hypothesis(
    hypothesis: Mapping[str, Any],
) -> AnonymousEditPrototype | None:
    """Derive an EP1 query prototype from one unverified edit hypothesis."""
    edits = tuple(
        edit
        for edit in hypothesis.get("edits") or ()
        if isinstance(edit, Mapping)
    )
    formed = []
    broken = []
    changed = []
    for edit in edits:
        pair = _hypothesis_element_pair(edit)
        if not pair:
            continue
        edit_type = str(edit.get("edit_type") or "")
        if edit_type == "formed":
            formed.append(pair)
        elif edit_type == "broken":
            broken.append(pair)
        elif edit_type == "order_changed":
            changed.append(pair)
    topology = hypothesis.get("topology")
    topology_value = topology if isinstance(topology, Mapping) else {}
    return _build_prototype(
        formed=tuple(sorted(formed)),
        broken=tuple(sorted(broken)),
        changed=tuple(sorted(changed)),
        topology=topology_value,
        event_count=_hypothesis_event_count(edits),
    )


def _multiset_jaccard(left: Tuple[str, ...], right: Tuple[str, ...]) -> float:
    left_counts = Counter(left)
    right_counts = Counter(right)
    union = sum((left_counts | right_counts).values())
    if not union:
        return 1.0
    return sum((left_counts & right_counts).values()) / union


def anonymous_edit_compatible(
    left: AnonymousEditPrototype,
    right: AnonymousEditPrototype,
) -> bool:
    """Apply hard structural gates before approximate prototype scoring."""
    left_formed = set(left.formed_element_pairs)
    right_formed = set(right.formed_element_pairs)
    if not left_formed.intersection(right_formed):
        return False
    left_broken = set(left.broken_element_pairs)
    right_broken = set(right.broken_element_pairs)
    if (
        left_broken
        and right_broken
        and not left_broken.intersection(right_broken)
    ):
        return False
    left_ring_sign = (left.ring_count_delta > 0) - (left.ring_count_delta < 0)
    right_ring_sign = (right.ring_count_delta > 0) - (right.ring_count_delta < 0)
    if left_ring_sign != right_ring_sign:
        return False
    if left_ring_sign and abs(left.ring_count_delta - right.ring_count_delta) > 1:
        return False
    return True


def anonymous_edit_similarity(
    left: AnonymousEditPrototype,
    right: AnonymousEditPrototype,
) -> float:
    """Score compatible edit graphs without family names or display labels."""
    if not anonymous_edit_compatible(left, right):
        return 0.0
    formed = _multiset_jaccard(
        left.formed_element_pairs,
        right.formed_element_pairs,
    )
    broken = _multiset_jaccard(
        left.broken_element_pairs,
        right.broken_element_pairs,
    )
    changed = _multiset_jaccard(
        left.order_changed_element_pairs,
        right.order_changed_element_pairs,
    )
    ring_delta = (
        1.0
        if left.ring_count_delta == right.ring_count_delta
        else 0.5
    )
    ring_sizes = _multiset_jaccard(
        tuple(str(value) for value in left.formed_ring_sizes),
        tuple(str(value) for value in right.formed_ring_sizes),
    )
    event = 1.0 if left.event_count == right.event_count else 0.5
    return round(
        0.35 * formed
        + 0.25 * broken
        + 0.10 * changed
        + 0.15 * ring_delta
        + 0.10 * ring_sizes
        + 0.05 * event,
        6,
    )


__all__ = [
    "ANONYMOUS_EDIT_PROTOTYPE_VERSION",
    "AnonymousEditPrototype",
    "anonymous_edit_compatible",
    "anonymous_edit_prototype",
    "anonymous_edit_prototype_from_hypothesis",
    "anonymous_edit_similarity",
]
