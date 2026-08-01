"""Observation-derived reaction-pattern interpretation.

Patterns consume normalized edits and topology. They never propose atom
correspondence, graph edits, or products.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    ReactionEdit,
    ReactionObservation,
    ReactionPatternMatch,
)


REACTION_PATTERN_DEFINITION_SCHEMA_VERSION = "1.0"
_DEFINITION_FILES = (
    "transformation_patterns.v1.json",
    "synthesis_patterns.v1.json",
)
_ORDER_RANK = {
    "NONE": 0.0,
    "SINGLE": 1.0,
    "AROMATIC": 1.5,
    "DOUBLE": 2.0,
    "TRIPLE": 3.0,
}


@lru_cache(maxsize=1)
def load_reaction_pattern_definitions() -> Tuple[dict[str, Any], ...]:
    """Load and validate optional reaction-pattern definitions."""
    root = Path(__file__).with_name("definitions")
    patterns = []
    identifiers: set[str] = set()
    forbidden = {
        "operator_id",
        "predicted_bond_changes",
        "reconstruction_rule_id",
        "slot_relationships",
        "slots",
    }
    for filename in _DEFINITION_FILES:
        with (root / filename).open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if payload.get("schema_version") != REACTION_PATTERN_DEFINITION_SCHEMA_VERSION:
            raise ValueError(f"Unsupported reaction-pattern schema: {filename}")
        for raw in payload.get("patterns") or ():
            pattern = dict(raw)
            pattern_id = str(pattern.get("id") or "")
            if not pattern_id or pattern_id in identifiers:
                raise ValueError(f"Invalid or duplicate reaction pattern: {pattern_id}")
            if forbidden.intersection(pattern):
                raise ValueError(
                    f"Reaction pattern contains structural execution data: {pattern_id}"
                )
            if pattern.get("tier") not in {"generic", "synthesis"}:
                raise ValueError(f"Invalid reaction-pattern tier: {pattern_id}")
            if not str(pattern.get("matcher") or ""):
                raise ValueError(f"Reaction pattern has no matcher: {pattern_id}")
            identifiers.add(pattern_id)
            patterns.append(pattern)
    return tuple(sorted(patterns, key=lambda item: str(item["id"])))


def _atom_key(edit: ReactionEdit, endpoint: int) -> tuple[int, int] | None:
    atom = edit.atom_1 if endpoint == 1 else edit.atom_2
    if atom is None or atom.side != "reactant":
        return None
    return atom.component_index, atom.atom_index


def _element_pair(edit: ReactionEdit) -> frozenset[str]:
    return frozenset(
        atom.element for atom in (edit.atom_1, edit.atom_2) if atom is not None
    )


def _indices(edits: Sequence[ReactionEdit], predicate: Callable[[ReactionEdit], bool]) -> tuple[int, ...]:
    return tuple(index for index, edit in enumerate(edits) if predicate(edit))


def _order_direction(edit: ReactionEdit) -> int:
    old = _ORDER_RANK.get(str(edit.old_order or "NONE").upper())
    new = _ORDER_RANK.get(str(edit.new_order or "NONE").upper())
    if old is None or new is None:
        return 0
    return (new > old) - (new < old)


def _net_substitution(observation: ReactionObservation) -> tuple[int, ...]:
    edits = observation.edits
    broken = _indices(edits, lambda edit: edit.edit_type == "broken")
    formed = _indices(edits, lambda edit: edit.edit_type == "formed")
    for broken_index in broken:
        broken_atoms = {_atom_key(edits[broken_index], 1), _atom_key(edits[broken_index], 2)}
        for formed_index in formed:
            formed_atoms = {_atom_key(edits[formed_index], 1), _atom_key(edits[formed_index], 2)}
            if (broken_atoms - {None}) & (formed_atoms - {None}):
                return tuple(sorted((broken_index, formed_index)))
    return ()


def _net_elimination(observation: ReactionObservation) -> tuple[int, ...]:
    increased = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "order_changed" and _order_direction(edit) > 0,
    )
    losses = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "broken"
        or (
            edit.edit_type == "hydrogen_change"
            and edit.old_order is not None
            and edit.new_order is None
        ),
    )
    return tuple(sorted(increased + losses)) if increased and losses else ()


def _net_addition(observation: ReactionObservation) -> tuple[int, ...]:
    decreased = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "order_changed" and _order_direction(edit) < 0,
    )
    gains = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed"
        or (
            edit.edit_type == "hydrogen_change"
            and edit.old_order is None
            and edit.new_order is not None
        ),
    )
    return tuple(sorted(decreased + gains)) if decreased and gains else ()


def _net_bond_order_change(observation: ReactionObservation) -> tuple[int, ...]:
    order_changes = _indices(
        observation.edits, lambda edit: edit.edit_type == "order_changed"
    )
    non_order = tuple(
        edit for edit in observation.edits if edit.edit_type != "order_changed"
    )
    return order_changes if order_changes and not non_order else ()


def _net_coupling(observation: ReactionObservation) -> tuple[int, ...]:
    if any(
        edit.edit_type == "order_changed" and _order_direction(edit) < 0
        for edit in observation.edits
    ):
        return ()
    formed = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed"
        and edit.atom_2 is not None
        and edit.atom_1.component_index != edit.atom_2.component_index,
    )
    return formed


def _net_ring_closure(observation: ReactionObservation) -> tuple[int, ...]:
    topology = observation.topology
    if topology is None or not (
        topology.formed_ring_sizes
        or (topology.ring_count_delta is not None and topology.ring_count_delta > 0)
    ):
        return ()
    return _indices(observation.edits, lambda edit: edit.edit_type == "formed")


def _reductive_amination_like(observation: ReactionObservation) -> tuple[int, ...]:
    edits = observation.edits
    broken_co = _indices(
        edits,
        lambda edit: edit.edit_type == "broken"
        and _element_pair(edit) == {"C", "O"}
        and str(edit.old_order).upper() == "DOUBLE",
    )
    formed_cn = _indices(
        edits,
        lambda edit: edit.edit_type == "formed" and _element_pair(edit) == {"C", "N"},
    )
    carbon_h_gain = _indices(
        edits,
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "C"
        and edit.old_order is None
        and edit.new_order is not None,
    )
    nitrogen_h_loss = _indices(
        edits,
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "N"
        and edit.old_order is not None
        and edit.new_order is None,
    )
    if not (broken_co and formed_cn and carbon_h_gain and nitrogen_h_loss):
        return ()
    carbonyl_centers = {
        key
        for index in broken_co
        for key in (_atom_key(edits[index], 1), _atom_key(edits[index], 2))
        if key is not None
        and (edits[index].atom_1 if _atom_key(edits[index], 1) == key else edits[index].atom_2).element == "C"
    }
    formed_centers = {
        key
        for index in formed_cn
        for key in (_atom_key(edits[index], 1), _atom_key(edits[index], 2))
        if key is not None
        and (edits[index].atom_1 if _atom_key(edits[index], 1) == key else edits[index].atom_2).element == "C"
    }
    if not carbonyl_centers.intersection(formed_centers):
        return ()
    return tuple(sorted(broken_co + formed_cn + carbon_h_gain + nitrogen_h_loss))


def _reactant_has_carbonyl_at(observation: ReactionObservation, key: tuple[int, int]) -> bool:
    component_index, atom_index = key
    component = next(
        (item for item in observation.reactants if item.component_index == component_index),
        None,
    )
    molecule = parse_smiles(component.input_smiles) if component is not None else None
    if molecule is None or atom_index >= molecule.GetNumAtoms():
        return False
    atom = molecule.GetAtomWithIdx(atom_index)
    return any(
        neighbor.GetSymbol() == "O"
        and str(molecule.GetBondBetweenAtoms(atom_index, int(neighbor.GetIdx())).GetBondType()).upper() == "DOUBLE"
        for neighbor in atom.GetNeighbors()
    )


def _amide_formation_like(observation: ReactionObservation) -> tuple[int, ...]:
    formed_cn = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed" and _element_pair(edit) == {"C", "N"},
    )
    for index in formed_cn:
        edit = observation.edits[index]
        for endpoint in (1, 2):
            atom = edit.atom_1 if endpoint == 1 else edit.atom_2
            key = _atom_key(edit, endpoint)
            if not (
                atom is not None
                and atom.element == "C"
                and key is not None
                and _reactant_has_carbonyl_at(observation, key)
            ):
                continue
            leaving = tuple(
                broken_index
                for broken_index, broken_edit in enumerate(observation.edits)
                if broken_edit.edit_type == "broken"
                and str(broken_edit.old_order or "").upper() == "SINGLE"
                and key
                in {
                    _atom_key(broken_edit, 1),
                    _atom_key(broken_edit, 2),
                }
            )
            if leaving:
                return tuple(sorted(set(leaving + (index,))))
    return ()


def _boron_transfer_coupling_like(observation: ReactionObservation) -> tuple[int, ...]:
    broken_cb = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "broken" and _element_pair(edit) == {"B", "C"},
    )
    formed_cc = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed" and _element_pair(edit) == {"C"},
    )
    broken_cx = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "broken"
        and "C" in _element_pair(edit)
        and bool(_element_pair(edit).intersection({"Cl", "Br", "I", "O"})),
    )
    return tuple(sorted(broken_cb + formed_cc + broken_cx)) if broken_cb and formed_cc and broken_cx else ()


def _heck_coupling_like(observation: ReactionObservation) -> tuple[int, ...]:
    substitution = _net_substitution(observation)
    formed_cc = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed" and _element_pair(edit) == {"C"},
    )
    h_loss = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.old_order is not None
        and edit.new_order is None,
    )
    return tuple(sorted(set(substitution + formed_cc + h_loss))) if substitution and formed_cc and h_loss else ()


def _cycloaddition_like(observation: ReactionObservation) -> tuple[int, ...]:
    ring = _net_ring_closure(observation)
    formed = _indices(observation.edits, lambda edit: edit.edit_type == "formed")
    order = _indices(observation.edits, lambda edit: edit.edit_type == "order_changed")
    return tuple(sorted(set(ring + formed + order))) if ring and len(formed) >= 2 and order else ()


def _is_carboxyl_carbon(
    observation: ReactionObservation,
    key: tuple[int, int],
) -> bool:
    component = next(
        (
            value
            for value in observation.reactants
            if value.component_index == key[0]
        ),
        None,
    )
    molecule = parse_smiles(component.input_smiles) if component is not None else None
    if molecule is None or not 0 <= key[1] < molecule.GetNumAtoms():
        return False
    atom = molecule.GetAtomWithIdx(key[1])
    if atom.GetSymbol() != "C":
        return False
    oxygen_orders = {
        str(molecule.GetBondBetweenAtoms(key[1], neighbor.GetIdx()).GetBondType()).upper()
        for neighbor in atom.GetNeighbors()
        if neighbor.GetSymbol() == "O"
    }
    return {"SINGLE", "DOUBLE"} <= oxygen_orders


def _decarboxylative_coupling_like(
    observation: ReactionObservation,
) -> tuple[int, ...]:
    """Match C–C formation coupled to loss of a carboxyl carbon."""
    formed_cc = _indices(
        observation.edits,
        lambda edit: edit.edit_type == "formed" and _element_pair(edit) == {"C"},
    )
    for broken_index, broken in enumerate(observation.edits):
        if broken.edit_type != "broken" or _element_pair(broken) != {"C"}:
            continue
        endpoint_keys = (_atom_key(broken, 1), _atom_key(broken, 2))
        for carboxyl_key, transfer_key in (
            endpoint_keys,
            endpoint_keys[::-1],
        ):
            if (
                carboxyl_key is None
                or transfer_key is None
                or not _is_carboxyl_carbon(observation, carboxyl_key)
            ):
                continue
            coupled = tuple(
                index
                for index in formed_cc
                if transfer_key
                in {
                    _atom_key(observation.edits[index], 1),
                    _atom_key(observation.edits[index], 2),
                }
            )
            if coupled:
                return tuple(sorted((broken_index, coupled[0])))
    return ()


_MATCHERS: Mapping[str, Callable[[ReactionObservation], tuple[int, ...]]] = {
    "net_substitution": _net_substitution,
    "net_elimination": _net_elimination,
    "net_addition": _net_addition,
    "net_bond_order_change": _net_bond_order_change,
    "net_coupling": _net_coupling,
    "net_ring_closure": _net_ring_closure,
    "reductive_amination_like": _reductive_amination_like,
    "amide_formation_like": _amide_formation_like,
    "boron_transfer_coupling_like": _boron_transfer_coupling_like,
    "heck_coupling_like": _heck_coupling_like,
    "cycloaddition_like": _cycloaddition_like,
    "decarboxylative_coupling_like": _decarboxylative_coupling_like,
}


def _confidence(observation: ReactionObservation) -> float:
    return max(0.0, min(1.0, float(observation.evidence_confidence)))


def match_reaction_patterns(
    observation: ReactionObservation,
) -> Tuple[ReactionPatternMatch, ...]:
    """Return all patterns supported by existing normalized evidence."""
    if not observation.valid or not observation.edits:
        return ()
    matches = []
    for definition in load_reaction_pattern_definitions():
        matcher_name = str(definition["matcher"])
        matcher = _MATCHERS.get(matcher_name)
        if matcher is None:
            raise ValueError(f"Unknown reaction-pattern matcher: {matcher_name}")
        indices = tuple(sorted(set(matcher(observation))))
        if not indices:
            continue
        matches.append(
            ReactionPatternMatch(
                pattern_id=str(definition["id"]),
                tier=str(definition["tier"]),  # type: ignore[arg-type]
                confidence=_confidence(observation),
                specificity=int(definition.get("specificity", 0)),
                display_importance=int(definition.get("display_importance", 0)),
                matched_edit_indices=indices,
                evidence=(observation.evidence_quality,),
                display_label=(
                    str(definition.get("display_label"))
                    if definition.get("display_label")
                    else None
                ),
                compatible_named_families=tuple(
                    definition.get("compatible_named_families") or ()
                ),
            )
        )
    return tuple(
        sorted(
            matches,
            key=lambda match: (
                -match.specificity,
                -match.display_importance,
                match.pattern_id,
            ),
        )
    )


def select_primary_pattern(
    matches: Iterable[ReactionPatternMatch],
) -> str | None:
    """Choose a display-only primary pattern from already supported matches."""
    values = tuple(matches)
    if not values:
        return None
    return min(
        values,
        key=lambda match: (
            -match.specificity,
            -match.display_importance,
            match.pattern_id,
        ),
    ).pattern_id


__all__ = [
    "REACTION_PATTERN_DEFINITION_SCHEMA_VERSION",
    "load_reaction_pattern_definitions",
    "match_reaction_patterns",
    "select_primary_pattern",
]
