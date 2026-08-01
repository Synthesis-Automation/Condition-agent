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
from .reaction_pattern_predicates import ReactionPatternContext, unique_indices


REACTION_PATTERN_DEFINITION_SCHEMA_VERSION = "2.0"
_DEFINITION_FILES = (
    "transformation_patterns.v2.json",
    "synthesis_patterns.v2.json",
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
        "matcher",
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
            predicate_id = str(pattern.get("predicate_id") or "")
            if not predicate_id:
                raise ValueError(f"Reaction pattern has no predicate: {pattern_id}")
            if predicate_id not in _PREDICATES:
                raise ValueError(
                    f"Unknown reaction-pattern predicate: {predicate_id}"
                )
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


def _organoboron_c_c_coupling_like(
    observation: ReactionObservation,
) -> tuple[int, ...]:
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


_LEAVING_ELEMENTS = frozenset({"F", "Cl", "Br", "I", "O", "S"})


def _centered_replacement(
    context: ReactionPatternContext,
    *,
    installed_element: str,
    carbon_kind: str,
) -> tuple[int, ...]:
    """Match C–leaving-group replacement at one carbon center."""
    for formed_index, formed in enumerate(context.edits):
        if (
            formed.edit_type != "formed"
            or context.element_pair(formed) != {"C", installed_element}
        ):
            continue
        carbon = context.endpoint_for_element(formed, "C")
        if carbon is None:
            continue
        _, carbon_reference, center_key = carbon
        if context.is_carbonyl_carbon(center_key):
            continue
        if carbon_kind == "sp2" and not context.is_sp2_carbon(carbon_reference):
            continue
        if carbon_kind == "sp3" and not context.is_sp3_carbon(carbon_reference):
            continue
        leaving = []
        for broken_index in context.edits_at(center_key, edit_type="broken"):
            broken = context.edits[broken_index]
            other_elements = context.element_pair(broken) - {"C"}
            if other_elements and other_elements <= _LEAVING_ELEMENTS:
                leaving.append(broken_index)
        if leaving:
            hydrogen = context.indices(
                lambda edit: edit.edit_type == "hydrogen_change"
                and edit.atom_1.element == installed_element
            )
            return unique_indices(leaving, (formed_index,), hydrogen)
    return ()


def _sp2_c_n_substitution_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _centered_replacement(
        context, installed_element="N", carbon_kind="sp2"
    )


def _sp3_c_n_substitution_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _centered_replacement(
        context, installed_element="N", carbon_kind="sp3"
    )


def _sp2_c_o_substitution_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _centered_replacement(
        context, installed_element="O", carbon_kind="sp2"
    )


def _sp3_c_o_substitution_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _centered_replacement(
        context, installed_element="O", carbon_kind="sp3"
    )


def _sulfonyl_transfer_like(
    context: ReactionPatternContext,
    installed_element: str,
) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if (
            formed.edit_type != "formed"
            or context.element_pair(formed) != {"S", installed_element}
        ):
            continue
        sulfur = context.endpoint_for_element(formed, "S")
        if sulfur is None or not context.is_sulfonyl_sulfur(sulfur[2]):
            continue
        leaving = tuple(
            index
            for index in context.edits_at(sulfur[2], edit_type="broken")
            if bool(context.element_pair(context.edits[index]) & {"F", "Cl", "Br", "I"})
        )
        if leaving:
            return unique_indices(leaving, (formed_index,))
    return ()


def _sulfonamide_formation_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    return _sulfonyl_transfer_like(context, "N")


def _o_sulfonylation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _sulfonyl_transfer_like(context, "O")


def _acyl_substitution(
    context: ReactionPatternContext,
    installed_element: str,
) -> tuple[tuple[int, ...], tuple[int, int] | None, tuple[str, ...]]:
    """Return acyl-substitution edits, center, and retained single neighbors."""
    for formed_index, formed in enumerate(context.edits):
        if (
            formed.edit_type != "formed"
            or context.element_pair(formed) != {"C", installed_element}
        ):
            continue
        carbon = context.endpoint_for_element(formed, "C")
        if carbon is None or not context.is_carbonyl_carbon(carbon[2]):
            continue
        leaving_indices = []
        leaving_atom_indices = []
        for broken_index in context.edits_at(carbon[2], edit_type="broken"):
            broken = context.edits[broken_index]
            if str(broken.old_order or "").upper() != "SINGLE":
                continue
            other = (
                broken.atom_2
                if context.atom_key(broken, 1) == carbon[2]
                else broken.atom_1
            )
            if other is not None:
                leaving_indices.append(broken_index)
                leaving_atom_indices.append(other.atom_index)
        if not leaving_indices:
            continue
        molecule = context.molecule("reactant", carbon[2][0])
        retained = []
        if molecule is not None:
            atom = molecule.GetAtomWithIdx(carbon[2][1])
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in leaving_atom_indices:
                    continue
                bond = molecule.GetBondBetweenAtoms(carbon[2][1], neighbor.GetIdx())
                if str(bond.GetBondType()).upper() == "SINGLE":
                    retained.append(str(neighbor.GetSymbol()))
        return (
            unique_indices(tuple(leaving_indices), (formed_index,)),
            carbon[2],
            tuple(sorted(retained)),
        )
    return (), None, ()


def _amide_formation(context: ReactionPatternContext) -> tuple[int, ...]:
    indices, _, retained = _acyl_substitution(context, "N")
    return indices if indices and not ({"O", "N"} & set(retained)) else ()


def _ester_formation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    indices, _, retained = _acyl_substitution(context, "O")
    return indices if indices and "O" not in retained else ()


def _carbamate_formation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    indices, _, retained = _acyl_substitution(context, "N")
    return indices if indices and "O" in retained else ()


def _carbonate_formation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    indices, _, retained = _acyl_substitution(context, "O")
    return indices if indices and "O" in retained else ()


def _urea_formation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    formed_cn = context.indices(
        lambda edit: edit.edit_type == "formed"
        and context.element_pair(edit) == {"C", "N"}
    )
    for formed_index in formed_cn:
        carbon = context.endpoint_for_element(context.edits[formed_index], "C")
        if carbon is None or not context.is_carbonyl_carbon(carbon[2]):
            continue
        if not context.has_neighbor(carbon[2], element="N", order="DOUBLE"):
            continue
        order_changes = tuple(
            index
            for index in context.edits_at(carbon[2], edit_type="order_changed")
            if context.element_pair(context.edits[index]) == {"C", "N"}
            and context.order_direction(context.edits[index]) < 0
        )
        if order_changes:
            return unique_indices((formed_index,), order_changes)
    return ()


def _organoboron_transfer_like(
    context: ReactionPatternContext,
    installed_element: str,
) -> tuple[int, ...]:
    broken_cb = context.indices(
        lambda edit: edit.edit_type == "broken"
        and context.element_pair(edit) == {"B", "C"}
    )
    for broken_index in broken_cb:
        carbon = context.endpoint_for_element(context.edits[broken_index], "C")
        if carbon is None:
            continue
        formed = tuple(
            index
            for index in context.edits_at(carbon[2], edit_type="formed")
            if context.element_pair(context.edits[index]) == {"C", installed_element}
        )
        if formed:
            return unique_indices((broken_index,), formed)
    return ()


def _organoboron_c_n_coupling_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    return _organoboron_transfer_like(context, "N")


def _organoboron_c_o_coupling_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    return _organoboron_transfer_like(context, "O")


def _organoboron_c_s_coupling_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    return _organoboron_transfer_like(context, "S")


def _sonogashira_coupling_like(context: ReactionPatternContext) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if formed.edit_type != "formed" or context.element_pair(formed) != {"C"}:
            continue
        endpoint_keys = (context.atom_key(formed, 1), context.atom_key(formed, 2))
        for alkynyl_key, electrophile_key in (endpoint_keys, endpoint_keys[::-1]):
            if alkynyl_key is None or electrophile_key is None:
                continue
            if not context.has_neighbor(alkynyl_key, element="C", order="TRIPLE"):
                continue
            leaving = tuple(
                index
                for index in context.edits_at(electrophile_key, edit_type="broken")
                if bool(context.element_pair(context.edits[index]) & {"F", "Cl", "Br", "I"})
            )
            if leaving:
                return unique_indices(leaving, (formed_index,))
    return ()


def _borylation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if formed.edit_type != "formed" or context.element_pair(formed) != {"B", "C"}:
            continue
        carbon = context.endpoint_for_element(formed, "C")
        if carbon is None:
            continue
        leaving = tuple(
            index
            for index in context.edits_at(carbon[2], edit_type="broken")
            if bool(context.element_pair(context.edits[index]) & {"F", "Cl", "Br", "I"})
        )
        if leaving:
            return unique_indices(leaving, (formed_index,))
    return ()


def _alkene_hydrogenation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    order = context.indices(
        lambda edit: edit.edit_type == "order_changed"
        and context.element_pair(edit) == {"C"}
        and str(edit.old_order or "").upper() == "DOUBLE"
        and context.order_direction(edit) < 0
    )
    hydrogen = context.indices(
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "C"
        and edit.old_order is None
        and edit.new_order is not None
    )
    return unique_indices(order, hydrogen) if order and hydrogen else ()


def _alkyne_hydrogenation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    order = context.indices(
        lambda edit: edit.edit_type == "order_changed"
        and context.element_pair(edit) == {"C"}
        and str(edit.old_order or "").upper() == "TRIPLE"
        and context.order_direction(edit) < 0
    )
    hydrogen = context.indices(
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "C"
        and edit.old_order is None
        and edit.new_order is not None
    )
    return unique_indices(order, hydrogen) if order and hydrogen else ()


def _carbonyl_reduction_like(context: ReactionPatternContext) -> tuple[int, ...]:
    for index, edit in enumerate(context.edits):
        if not (
            edit.edit_type == "order_changed"
            and context.element_pair(edit) == {"C", "O"}
            and context.order_direction(edit) < 0
        ):
            continue
        carbon = context.endpoint_for_element(edit, "C")
        if carbon is None:
            continue
        if any(
            context.edits[related_index].edit_type in {"formed", "broken"}
            for related_index in context.edits_at(carbon[2])
        ):
            continue
        if not (
            context.has_neighbor(carbon[2], element="O", order="SINGLE")
            or context.has_neighbor(carbon[2], element="N", order="SINGLE")
        ):
            return (index,)
    return ()


def _alcohol_oxidation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    order = context.indices(
        lambda edit: edit.edit_type == "order_changed"
        and context.element_pair(edit) == {"C", "O"}
        and context.order_direction(edit) > 0
    )
    return order


def _nitro_reduction_like(context: ReactionPatternContext) -> tuple[int, ...]:
    broken_no = context.indices(
        lambda edit: edit.edit_type == "broken"
        and context.element_pair(edit) == {"N", "O"}
    )
    hydrogen_n = context.indices(
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "N"
        and edit.old_order is None
        and edit.new_order is not None
    )
    return unique_indices(broken_no, hydrogen_n) if broken_no and hydrogen_n else ()


def _heteroatom_deprotection_like(
    context: ReactionPatternContext,
    element: str,
) -> tuple[int, ...]:
    hydrogen = context.indices(
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == element
        and edit.old_order is None
        and edit.new_order is not None
    )
    for hydrogen_index in hydrogen:
        key = context.atom_key(context.edits[hydrogen_index], 1)
        if key is None:
            continue
        if element == "O":
            molecule = context.molecule("reactant", key[0])
            if molecule is not None:
                oxygen = molecule.GetAtomWithIdx(key[1])
                if any(
                    neighbor.GetSymbol() == "C"
                    and context.is_carbonyl_carbon((key[0], neighbor.GetIdx()))
                    for neighbor in oxygen.GetNeighbors()
                ):
                    continue
        broken = tuple(
            index
            for index in context.edits_at(key, edit_type="broken")
            if bool(context.element_pair(context.edits[index]) & {"C", "Si", "S"})
        )
        formed = context.edits_at(key, edit_type="formed")
        if broken and not formed:
            return unique_indices(broken, (hydrogen_index,))
    return ()


def _amine_deprotection_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _heteroatom_deprotection_like(context, "N")


def _alcohol_deprotection_like(context: ReactionPatternContext) -> tuple[int, ...]:
    return _heteroatom_deprotection_like(context, "O")


def _ester_hydrolysis_like(context: ReactionPatternContext) -> tuple[int, ...]:
    """Match ester cleavage that restores a carboxylic-acid O–H bond."""
    hydrogen_o = context.indices(
        lambda edit: edit.edit_type == "hydrogen_change"
        and edit.atom_1.element == "O"
        and edit.old_order is None
        and edit.new_order is not None
    )
    for hydrogen_index in hydrogen_o:
        oxygen_key = context.atom_key(context.edits[hydrogen_index], 1)
        if oxygen_key is None:
            continue
        molecule = context.molecule("reactant", oxygen_key[0])
        if molecule is None:
            continue
        oxygen = molecule.GetAtomWithIdx(oxygen_key[1])
        carbonyl_neighbors = tuple(
            neighbor.GetIdx()
            for neighbor in oxygen.GetNeighbors()
            if neighbor.GetSymbol() == "C"
            and context.is_carbonyl_carbon((oxygen_key[0], neighbor.GetIdx()))
        )
        if not carbonyl_neighbors:
            continue
        broken = context.edits_at(oxygen_key, edit_type="broken")
        if broken:
            return unique_indices(broken, (hydrogen_index,))
    return ()


def _acid_chloride_formation_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if formed.edit_type != "formed" or context.element_pair(formed) != {"C", "Cl"}:
            continue
        carbon = context.endpoint_for_element(formed, "C")
        if carbon is None or not context.is_carbonyl_carbon(carbon[2]):
            continue
        broken = tuple(
            index
            for index in context.edits_at(carbon[2], edit_type="broken")
            if "O" in context.element_pair(context.edits[index])
        )
        if broken:
            return unique_indices(broken, (formed_index,))
    return ()


def _alcohol_to_halide_like(context: ReactionPatternContext) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if formed.edit_type != "formed" or not (
            "C" in context.element_pair(formed)
            and context.element_pair(formed) & {"F", "Cl", "Br", "I"}
        ):
            continue
        carbon = context.endpoint_for_element(formed, "C")
        if carbon is None or context.is_carbonyl_carbon(carbon[2]):
            continue
        broken = tuple(
            index
            for index in context.edits_at(carbon[2], edit_type="broken")
            if "O" in context.element_pair(context.edits[index])
        )
        if broken:
            return unique_indices(broken, (formed_index,))
    return ()


def _cyanation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    for formed_index, formed in enumerate(context.edits):
        if formed.edit_type != "formed" or context.element_pair(formed) != {"C"}:
            continue
        keys = (context.atom_key(formed, 1), context.atom_key(formed, 2))
        for nitrile_key, electrophile_key in (keys, keys[::-1]):
            if nitrile_key is None or electrophile_key is None:
                continue
            if not context.has_neighbor(nitrile_key, element="N", order="TRIPLE"):
                continue
            leaving = tuple(
                index
                for index in context.edits_at(electrophile_key, edit_type="broken")
                if bool(context.element_pair(context.edits[index]) & _LEAVING_ELEMENTS)
            )
            if leaving:
                return unique_indices(leaving, (formed_index,))
    return ()


def _sulfide_oxidation_like(context: ReactionPatternContext) -> tuple[int, ...]:
    formed_so = context.indices(
        lambda edit: edit.edit_type == "formed"
        and context.element_pair(edit) == {"O", "S"}
        and str(edit.new_order or "").upper() in {"SINGLE", "DOUBLE"}
    )
    for index in formed_so:
        sulfur = context.endpoint_for_element(context.edits[index], "S")
        if sulfur is not None and not context.is_sulfonyl_sulfur(sulfur[2]):
            return (index,)
    return ()


def _carboxyl_derivative_reduction_like(
    context: ReactionPatternContext,
) -> tuple[int, ...]:
    order = context.indices(
        lambda edit: edit.edit_type == "order_changed"
        and context.element_pair(edit) == {"C", "O"}
        and context.order_direction(edit) < 0
    )
    for index in order:
        carbon = context.endpoint_for_element(context.edits[index], "C")
        if carbon is not None and context.has_neighbor(
            carbon[2], element="O", order="SINGLE"
        ):
            related = context.edits_at(carbon[2])
            return unique_indices((index,), related)
    return ()


def _net_bond_cleavage(observation: ReactionObservation) -> tuple[int, ...]:
    broken = _indices(observation.edits, lambda edit: edit.edit_type == "broken")
    formed = _indices(observation.edits, lambda edit: edit.edit_type == "formed")
    return broken if broken and not formed else ()


def _net_ring_opening(observation: ReactionObservation) -> tuple[int, ...]:
    topology = observation.topology
    if topology is None or not (
        topology.ring_count_delta is not None and topology.ring_count_delta < 0
    ):
        return ()
    return _indices(observation.edits, lambda edit: edit.edit_type == "broken")


def _net_bond_order_direction(
    observation: ReactionObservation,
    direction: int,
) -> tuple[int, ...]:
    return _indices(
        observation.edits,
        lambda edit: edit.edit_type == "order_changed"
        and _order_direction(edit) == direction,
    )


_PREDICATES: Mapping[str, Callable[[ReactionPatternContext], tuple[int, ...]]] = {
    "net_substitution": lambda context: _net_substitution(context.observation),
    "net_elimination": lambda context: _net_elimination(context.observation),
    "net_addition": lambda context: _net_addition(context.observation),
    "net_bond_order_change": lambda context: _net_bond_order_change(
        context.observation
    ),
    "net_coupling": lambda context: _net_coupling(context.observation),
    "net_ring_closure": lambda context: _net_ring_closure(context.observation),
    "net_bond_cleavage": lambda context: _net_bond_cleavage(
        context.observation
    ),
    "net_ring_opening": lambda context: _net_ring_opening(context.observation),
    "net_bond_order_increase": lambda context: _net_bond_order_direction(
        context.observation, 1
    ),
    "net_bond_order_decrease": lambda context: _net_bond_order_direction(
        context.observation, -1
    ),
    "reductive_amination_like": lambda context: _reductive_amination_like(
        context.observation
    ),
    "amide_formation_like": _amide_formation,
    "ester_formation_like": _ester_formation_like,
    "carbamate_formation_like": _carbamate_formation_like,
    "carbonate_formation_like": _carbonate_formation_like,
    "urea_formation_like": _urea_formation_like,
    "sp2_c_n_substitution_like": _sp2_c_n_substitution_like,
    "sp3_c_n_substitution_like": _sp3_c_n_substitution_like,
    "sp2_c_o_substitution_like": _sp2_c_o_substitution_like,
    "sp3_c_o_substitution_like": _sp3_c_o_substitution_like,
    "sulfonamide_formation_like": _sulfonamide_formation_like,
    "o_sulfonylation_like": _o_sulfonylation_like,
    "organoboron_c_c_coupling_like": lambda context: _organoboron_c_c_coupling_like(
        context.observation
    ),
    "organoboron_c_n_coupling_like": _organoboron_c_n_coupling_like,
    "organoboron_c_o_coupling_like": _organoboron_c_o_coupling_like,
    "organoboron_c_s_coupling_like": _organoboron_c_s_coupling_like,
    "sonogashira_coupling_like": _sonogashira_coupling_like,
    "borylation_like": _borylation_like,
    "alkene_hydrogenation_like": _alkene_hydrogenation_like,
    "alkyne_hydrogenation_like": _alkyne_hydrogenation_like,
    "carbonyl_reduction_like": _carbonyl_reduction_like,
    "alcohol_oxidation_like": _alcohol_oxidation_like,
    "nitro_reduction_like": _nitro_reduction_like,
    "amine_deprotection_like": _amine_deprotection_like,
    "alcohol_deprotection_like": _alcohol_deprotection_like,
    "ester_hydrolysis_like": _ester_hydrolysis_like,
    "acid_chloride_formation_like": _acid_chloride_formation_like,
    "alcohol_to_halide_like": _alcohol_to_halide_like,
    "cyanation_like": _cyanation_like,
    "sulfide_oxidation_like": _sulfide_oxidation_like,
    "carboxyl_derivative_reduction_like": _carboxyl_derivative_reduction_like,
    "heck_coupling_like": lambda context: _heck_coupling_like(
        context.observation
    ),
    "cycloaddition_like": lambda context: _cycloaddition_like(
        context.observation
    ),
    "decarboxylative_coupling_like": lambda context: _decarboxylative_coupling_like(
        context.observation
    ),
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
    context = ReactionPatternContext(observation)
    for definition in load_reaction_pattern_definitions():
        predicate_id = str(definition["predicate_id"])
        predicate = _PREDICATES.get(predicate_id)
        if predicate is None:
            raise ValueError(f"Unknown reaction-pattern predicate: {predicate_id}")
        indices = tuple(sorted(set(predicate(context))))
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
                requires_condition_evidence=bool(
                    definition.get("requires_condition_evidence", False)
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
