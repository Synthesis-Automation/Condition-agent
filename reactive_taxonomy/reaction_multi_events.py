"""Exact reconstruction of balanced, unmapped multi-event reactions."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Dict, Iterable, Sequence, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .connectivity_rewrite import connectivity_rewrite_by_id
from .reaction_models import ReactionComponent, ReactionSiteReference
from .reaction_site_interfaces import normalize_reaction_assignment


RawCandidate = Tuple[Dict[str, Any], Dict[str, ReactionSiteReference]]
_Coordinate = Tuple[int, int]


@dataclass(frozen=True)
class _EventProgram:
    """One composable release-and-connect program in reactant coordinates."""

    participants: Tuple[ReactionSiteReference, ...]
    join: Tuple[_Coordinate, _Coordinate]
    removals: Tuple[Tuple[int, int, int], ...]


def _instruction_atom(
    selector: object,
    normalized: Dict[str, Any],
) -> _Coordinate | None:
    """Resolve the instruction selectors needed by multi-event composition."""
    parts = str(selector or "").split(".")
    if len(parts) != 3 or parts[0] not in normalized:
        return None
    role, interface, attribute = parts
    view = normalized[role]
    if interface == "reactive_link" and len(view.reactive_links) == 1:
        link = view.reactive_links[0]
        endpoint = getattr(link, attribute, None)
    elif (
        interface == "connection_endpoint"
        and attribute == "atom"
        and len(view.connection_endpoints) == 1
    ):
        endpoint = view.connection_endpoints[0].endpoint
    else:
        return None
    atom_index = getattr(endpoint, "atom_index", None)
    component_index = getattr(endpoint, "component_index", None)
    if atom_index is None or component_index is None:
        return None
    return int(component_index), int(atom_index)


def _instruction_event_program(
    candidate: RawCandidate,
    reactants: Tuple[ReactionComponent, ...],
) -> _EventProgram | None:
    """Compile one declarative release-and-connect variant for composition."""
    rule, assignment = candidate
    rewrite = connectivity_rewrite_by_id(str(rule.get("operator_id") or ""))
    if rewrite is None or len(rewrite.variants) != 1:
        return None
    bindings = rule.get("operator_slot_bindings") or {}
    execution_assignment = {
        str(operator_slot): assignment[str(rule_slot)]
        for operator_slot, rule_slot in bindings.items()
    }
    try:
        normalized = normalize_reaction_assignment(execution_assignment, reactants)
    except (KeyError, StopIteration, ValueError):
        return None
    joins = []
    removals = []
    for instruction in rewrite.variants[0].instructions:
        operation = str(instruction.get("op") or "")
        if operation == "change_localized_bond_state" and (
            str(instruction.get("before") or "").upper() == "NONE"
            and str(instruction.get("after") or "").upper() == "SINGLE"
        ):
            endpoints = tuple(instruction.get("endpoints") or ())
            if len(endpoints) != 2:
                return None
            left = _instruction_atom(endpoints[0], normalized)
            right = _instruction_atom(endpoints[1], normalized)
            if left is None or right is None:
                return None
            joins.append((left, right))
        elif operation == "declare_projection_discardable_attachment":
            retained = _instruction_atom(instruction.get("retained"), normalized)
            discarded = _instruction_atom(
                instruction.get("discarded"), normalized
            )
            if (
                retained is None
                or discarded is None
                or retained[0] != discarded[0]
            ):
                return None
            removals.append((retained[0], retained[1], discarded[1]))
    if len(joins) != 1 or not removals:
        return None
    participants = tuple(
        sorted(
            assignment.values(),
            key=lambda site: (
                site.component_index,
                site.site_id,
                site.canonical_signature,
            ),
        )
    )
    return _EventProgram(
        participants=participants,
        join=joins[0],
        removals=tuple(sorted(set(removals))),
    )


def _event_sites(
    candidate: RawCandidate,
    reactants: Tuple[ReactionComponent, ...],
) -> _EventProgram | None:
    return _instruction_event_program(candidate, reactants)


def _operation_key(candidate: RawCandidate) -> Tuple[Any, ...]:
    rule, assignment = candidate
    return (
        rule["id"],
        tuple(
            sorted(
                (
                    site.component_index,
                    site.site_id,
                    site.canonical_signature,
                )
                for site in assignment.values()
            )
        ),
    )


def _component_by_index(
    components: Tuple[ReactionComponent, ...],
    component_index: int,
) -> ReactionComponent:
    return next(
        component
        for component in components
        if component.component_index == component_index
    )


def _fragment_to_remove(
    components: Tuple[ReactionComponent, ...],
    component_index: int,
    retained_atom: int,
    discarded_atom: int,
) -> set[int]:
    from rdkit import Chem

    molecule = parse_smiles(
        _component_by_index(components, component_index).input_smiles
    )
    if molecule is None:
        return set()
    rw = Chem.RWMol(molecule)
    rw.RemoveBond(retained_atom, discarded_atom)
    return next(
        (
            set(int(index) for index in fragment)
            for fragment in Chem.GetMolFrags(
                rw.GetMol(),
                asMols=False,
                sanitizeFrags=False,
            )
            if discarded_atom in fragment
        ),
        set(),
    )


def apply_rewrite_sequence(
    operations: Sequence[RawCandidate],
    components: Tuple[ReactionComponent, ...],
) -> str | None:
    """Apply compatible release-and-connect programs as one edit graph."""
    from rdkit import Chem

    if len(operations) < 2:
        return None
    resolved = [
        _event_sites(candidate, components) for candidate in operations
    ]
    if any(event is None for event in resolved):
        return None
    events = [event for event in resolved if event is not None]
    participant_indices: set[int] = set()
    removals: Dict[int, set[int]] = {}
    joins: list[Tuple[Tuple[int, int], Tuple[int, int]]] = []
    used_join_atoms: set[_Coordinate] = set()
    for event in events:
        left, right = event.join
        if left in used_join_atoms or right in used_join_atoms:
            return None
        used_join_atoms.update((left, right))
        for component_index, retained, discarded in event.removals:
            fragment = _fragment_to_remove(
                components,
                component_index,
                retained,
                discarded,
            )
            if not fragment:
                return None
            removals.setdefault(component_index, set()).update(fragment)
        joins.append(event.join)
        participant_indices.update(
            site.component_index for site in event.participants
        )

    used_indices = sorted(participant_indices)
    molecules = []
    offsets: Dict[int, int] = {}
    total = 0
    for component_index in used_indices:
        molecule = parse_smiles(
            _component_by_index(
                components, component_index
            ).input_smiles
        )
        if molecule is None:
            return None
        offsets[component_index] = total
        total += molecule.GetNumAtoms()
        molecules.append(molecule)
    combined = molecules[0]
    for molecule in molecules[1:]:
        combined = Chem.CombineMols(combined, molecule)
    removed_global = {
        offsets[component_index] + atom_index
        for component_index, atom_indices in removals.items()
        for atom_index in atom_indices
    }
    global_joins = [
        (
            offsets[left_component] + left_atom,
            offsets[right_component] + right_atom,
        )
        for (left_component, left_atom), (
            right_component,
            right_atom,
        ) in joins
    ]
    if any(
        left in removed_global or right in removed_global
        for left, right in global_joins
    ):
        return None
    rw = Chem.RWMol(combined)
    for atom_index in sorted(removed_global, reverse=True):
        rw.RemoveAtom(atom_index)

    removed_ascending = tuple(sorted(removed_global))

    def shifted(index: int) -> int:
        return index - sum(
            removed < index for removed in removed_ascending
        )

    for left_global, right_global in global_joins:
        left, right = shifted(left_global), shifted(right_global)
        if (
            left == right
            or rw.GetBondBetweenAtoms(left, right) is not None
        ):
            return None
        rw.AddBond(left, right, Chem.BondType.SINGLE)
    product = rw.GetMol()
    try:
        product.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(product)
        return Chem.MolToSmiles(
            product, canonical=True, isomericSmiles=True
        )
    except Exception:
        return None


def _interpretation_key(candidates: Sequence[RawCandidate]) -> Tuple[Any, ...]:
    return tuple(
        sorted(
            (
                rule["id"],
                tuple(
                    sorted(
                        site.canonical_signature for site in assignment.values()
                    )
                ),
            )
            for rule, assignment in candidates
        )
    )


def exact_multi_event_reconstructions(
    raw_candidates: Iterable[RawCandidate],
    reactants: Tuple[ReactionComponent, ...],
    observed_products: set[str],
    *,
    max_events: int = 4,
    max_combinations: int = 5000,
) -> Tuple[Tuple[RawCandidate, ...], ...]:
    """Return composite rewrite sets that exactly reconstruct one product."""
    eligible = tuple(
        sorted(
            (
                candidate
                for candidate in raw_candidates
                if _event_sites(candidate, reactants) is not None
            ),
            key=_operation_key,
        )
    )
    if len(eligible) < 2 or not observed_products:
        return ()
    exact = []
    attempted = 0
    upper = min(max_events, len(eligible))
    for event_count in range(2, upper + 1):
        for selected in combinations(eligible, event_count):
            attempted += 1
            if attempted > max_combinations:
                return tuple(exact)
            predicted = apply_rewrite_sequence(selected, reactants)
            if predicted in observed_products:
                exact.append(tuple(selected))
        if exact:
            break
    if not exact:
        return ()
    return tuple(
        sorted(
            exact,
            key=lambda selected: (
                _interpretation_key(selected),
                tuple(_operation_key(candidate) for candidate in selected),
            ),
        )
    )


def equivalent_multi_event_interpretations(
    reconstructions: Sequence[Sequence[RawCandidate]],
) -> bool:
    """Return whether exact alternatives differ only by equivalent sites."""
    return len({_interpretation_key(selected) for selected in reconstructions}) <= 1


__all__ = [
    "RawCandidate",
    "apply_rewrite_sequence",
    "equivalent_multi_event_interpretations",
    "exact_multi_event_reconstructions",
]
