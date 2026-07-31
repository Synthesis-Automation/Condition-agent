"""Explicit atom-property and stereochemical changes in reaction cores."""

from __future__ import annotations

from typing import Any, Iterable, Tuple

from ..reaction_models import ReactionStereoChange
from .keys import digest
from .models import ReactionCoreAtomTransition, ReactionCoreStateChange


def _text(value: Any) -> str:
    if isinstance(value, bool):
        return "aromatic" if value else "aliphatic"
    return str(value)


def _state_change(
    *,
    change_type: str,
    transition: ReactionCoreAtomTransition,
    before_value: Any,
    after_value: Any,
    evidence: str,
) -> ReactionCoreStateChange:
    before = transition.before_state
    after = transition.after_state
    elements = tuple(
        dict.fromkeys(
            state.element for state in (before, after) if state is not None
        )
    )
    payload = (
        change_type,
        elements,
        _text(before_value),
        _text(after_value),
    )
    return ReactionCoreStateChange(
        change_id=digest("RSC1", payload, length=24),
        change_type=change_type,  # type: ignore[arg-type]
        atom_map_numbers=(
            (int(transition.atom_map_number),)
            if transition.atom_map_number is not None
            else ()
        ),
        elements=elements,
        before_value=_text(before_value),
        after_value=_text(after_value),
        evidence=evidence,
    )


def _transition_changes(
    transition: ReactionCoreAtomTransition,
    *,
    evidence: str,
) -> Iterable[ReactionCoreStateChange]:
    before = transition.before_state
    after = transition.after_state
    if before is None or after is None:
        return ()
    fields = (
        ("hydrogen", before.total_hydrogens, after.total_hydrogens),
        ("formal_charge", before.formal_charge, after.formal_charge),
        ("radical", before.radical_electrons, after.radical_electrons),
        ("isotope", before.isotope, after.isotope),
        ("aromaticity", before.aromatic, after.aromatic),
        ("hybridization", before.hybridization, after.hybridization),
    )
    return tuple(
        _state_change(
            change_type=change_type,
            transition=transition,
            before_value=before_value,
            after_value=after_value,
            evidence=evidence,
        )
        for change_type, before_value, after_value in fields
        if before_value != after_value
    )


def _stereo_change(
    change: ReactionStereoChange,
) -> ReactionCoreStateChange:
    references = tuple(
        reference
        for reference in (change.atom_1, change.atom_2)
        if reference is not None
    )
    change_type = (
        "atom_stereochemistry"
        if change.stereo_type == "atom"
        else "bond_stereochemistry"
    )
    elements = tuple(reference.element for reference in references)
    before_value = str(change.old_descriptor or "unspecified")
    after_value = str(change.new_descriptor or "unspecified")
    payload = (change_type, elements, before_value, after_value)
    return ReactionCoreStateChange(
        change_id=digest("RSC1", payload, length=24),
        change_type=change_type,  # type: ignore[arg-type]
        atom_map_numbers=tuple(
            sorted(
                int(reference.atom_map_number)
                for reference in references
                if reference.atom_map_number is not None
            )
        ),
        elements=elements,
        before_value=before_value,
        after_value=after_value,
        evidence=change.evidence,
    )


def build_state_changes(
    transitions: Iterable[ReactionCoreAtomTransition],
    stereo_changes: Iterable[ReactionStereoChange],
    *,
    evidence: str,
) -> Tuple[ReactionCoreStateChange, ...]:
    """Return deterministic, explicit property changes for an observed core."""
    values = [
        value
        for transition in transitions
        for value in _transition_changes(transition, evidence=evidence)
    ]
    values.extend(_stereo_change(change) for change in stereo_changes)
    return tuple(
        sorted(
            values,
            key=lambda value: (
                value.change_type,
                value.elements,
                value.before_value,
                value.after_value,
                value.atom_map_numbers,
                value.change_id,
            ),
        )
    )


__all__ = ["build_state_changes"]
