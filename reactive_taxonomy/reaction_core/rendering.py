"""Chemist-facing labels for graph-derived reaction cores.

These functions render already-observed chemistry.  Their output is never used
to construct reaction-core identities or decide recommendation eligibility.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping, Optional, Sequence

from ..reaction_models import (
    ReactionAtomReference,
    ReactionCoreAtomTransition,
    ReactionCorePresentation,
    ReactionCoreQuality,
    ReactionCoreRemoteClass,
    ReactionCoreRemoteSubgraph,
    ReactionCoreStateChange,
    ReactionEdit,
)
from .common import AtomIdentity, Coordinate, atom_identity


def _remote_display(remote_class: ReactionCoreRemoteClass) -> str:
    return {
        "aryl": "Ar",
        "heteroaryl": "HetAr",
        "alkyl": "R",
        "alkenyl": "Alkenyl",
        "alkynyl": "Alkynyl",
        "acyl": "Acyl",
        "ring_aliphatic": "Cycloalkyl",
        "heteroatom": "X",
        "generic_R": "R",
    }[remote_class]


def _bond_prefix(order: str) -> str:
    return {
        "SINGLE": "",
        "DOUBLE": "=",
        "TRIPLE": "#",
        "AROMATIC": ":",
    }.get(str(order).upper(), "~")


def _edit_atom_display(atom: ReactionAtomReference) -> str:
    """Render an edit endpoint without duplicating its aromatic scaffold."""
    if atom.aromatic:
        return "Ar" if atom.element == "C" else "HetAr"
    return atom.element


def _edit_bond_display(
    atom_1: ReactionAtomReference,
    atom_2: Optional[ReactionAtomReference],
    order: Optional[str],
) -> str:
    """Render one schema-level bond state for an event equation."""
    left = _edit_atom_display(atom_1)
    right = "H" if atom_2 is None else _edit_atom_display(atom_2)
    bond = {
        "SINGLE": "–",
        "DOUBLE": "=",
        "TRIPLE": "≡",
        "AROMATIC": ":",
    }.get(str(order or "SINGLE").upper(), "~")
    endpoints = sorted((left, right))
    return f"{endpoints[0]}{bond}{endpoints[1]}"


def _counted_edit_terms(values: Sequence[str]) -> str:
    counts = Counter(values)
    return " + ".join(
        f"{count} × {term}" if count > 1 else term
        for term, count in sorted(counts.items())
    )


def multi_center_edit_label(edits: Sequence[ReactionEdit]) -> str:
    """Render a multi-center event once instead of once per center atom."""
    before = []
    after = []
    heavy_edits = tuple(edit for edit in edits if edit.atom_2 is not None)
    suppress_hydrogen = bool(heavy_edits) and all(
        edit.edit_type == "order_changed" for edit in heavy_edits
    )
    for edit in edits:
        if edit.edit_type == "broken":
            before.append(
                _edit_bond_display(edit.atom_1, edit.atom_2, edit.old_order)
            )
        elif edit.edit_type == "formed":
            after.append(
                _edit_bond_display(edit.atom_1, edit.atom_2, edit.new_order)
            )
        elif edit.edit_type == "order_changed":
            before.append(
                _edit_bond_display(edit.atom_1, edit.atom_2, edit.old_order)
            )
            after.append(
                _edit_bond_display(edit.atom_1, edit.atom_2, edit.new_order)
            )
        elif edit.edit_type == "hydrogen_change" and not suppress_hydrogen:
            term = _edit_bond_display(edit.atom_1, None, "SINGLE")
            if edit.old_order is not None:
                before.append(term)
            if edit.new_order is not None:
                after.append(term)
    return (
        f"{_counted_edit_terms(before) or '∅'}"
        " → "
        f"{_counted_edit_terms(after) or '∅'}"
    )


def single_center_transition_label(
    *,
    primary_identity: AtomIdentity,
    transition_by_identity: Mapping[
        AtomIdentity,
        ReactionCoreAtomTransition,
    ],
    edits: Sequence[ReactionEdit],
) -> str:
    """Render one center transition with external formed-bond partners."""
    primary = transition_by_identity[primary_identity]
    before_label = (
        primary.before_state.concise_label
        if primary.before_state is not None
        else "∅"
    )
    after_label = (
        primary.after_state.concise_label
        if primary.after_state is not None
        else "∅"
    )
    if primary.before_state is None:
        return f"{before_label} → {after_label}"

    partner_labels = set()
    primary_component = primary.before_state.component_index
    for edit in edits:
        if edit.edit_type != "formed" or edit.atom_2 is None:
            continue
        identities = (atom_identity(edit.atom_1), atom_identity(edit.atom_2))
        if primary_identity not in identities:
            continue
        partner_identity = (
            identities[1]
            if identities[0] == primary_identity
            else identities[0]
        )
        partner = transition_by_identity.get(partner_identity)
        if (
            partner is None
            or partner.before_state is None
            or partner.after_state is None
            or partner.before_state.component_index == primary_component
        ):
            continue
        partner_labels.add(partner.before_state.concise_label)

    reactant_terms = [before_label]
    if partner_labels:
        reactant_terms.extend(sorted(partner_labels))
    return f"{' + '.join(reactant_terms)} → {after_label}"


def _active_neighbor_display(
    molecule: Any,
    center_index: int,
    neighbor_index: int,
) -> str:
    neighbor = molecule.GetAtomWithIdx(int(neighbor_index))
    symbol = neighbor.GetSymbol()
    if neighbor.GetIsAromatic():
        return "ArC" if symbol == "C" else symbol
    external = [
        atom
        for atom in neighbor.GetNeighbors()
        if atom.GetAtomicNum() > 1 and int(atom.GetIdx()) != int(center_index)
    ]
    if symbol == "C":
        carbonyl_atoms = tuple(
            atom
            for atom in external
            if str(
                molecule.GetBondBetweenAtoms(
                    int(neighbor_index),
                    int(atom.GetIdx()),
                ).GetBondType()
            ).upper()
            == "DOUBLE"
            and atom.GetSymbol() in {"N", "O", "S"}
        )
        if carbonyl_atoms:
            substituents = []
            for atom in external:
                bond = molecule.GetBondBetweenAtoms(
                    int(neighbor_index),
                    int(atom.GetIdx()),
                )
                prefix = _bond_prefix(str(bond.GetBondType()).upper())
                value = atom.GetSymbol()
                onward = tuple(
                    other
                    for other in atom.GetNeighbors()
                    if (
                        other.GetAtomicNum() > 1
                        and int(other.GetIdx()) != int(neighbor_index)
                    )
                )
                if onward and all(other.GetAtomicNum() == 6 for other in onward):
                    value = f"{value}-R"
                elif int(atom.GetTotalNumHs(includeNeighbors=True)) > 0:
                    value = f"{value}-H"
                substituents.append(f"({prefix}{value})")
            return f"C{''.join(sorted(substituents))}"
        return "R"
    if symbol not in {"N", "O", "P", "S"}:
        return symbol
    if external:
        if all(atom.GetAtomicNum() == 6 for atom in external):
            return f"{symbol}-R"
        suffix = ",".join(
            sorted(
                "R" if atom.GetAtomicNum() == 6 else atom.GetSymbol()
                for atom in external
            )
        )
        return f"{symbol}-({suffix})"
    if int(neighbor.GetTotalNumHs(includeNeighbors=True)) > 0:
        return f"{symbol}-H"
    return symbol


def state_label(
    *,
    molecule: Any,
    component_index: int,
    atom_index: int,
    active_coordinates: set[Coordinate],
    remote_classes: Mapping[
        tuple[str, int, int, int],
        ReactionCoreRemoteClass,
    ],
    side: str,
) -> str:
    """Render one active atom and its local/remote substituent context."""
    atom = molecule.GetAtomWithIdx(int(atom_index))
    tokens = ["H"] * int(atom.GetTotalNumHs(includeNeighbors=True))
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() <= 1:
            continue
        neighbor_index = int(neighbor.GetIdx())
        order = str(
            molecule.GetBondBetweenAtoms(atom_index, neighbor_index).GetBondType()
        ).upper()
        if (component_index, neighbor_index) in active_coordinates:
            value = _active_neighbor_display(molecule, atom_index, neighbor_index)
        else:
            remote_class = remote_classes.get(
                (side, component_index, atom_index, neighbor_index),
                "generic_R",
            )
            value = (
                neighbor.GetSymbol()
                if remote_class == "heteroatom"
                else _remote_display(remote_class)
            )
        tokens.append(f"{_bond_prefix(order)}{value}")

    def token_order(token: str) -> tuple[int, str]:
        plain = token.lstrip("=#:~")
        if plain == "H":
            return 0, token
        if plain in {
            "Ar",
            "HetAr",
            "R",
            "Alkenyl",
            "Alkynyl",
            "Acyl",
        }:
            return 1, token
        if token.startswith(("=", "#", ":", "~")):
            return 3, token
        return 2, token

    counts = Counter(tokens)
    rendered = []
    for token in sorted(counts, key=token_order):
        count = counts[token]
        rendered.append(f"({token}){count if count > 1 else ''}")
    return f"{atom.GetSymbol()}{''.join(rendered)}"


def _bond_change_label(edit: ReactionEdit) -> str | None:
    if edit.atom_2 is None:
        return None
    before = _edit_bond_display(edit.atom_1, edit.atom_2, edit.old_order)
    after = _edit_bond_display(edit.atom_1, edit.atom_2, edit.new_order)
    if edit.edit_type == "formed":
        return f"formed: {after}"
    if edit.edit_type == "broken":
        return f"broken: {before}"
    if edit.edit_type == "order_changed":
        return f"order changed: {before} → {after}"
    return None


def _state_change_label(change: ReactionCoreStateChange) -> str:
    center = "–".join(change.elements) or "atom"
    names = {
        "hydrogen": "hydrogen count",
        "formal_charge": "formal charge",
        "radical": "radical electrons",
        "isotope": "isotope",
        "aromaticity": "aromaticity",
        "hybridization": "hybridization",
        "atom_stereochemistry": "stereochemistry",
        "bond_stereochemistry": "bond stereochemistry",
    }
    return (
        f"{center} {names[change.change_type]}: "
        f"{change.before_value} → {change.after_value}"
    )


def _remote_context_label(subgraph: ReactionCoreRemoteSubgraph) -> str:
    label = _remote_display(subgraph.remote_class)
    details = []
    if subgraph.fragment_smiles:
        details.append(subgraph.fragment_smiles)
    if subgraph.functional_group_ids:
        details.append("groups=" + ",".join(subgraph.functional_group_ids))
    return f"{label} ({'; '.join(details)})" if details else label


def build_core_presentation(
    *,
    equation: str,
    edits: Sequence[ReactionEdit],
    state_changes: Sequence[ReactionCoreStateChange],
    remote_subgraphs: Sequence[ReactionCoreRemoteSubgraph],
    evidence_status: str,
    quality: ReactionCoreQuality,
) -> ReactionCorePresentation:
    """Render a concise audit view without affecting chemical identity."""
    bond_changes = tuple(
        label
        for edit in edits
        if (label := _bond_change_label(edit)) is not None
    )
    contexts = {
        continuity: tuple(
            sorted(
                {
                    _remote_context_label(subgraph)
                    for subgraph in remote_subgraphs
                    if subgraph.continuity == continuity
                }
            )
        )
        for continuity in ("retained", "departing", "appearing")
    }
    evidence_label = {
        "verified": "verified structural evidence",
        "inferred": "inferred structural correspondence",
        "external": "external atom mapping",
        "hypothesis": "review-only structural hypothesis",
    }.get(evidence_status, evidence_status)
    quality_label = (
        "checks passed"
        if quality.status == "pass"
        else "review: " + ", ".join(quality.review_reasons)
        if quality.status == "review"
        else "blocked: " + ", ".join(quality.blocking_reasons)
    )
    return ReactionCorePresentation(
        equation=equation,
        bond_changes=bond_changes,
        atom_state_changes=tuple(
            _state_change_label(change) for change in state_changes
        ),
        retained_context=contexts["retained"],
        departing_context=contexts["departing"],
        appearing_context=contexts["appearing"],
        evidence_label=evidence_label,
        quality_label=quality_label,
    )


__all__ = [
    "build_core_presentation",
    "multi_center_edit_label",
    "single_center_transition_label",
    "state_label",
]
