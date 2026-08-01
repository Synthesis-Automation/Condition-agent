"""Reusable graph predicates for optional reaction-pattern interpretation.

The helpers in this module consume only an existing ``ReactionObservation``.
They cannot propose atom correspondence, edits, products, or reaction cores.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import ReactionAtomReference, ReactionEdit, ReactionObservation


_ORDER_RANK = {
    "NONE": 0.0,
    "SINGLE": 1.0,
    "AROMATIC": 1.5,
    "DOUBLE": 2.0,
    "TRIPLE": 3.0,
}


AtomKey = tuple[int, int]
EditPredicate = Callable[[ReactionEdit], bool]


@dataclass
class ReactionPatternContext:
    """Cached structural queries over one normalized reaction observation."""

    observation: ReactionObservation
    _molecules: dict[tuple[str, int], Any] = field(default_factory=dict)

    @property
    def edits(self) -> tuple[ReactionEdit, ...]:
        """Return the observation's normalized edits."""
        return self.observation.edits

    def indices(self, predicate: EditPredicate) -> tuple[int, ...]:
        """Return sorted edit indices satisfying ``predicate``."""
        return tuple(
            index for index, edit in enumerate(self.edits) if predicate(edit)
        )

    def molecule(self, side: str, component_index: int) -> Any | None:
        """Return a lazily parsed component molecule."""
        key = (side, component_index)
        if key not in self._molecules:
            components = {
                "reactant": self.observation.reactants,
                "agent": self.observation.agents,
                "product": self.observation.products,
            }.get(side, ())
            component = next(
                (
                    item
                    for item in components
                    if item.component_index == component_index
                ),
                None,
            )
            self._molecules[key] = (
                parse_smiles(component.input_smiles) if component is not None else None
            )
        return self._molecules[key]

    def atom(self, reference: ReactionAtomReference | None) -> Any | None:
        """Resolve an atom reference into its source molecule."""
        if reference is None:
            return None
        molecule = self.molecule(reference.side, reference.component_index)
        if molecule is None or not 0 <= reference.atom_index < molecule.GetNumAtoms():
            return None
        return molecule.GetAtomWithIdx(reference.atom_index)

    def atom_key(
        self,
        edit: ReactionEdit,
        endpoint: int,
        *,
        side: str = "reactant",
    ) -> AtomKey | None:
        """Return a stable component/atom key for one edit endpoint."""
        atom = edit.atom_1 if endpoint == 1 else edit.atom_2
        if atom is None or atom.side != side:
            return None
        return atom.component_index, atom.atom_index

    @staticmethod
    def element_pair(edit: ReactionEdit) -> frozenset[str]:
        """Return the unordered element set for one edit."""
        return frozenset(
            atom.element for atom in (edit.atom_1, edit.atom_2) if atom is not None
        )

    @staticmethod
    def order_direction(edit: ReactionEdit) -> int:
        """Return -1, 0, or +1 for a normalized bond-order edit."""
        old = _ORDER_RANK.get(str(edit.old_order or "NONE").upper())
        new = _ORDER_RANK.get(str(edit.new_order or "NONE").upper())
        if old is None or new is None:
            return 0
        return (new > old) - (new < old)

    def endpoint_for_element(
        self,
        edit: ReactionEdit,
        element: str,
    ) -> tuple[int, ReactionAtomReference, AtomKey] | None:
        """Resolve the unique reactant endpoint having ``element``."""
        matches = []
        for endpoint, atom in ((1, edit.atom_1), (2, edit.atom_2)):
            key = self.atom_key(edit, endpoint)
            if atom is not None and atom.element == element and key is not None:
                matches.append((endpoint, atom, key))
        return matches[0] if len(matches) == 1 else None

    def has_neighbor(
        self,
        key: AtomKey,
        *,
        element: str,
        order: str | None = None,
        exclude_atom_indices: Sequence[int] = (),
    ) -> bool:
        """Whether a reactant atom has a matching bonded neighbor."""
        molecule = self.molecule("reactant", key[0])
        if molecule is None or not 0 <= key[1] < molecule.GetNumAtoms():
            return False
        atom = molecule.GetAtomWithIdx(key[1])
        excluded = set(exclude_atom_indices)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in excluded or neighbor.GetSymbol() != element:
                continue
            bond = molecule.GetBondBetweenAtoms(key[1], neighbor.GetIdx())
            if order is None or str(bond.GetBondType()).upper() == order.upper():
                return True
        return False

    def neighbor_elements(
        self,
        key: AtomKey,
        *,
        order: str | None = None,
    ) -> tuple[str, ...]:
        """Return sorted reactant-neighbor elements, optionally by bond order."""
        molecule = self.molecule("reactant", key[0])
        if molecule is None or not 0 <= key[1] < molecule.GetNumAtoms():
            return ()
        elements = []
        for neighbor in molecule.GetAtomWithIdx(key[1]).GetNeighbors():
            bond = molecule.GetBondBetweenAtoms(key[1], neighbor.GetIdx())
            if order is None or str(bond.GetBondType()).upper() == order.upper():
                elements.append(str(neighbor.GetSymbol()))
        return tuple(sorted(elements))

    def is_carbonyl_carbon(self, key: AtomKey) -> bool:
        """Whether ``key`` is a reactant carbon double-bonded to oxygen."""
        molecule = self.molecule("reactant", key[0])
        if molecule is None or not 0 <= key[1] < molecule.GetNumAtoms():
            return False
        atom = molecule.GetAtomWithIdx(key[1])
        return atom.GetSymbol() == "C" and self.has_neighbor(
            key, element="O", order="DOUBLE"
        )

    def is_sulfonyl_sulfur(self, key: AtomKey) -> bool:
        """Whether ``key`` is sulfur with at least two S=O bonds."""
        molecule = self.molecule("reactant", key[0])
        if molecule is None or not 0 <= key[1] < molecule.GetNumAtoms():
            return False
        atom = molecule.GetAtomWithIdx(key[1])
        if atom.GetSymbol() != "S":
            return False
        double_oxygen_count = 0
        for neighbor in atom.GetNeighbors():
            bond = molecule.GetBondBetweenAtoms(key[1], neighbor.GetIdx())
            if (
                neighbor.GetSymbol() == "O"
                and str(bond.GetBondType()).upper() == "DOUBLE"
            ):
                double_oxygen_count += 1
        return double_oxygen_count >= 2

    def is_sp2_carbon(self, reference: ReactionAtomReference) -> bool:
        """Whether a carbon endpoint is aromatic, vinylic, or otherwise sp2."""
        return reference.element == "C" and (
            reference.aromatic or reference.hybridization.upper() == "SP2"
        )

    def is_sp3_carbon(self, reference: ReactionAtomReference) -> bool:
        """Whether a carbon endpoint is sp3."""
        return reference.element == "C" and reference.hybridization.upper() == "SP3"

    def edits_at(self, key: AtomKey, *, edit_type: str | None = None) -> tuple[int, ...]:
        """Return normalized edits incident to one reactant atom."""
        return tuple(
            index
            for index, edit in enumerate(self.edits)
            if (edit_type is None or edit.edit_type == edit_type)
            and key
            in {
                self.atom_key(edit, 1),
                self.atom_key(edit, 2),
            }
        )


def unique_indices(*groups: Sequence[int]) -> tuple[int, ...]:
    """Return deterministic unique indices from several match clauses."""
    return tuple(sorted({index for group in groups for index in group}))


__all__ = ["AtomKey", "ReactionPatternContext", "unique_indices"]
