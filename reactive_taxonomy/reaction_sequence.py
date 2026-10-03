"""Exact structural continuity for supplied, forward reaction sequences.

This validates connectivity only, not reaction feasibility or atom mapping.
No salts, protonation states, tautomers or stereoisomers are conflated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from rdkit import Chem

from .chemistry.rdkit_utils import parse_smiles


@dataclass(frozen=True)
class ConnectedReactionStep:
    """Canonical components and selected occurrences in a linear route."""

    reactants: tuple[str, ...]
    agents: tuple[str, ...]
    products: tuple[str, ...]
    carried_reactant_index: int | None
    product_index: int


def _components(section: str) -> tuple[str, ...]:
    if not section:
        return ()
    result = []
    for text in section.split("."):
        mol = parse_smiles(text)
        if mol is None or not mol.GetNumAtoms() or any(c.isspace() for c in text):
            raise ValueError(f"Invalid SMILES component: {text!r}")
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        result.append(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    return tuple(result)


def _index(index: int | None, count: int, label: str) -> None:
    if index is not None and (type(index) is not int or not 0 <= index < count):
        raise ValueError(f"Invalid {label}: {index!r}")


def connect_reaction_sequence(
    reactions: Sequence[str],
    *,
    main_reactant_index: int | None = None,
    product_indices: Sequence[int | None] | None = None,
) -> tuple[ConnectedReactionStep, ...]:
    """Connect exact product/reactant identities, ignoring map numbers.

    The first starting-material block includes all reactants unless an explicit
    zero-based main reactant index is supplied. Multiproduct steps require a
    unique next-step match or an explicit product index. Other products remain
    in the result. Duplicate matching reactant occurrences are ambiguous.
    """
    if isinstance(reactions, str) or not reactions:
        raise ValueError("Supply a nonempty ordered sequence of reaction SMILES")
    selections = (
        tuple(product_indices)
        if product_indices is not None
        else (None,) * len(reactions)
    )
    if len(selections) != len(reactions):
        raise ValueError("product_indices must have one entry per reaction")
    parsed = []
    for number, reaction in enumerate(reactions, 1):
        if not isinstance(reaction, str):
            raise ValueError(f"Step {number}: reaction SMILES must be text")
        sections = reaction.strip().split(">")
        if len(sections) != 3 or not sections[0] or not sections[2]:
            raise ValueError(f"Step {number}: expected reactants>agents>products")
        try:
            parsed.append(tuple(_components(section) for section in sections))
        except ValueError as exc:
            raise ValueError(f"Step {number}: {exc}") from exc
    _index(main_reactant_index, len(parsed[0][0]), "main_reactant_index")
    output = []
    previous = None
    for i, (reactants, agents, products) in enumerate(parsed):
        carried = main_reactant_index if i == 0 else None
        if previous is not None:
            matches = [j for j, smiles in enumerate(reactants) if smiles == previous]
            if len(matches) != 1:
                reason = (
                    "disconnected" if not matches else "ambiguous duplicate reactant"
                )
                raise ValueError(
                    f"Step {i + 1}: {reason}; previous product must match exactly one reactant"
                )
            carried = matches[0]
        selected = selections[i]
        _index(selected, len(products), f"product index at step {i + 1}")
        if selected is None:
            candidates = list(range(len(products)))
            if len(products) > 1 and i + 1 < len(parsed):
                candidates = [
                    j for j, smiles in enumerate(products) if smiles in parsed[i + 1][0]
                ]
            if len(candidates) != 1:
                raise ValueError(
                    f"Step {i + 1}: ambiguous product selection; supply product_indices"
                )
            selected = candidates[0]
        output.append(
            ConnectedReactionStep(reactants, agents, products, carried, selected)
        )
        previous = products[selected]
    return tuple(output)
