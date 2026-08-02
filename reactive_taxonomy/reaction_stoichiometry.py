"""Conservative normalization of omitted reactant stoichiometric copies."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from itertools import combinations_with_replacement
from typing import Iterable, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import ReactionComponent


@dataclass(frozen=True)
class ReactantMultiplicityNormalization:
    """Reactant inventory plus any uniquely inferred whole-component copies."""

    reactants: Tuple[ReactionComponent, ...]
    inferred_copy_count: int = 0
    warnings: Tuple[str, ...] = ()


def _element_counts(
    components: Iterable[ReactionComponent],
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        counts.update(
            str(atom.GetSymbol())
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() > 1
        )
    return counts


def _positive_difference(
    left: Counter[str], right: Counter[str]
) -> Counter[str]:
    return Counter(
        {
            element: count - right.get(element, 0)
            for element, count in left.items()
            if count > right.get(element, 0)
        }
    )


def _composition(component: ReactionComponent) -> Counter[str]:
    return _element_counts((component,))


def infer_reactant_multiplicity(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_inferred_copies: int = 3,
) -> ReactantMultiplicityNormalization:
    """Infer uniquely required copies of an entirely supplied component.

    Reaction SMILES often omit stoichiometric repetition.  A copy is added only
    when one bounded multiset of complete, unmapped reactant components exactly
    accounts for the product-side heavy-element deficit.  This does not infer a
    reaction family or atom correspondence; the normal correspondence stage
    must still independently validate the expanded inventory.
    """
    if (
        not reactants
        or not products
        or max_inferred_copies < 1
        or any(component.atom_mapped for component in reactants + products)
    ):
        return ReactantMultiplicityNormalization(reactants)

    deficit = _positive_difference(
        _element_counts(products),
        _element_counts(reactants),
    )
    if not deficit:
        return ReactantMultiplicityNormalization(reactants)

    candidates_by_structure: dict[str, ReactionComponent] = {}
    for component in reactants:
        composition = _composition(component)
        if not composition or not all(
            element in deficit and count <= deficit[element]
            for element, count in composition.items()
        ):
            continue
        candidates_by_structure.setdefault(component.canonical_smiles, component)
    candidates = tuple(
        candidates_by_structure[key] for key in sorted(candidates_by_structure)
    )
    solutions: set[tuple[str, ...]] = set()
    for copy_count in range(1, max_inferred_copies + 1):
        for selected in combinations_with_replacement(candidates, copy_count):
            total: Counter[str] = Counter()
            for component in selected:
                total.update(_composition(component))
            if total == deficit:
                solutions.add(
                    tuple(sorted(component.canonical_smiles for component in selected))
                )
        if solutions:
            break
    if len(solutions) != 1:
        return ReactantMultiplicityNormalization(reactants)

    source_structures = next(iter(solutions))
    by_index = {component.component_index: component for component in reactants}
    next_index = max(by_index, default=-1) + 1
    copies = []
    warnings = {"INFERRED_REACTANT_MULTIPLICITY"}
    for source_structure in source_structures:
        source = candidates_by_structure[source_structure]
        source_index = source.component_index
        copy_index = next_index + len(copies)
        copies.append(
            replace(
                source,
                component_index=copy_index,
                inferred_copy_of_component_index=source_index,
            )
        )
        warnings.add(f"INFERRED_REACTANT_COPY:c{source_index}->c{copy_index}")
    return ReactantMultiplicityNormalization(
        reactants=reactants + tuple(copies),
        inferred_copy_count=len(copies),
        warnings=tuple(sorted(warnings)),
    )


__all__ = [
    "ReactantMultiplicityNormalization",
    "infer_reactant_multiplicity",
]
