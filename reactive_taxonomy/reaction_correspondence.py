"""Conservative graph correspondence for unmapped reaction transformations."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .chemistry.smarts_cache import compile_smarts
from .reaction_models import ReactionComponent


AtomPair = Tuple[int, int, int, int]
REACTION_CORRESPONDENCE_VERSION = "2.2"


@dataclass(frozen=True)
class ScaffoldCorrespondenceCandidates:
    """Candidate heavy-atom correspondences for one conserved scaffold."""

    candidates: Tuple[Tuple[AtomPair, ...], ...]
    warnings: Tuple[str, ...] = ()
    valid: bool = False


def _heavy_atom_count(component: ReactionComponent) -> int:
    mol = parse_smiles(component.input_smiles)
    return int(mol.GetNumHeavyAtoms()) if mol is not None else 0


def _element_counts(
    components: Tuple[ReactionComponent, ...],
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        counts.update(
            atom.GetSymbol()
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() > 1
        )
    return counts


def _component_correspondence_candidates(
    reactant: ReactionComponent,
    product: ReactionComponent,
    *,
    max_matches: int,
) -> tuple[Tuple[Tuple[AtomPair, ...], ...], Tuple[str, ...]]:
    """Return best local mappings from one reactant component into a product."""
    reactant_mol = parse_smiles(reactant.input_smiles)
    product_mol = parse_smiles(product.input_smiles)
    if reactant_mol is None or product_mol is None:
        return (), ("GLOBAL_CORRESPONDENCE_PARSE_FAILED",)
    reactant_heavy_count = int(reactant_mol.GetNumHeavyAtoms())
    product_heavy_count = int(product_mol.GetNumHeavyAtoms())
    if not reactant_heavy_count or not product_heavy_count:
        return (), ()

    if reactant_heavy_count == 1:
        reactant_atom = next(
            atom for atom in reactant_mol.GetAtoms() if atom.GetAtomicNum() > 1
        )
        candidates = tuple(
            (
                (
                    reactant.component_index,
                    int(reactant_atom.GetIdx()),
                    product.component_index,
                    int(product_atom.GetIdx()),
                ),
            )
            for product_atom in product_mol.GetAtoms()
            if product_atom.GetAtomicNum() == reactant_atom.GetAtomicNum()
        )
        if len(candidates) > max_matches:
            return (), ("GLOBAL_CORRESPONDENCE_COMPONENT_MATCH_LIMIT",)
        return candidates, ()

    from rdkit.Chem import rdFMCS

    result = rdFMCS.FindMCS(
        [reactant_mol, product_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        matchValences=False,
        timeout=5,
    )
    if result.canceled:
        return (), ("GLOBAL_CORRESPONDENCE_TIMEOUT",)
    minimum_conservation = min(reactant_heavy_count, product_heavy_count)
    if result.numAtoms < 2 or (
        result.numAtoms / max(minimum_conservation, 1) < 0.4
        and result.numAtoms < 4
    ):
        return (), ()
    query = compile_smarts(result.smartsString, validate=False)
    if query is None:
        return (), ("GLOBAL_CORRESPONDENCE_QUERY_FAILED",)
    reactant_matches = reactant_mol.GetSubstructMatches(
        query, uniquify=False, maxMatches=max_matches + 1
    )
    product_matches = product_mol.GetSubstructMatches(
        query, uniquify=False, maxMatches=max_matches + 1
    )
    if len(reactant_matches) > max_matches or len(product_matches) > max_matches:
        return (), ("GLOBAL_CORRESPONDENCE_COMPONENT_MATCH_LIMIT",)
    candidates = set()
    for reactant_match in reactant_matches:
        for product_match in product_matches:
            mapping = tuple(
                sorted(
                    (
                        reactant.component_index,
                        int(reactant_atom),
                        product.component_index,
                        int(product_atom),
                    )
                    for reactant_atom, product_atom in zip(
                        reactant_match, product_match
                    )
                    if reactant_mol.GetAtomWithIdx(
                        int(reactant_atom)
                    ).GetAtomicNum()
                    > 1
                    and product_mol.GetAtomWithIdx(
                        int(product_atom)
                    ).GetAtomicNum()
                    > 1
                )
            )
            if mapping:
                candidates.add(mapping)
    if not candidates:
        return (), ()
    return (
        tuple(sorted(candidates)),
        (),
    )


def infer_global_correspondence_candidates(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_component_matches: int = 64,
    max_combinations: int = 20_000,
    max_candidates: int = 512,
    max_product_heavy_atoms: int = 80,
    max_reactant_components: int = 6,
) -> ScaffoldCorrespondenceCandidates:
    """Map multiple reactant components into one product without family routing.

    Every product heavy atom must be accounted for, at least two reactant
    components must contribute, and overlapping product assignments are
    forbidden. Ambiguous chemistry is resolved later by edit-set consensus.
    """
    if any(component.atom_mapped for component in reactants + products):
        return ScaffoldCorrespondenceCandidates(
            (), ("GLOBAL_CORRESPONDENCE_SKIPPED_MAPPED_INPUT",), False
        )
    substantial_products = tuple(
        component for component in products if _heavy_atom_count(component) >= 2
    )
    if len(substantial_products) != 1 or any(
        _heavy_atom_count(component) > 0
        for component in products
        if component is not substantial_products[0]
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("GLOBAL_CORRESPONDENCE_REQUIRES_ONE_PRODUCT",), False
        )
    product = substantial_products[0]
    product_mol = parse_smiles(product.input_smiles)
    if product_mol is None:
        return ScaffoldCorrespondenceCandidates(
            (), ("GLOBAL_CORRESPONDENCE_PARSE_FAILED",), False
        )
    if (
        int(product_mol.GetNumHeavyAtoms()) > max_product_heavy_atoms
        or len(reactants) > max_reactant_components
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("GLOBAL_CORRESPONDENCE_COMPLEXITY_LIMIT",), False
        )
    reactant_elements = _element_counts(reactants)
    product_elements = _element_counts((product,))
    if any(
        count > reactant_elements.get(element, 0)
        for element, count in product_elements.items()
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("GLOBAL_CORRESPONDENCE_UNACCOUNTED_PRODUCT_ATOMS",), False
        )
    product_heavy_atoms = frozenset(
        int(atom.GetIdx())
        for atom in product_mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    )
    options: Dict[
        int, Tuple[Tuple[AtomPair, ...], ...]
    ] = {}
    warnings = []
    for reactant in reactants:
        candidates, candidate_warnings = _component_correspondence_candidates(
            reactant,
            product,
            max_matches=max_component_matches,
        )
        warnings.extend(candidate_warnings)
        if candidates:
            options[reactant.component_index] = candidates
    if any(
        warning in {
            "GLOBAL_CORRESPONDENCE_COMPONENT_MATCH_LIMIT",
            "GLOBAL_CORRESPONDENCE_TIMEOUT",
        }
        for warning in warnings
    ):
        return ScaffoldCorrespondenceCandidates(
            (), tuple(sorted(set(warnings))), False
        )
    if len(options) < 2:
        return ScaffoldCorrespondenceCandidates(
            (),
            tuple(
                sorted(
                    set(warnings).union(
                        {"GLOBAL_CORRESPONDENCE_REQUIRES_MULTIPLE_PARTNERS"}
                    )
                )
            ),
            False,
        )
    ordered_components = tuple(
        sorted(
            options,
            key=lambda component_index: (
                -max(len(mapping) for mapping in options[component_index]),
                component_index,
            ),
        )
    )
    remaining_coverage = []
    running: set[int] = set()
    for component_index in reversed(ordered_components):
        running = running.union(
            product_atom
            for mapping in options[component_index]
            for _, _, _, product_atom in mapping
        )
        remaining_coverage.append(frozenset(running))
    remaining_coverage.reverse()

    results = set()
    attempted = 0
    limit_reached = False

    def search(
        position: int,
        covered: frozenset[int],
        selected: Tuple[Tuple[AtomPair, ...], ...],
    ) -> None:
        nonlocal attempted, limit_reached
        if limit_reached:
            return
        attempted += 1
        if attempted > max_combinations:
            limit_reached = True
            return
        if covered == product_heavy_atoms:
            if len(selected) >= 2:
                results.add(tuple(sorted(pair for mapping in selected for pair in mapping)))
                if len(results) > max_candidates:
                    limit_reached = True
            return
        if position >= len(ordered_components):
            return
        if not product_heavy_atoms <= covered.union(
            remaining_coverage[position]
        ):
            return
        component_index = ordered_components[position]
        search(position + 1, covered, selected)
        for mapping in options[component_index]:
            mapped_product_atoms = frozenset(
                product_atom for _, _, _, product_atom in mapping
            )
            if covered.intersection(mapped_product_atoms):
                continue
            search(
                position + 1,
                covered.union(mapped_product_atoms),
                selected + (mapping,),
            )

    search(0, frozenset(), ())
    if limit_reached:
        return ScaffoldCorrespondenceCandidates(
            (),
            tuple(
                sorted(
                    set(warnings).union(
                        {"GLOBAL_CORRESPONDENCE_SEARCH_LIMIT"}
                    )
                )
            ),
            False,
        )
    if not results:
        return ScaffoldCorrespondenceCandidates(
            (),
            tuple(
                sorted(
                    set(warnings).union({"GLOBAL_CORRESPONDENCE_NOT_FOUND"})
                )
            ),
            False,
        )
    return ScaffoldCorrespondenceCandidates(
        tuple(sorted(results)),
        tuple(sorted(set(warnings))),
        True,
    )


def infer_scaffold_correspondence_candidates(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_candidates: int = 512,
) -> ScaffoldCorrespondenceCandidates:
    """Enumerate mappings only for a uniquely conserved, single scaffold.

    The fallback intentionally excludes multi-substrate assembly. Every product
    heavy atom must belong to the conserved scaffold, while other reactant
    components may contain at most one heavy atom (for example H2, Na, or O).
    """
    if any(component.atom_mapped for component in reactants + products):
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_SKIPPED_MAPPED_INPUT",), False
        )
    substantial_reactants = tuple(
        component for component in reactants if _heavy_atom_count(component) >= 2
    )
    substantial_products = tuple(
        component for component in products if _heavy_atom_count(component) >= 2
    )
    if len(substantial_reactants) != 1:
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_REQUIRES_ONE_SUBSTRATE",), False
        )
    if len(substantial_products) != 1 or any(
        _heavy_atom_count(component) > 0
        for component in products
        if component is not substantial_products[0]
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_REQUIRES_ONE_PRODUCT",), False
        )
    reactant = substantial_reactants[0]
    product = substantial_products[0]
    reactant_mol = parse_smiles(reactant.input_smiles)
    product_mol = parse_smiles(product.input_smiles)
    if reactant_mol is None or product_mol is None:
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_PARSE_FAILED",), False
        )

    from rdkit.Chem import rdFMCS

    result = rdFMCS.FindMCS(
        [reactant_mol, product_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        matchValences=False,
        timeout=5,
    )
    if result.canceled:
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_TIMEOUT",), False
        )
    query = compile_smarts(result.smartsString, validate=False)
    if query is None:
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_QUERY_FAILED",), False
        )
    product_heavy_atoms = {
        atom.GetIdx() for atom in product_mol.GetAtoms() if atom.GetAtomicNum() > 1
    }
    reactant_heavy_count = int(reactant_mol.GetNumHeavyAtoms())
    if result.numAtoms < len(product_heavy_atoms) or (
        result.numAtoms / max(reactant_heavy_count, 1) < 0.6
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_INSUFFICIENT_CONSERVATION",), False
        )
    reactant_matches = reactant_mol.GetSubstructMatches(
        query, uniquify=False, maxMatches=max_candidates + 1
    )
    product_matches = product_mol.GetSubstructMatches(
        query, uniquify=False, maxMatches=max_candidates + 1
    )
    candidates = set()
    for reactant_match in reactant_matches:
        for product_match in product_matches:
            if not product_heavy_atoms <= set(product_match):
                continue
            mapping = tuple(
                sorted(
                    (
                        reactant.component_index,
                        int(reactant_atom),
                        product.component_index,
                        int(product_atom),
                    )
                    for reactant_atom, product_atom in zip(
                        reactant_match, product_match
                    )
                    if reactant_mol.GetAtomWithIdx(reactant_atom).GetAtomicNum() > 1
                    and product_mol.GetAtomWithIdx(product_atom).GetAtomicNum() > 1
                )
            )
            candidates.add(mapping)
            if len(candidates) > max_candidates:
                return ScaffoldCorrespondenceCandidates(
                    (), ("SCAFFOLD_CORRESPONDENCE_CANDIDATE_LIMIT",), False
                )
    if not candidates:
        return ScaffoldCorrespondenceCandidates(
            (), ("SCAFFOLD_CORRESPONDENCE_NOT_FOUND",), False
        )
    return ScaffoldCorrespondenceCandidates(
        tuple(sorted(candidates)), (), True
    )


def infer_partial_scaffold_correspondence_candidates(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_candidates: int = 512,
) -> ScaffoldCorrespondenceCandidates:
    """Map a conserved scaffold while allowing unmatched terminal attachments.

    Unlike verified scaffold correspondence, this observation-only path does
    not require every product atom to have a reported reactant source. It is
    intentionally limited to one substantial substrate and one product; the
    caller must validate any inferred local transformation across all mappings.
    """
    substantial_reactants = tuple(
        component for component in reactants if _heavy_atom_count(component) >= 2
    )
    substantial_products = tuple(
        component for component in products if _heavy_atom_count(component) >= 2
    )
    if len(substantial_reactants) != 1:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_REQUIRES_ONE_SUBSTRATE",), False
        )
    if len(substantial_products) != 1 or any(
        _heavy_atom_count(component) > 0
        for component in products
        if component is not substantial_products[0]
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_REQUIRES_ONE_PRODUCT",), False
        )
    reactant = substantial_reactants[0]
    product = substantial_products[0]
    reactant_mol = parse_smiles(reactant.input_smiles)
    product_mol = parse_smiles(product.input_smiles)
    if reactant_mol is None or product_mol is None:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_PARSE_FAILED",), False
        )
    reactant_elements = _element_counts((reactant,))
    product_elements = _element_counts((product,))
    removed_elements = tuple(
        element
        for element, count in sorted(reactant_elements.items())
        for _ in range(max(count - product_elements.get(element, 0), 0))
    )
    added_elements = tuple(
        element
        for element, count in sorted(product_elements.items())
        for _ in range(max(count - reactant_elements.get(element, 0), 0))
    )
    if len(removed_elements) != 1 or len(added_elements) != 1:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_REQUIRES_ONE_ATOM_EXCHANGE",), False
        )

    from rdkit import Chem

    def terminal_atoms(molecule: object, element: str) -> Tuple[int, ...]:
        return tuple(
            int(atom.GetIdx())
            for atom in molecule.GetAtoms()
            if atom.GetSymbol() == element
            and atom.GetAtomicNum() > 1
            and sum(
                neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors()
            )
            == 1
        )

    def delete_atom(
        molecule: object, atom_index: int
    ) -> tuple[object, Tuple[int, ...], str] | None:
        original_indices = tuple(
            index
            for index in range(int(molecule.GetNumAtoms()))
            if index != atom_index
        )
        editable = Chem.RWMol(molecule)
        editable.RemoveAtom(int(atom_index))
        reduced = editable.GetMol()
        try:
            Chem.SanitizeMol(reduced)
            for atom in reduced.GetAtoms():
                atom.SetAtomMapNum(0)
            Chem.AssignStereochemistry(reduced, cleanIt=True, force=True)
            canonical = Chem.MolToSmiles(
                reduced, canonical=True, isomericSmiles=True
            )
        except Exception:
            return None
        return reduced, original_indices, canonical

    candidates = set()
    removed_candidates = terminal_atoms(reactant_mol, removed_elements[0])
    added_candidates = terminal_atoms(product_mol, added_elements[0])
    for removed_index in removed_candidates:
        reactant_reduced = delete_atom(reactant_mol, removed_index)
        if reactant_reduced is None:
            continue
        reactant_graph, reactant_original, reactant_canonical = reactant_reduced
        for added_index in added_candidates:
            product_reduced = delete_atom(product_mol, added_index)
            if product_reduced is None:
                continue
            product_graph, product_original, product_canonical = product_reduced
            if reactant_canonical != product_canonical:
                continue
            matches = product_graph.GetSubstructMatches(
                reactant_graph,
                uniquify=False,
                useChirality=True,
                maxMatches=max_candidates + 1,
            )
            if len(matches) > max_candidates:
                return ScaffoldCorrespondenceCandidates(
                    (), ("PARTIAL_CORRESPONDENCE_CANDIDATE_LIMIT",), False
                )
            for match in matches:
                mapping = tuple(
                    sorted(
                        (
                            reactant.component_index,
                            int(reactant_original[reactant_atom]),
                            product.component_index,
                            int(product_original[product_atom]),
                        )
                        for reactant_atom, product_atom in enumerate(match)
                        if reactant_graph.GetAtomWithIdx(
                            reactant_atom
                        ).GetAtomicNum()
                        > 1
                    )
                )
                if len(mapping) >= 3:
                    candidates.add(mapping)
                    if len(candidates) > max_candidates:
                        return ScaffoldCorrespondenceCandidates(
                            (),
                            ("PARTIAL_CORRESPONDENCE_CANDIDATE_LIMIT",),
                            False,
                        )
    if not candidates:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_NOT_FOUND",), False
        )
    return ScaffoldCorrespondenceCandidates(
        tuple(sorted(candidates)), (), True
    )


__all__ = [
    "AtomPair",
    "REACTION_CORRESPONDENCE_VERSION",
    "ScaffoldCorrespondenceCandidates",
    "infer_global_correspondence_candidates",
    "infer_partial_scaffold_correspondence_candidates",
    "infer_scaffold_correspondence_candidates",
]
