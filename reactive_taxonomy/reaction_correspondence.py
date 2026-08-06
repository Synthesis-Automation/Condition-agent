"""Conservative graph correspondence for unmapped reaction transformations."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .chemistry.smarts_cache import compile_smarts
from .reaction_models import ReactionComponent


AtomPair = Tuple[int, int, int, int]
REACTION_CORRESPONDENCE_VERSION = "2.6"


def _has_unsupported_mcs_bond(molecule: Any) -> bool:
    """Return whether RDKit MCS is unsafe for this coordination graph."""
    from rdkit.Chem.rdchem import BondType

    unsupported_bond_types = {
        BondType.DATIVE,
        BondType.DATIVEL,
        BondType.DATIVEONE,
        BondType.DATIVER,
    }
    return any(
        bond.GetBondType() in unsupported_bond_types
        for bond in molecule.GetBonds()
    )


def _symmetry_distinct_matches(
    molecule: Any,
    query: Any,
    *,
    max_matches: int,
) -> tuple[Tuple[Tuple[int, ...], ...], bool]:
    """Return atom-distinct matches after collapsing molecular automorphisms.

    ``uniquify=True`` is too aggressive for reaction correspondence because it
    can collapse alternative orientations that place an edit at chemically
    different atoms.  Conversely, raw RDKit matches enumerate every
    permutation of equivalent groups such as tert-butyl and trimethylsilyl
    methyls.  Canonical symmetry-class sequences retain the former while
    collapsing the latter before the chemistry search limit is applied.

    The second return value reports bounded-search overflow.  We continue to
    abstain instead of silently accepting a truncated correspondence set.
    """
    from rdkit import Chem

    raw_match_limit = max(4_096, max_matches * 128)
    raw_matches = molecule.GetSubstructMatches(
        query,
        uniquify=False,
        maxMatches=raw_match_limit + 1,
    )
    if len(raw_matches) > raw_match_limit:
        return (), True
    symmetry_classes = tuple(
        int(value)
        for value in Chem.CanonicalRankAtoms(
            molecule,
            breakTies=False,
            includeChirality=True,
        )
    )
    representatives: Dict[
        Tuple[Tuple[int, ...], Tuple[int, ...]], Tuple[int, ...]
    ] = {}
    for raw_match in raw_matches:
        match = tuple(int(atom_index) for atom_index in raw_match)
        # The atom set distinguishes separate copies/sites in the target.  The
        # symmetry-class sequence additionally retains inequivalent query
        # orientations over the same atom set (for example the two possible
        # ends of an elimination scaffold).
        signature = (
            tuple(sorted(match)),
            tuple(symmetry_classes[atom_index] for atom_index in match),
        )
        previous = representatives.get(signature)
        if previous is None or match < previous:
            representatives[signature] = match
    if len(representatives) > max_matches:
        return (), True
    return tuple(sorted(representatives.values())), False


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

    if _has_unsupported_mcs_bond(reactant_mol) or _has_unsupported_mcs_bond(
        product_mol
    ):
        return (), ("GLOBAL_CORRESPONDENCE_UNSUPPORTED_DATIVE_BOND",)

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
    reactant_matches, reactant_limit_reached = _symmetry_distinct_matches(
        reactant_mol,
        query,
        max_matches=max_matches,
    )
    product_matches, product_limit_reached = _symmetry_distinct_matches(
        product_mol,
        query,
        max_matches=max_matches,
    )
    if reactant_limit_reached or product_limit_reached:
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


def _induced_molecule(molecule: Any, atom_indices: Tuple[int, ...]) -> Any:
    """Return an index-traceable induced RDKit molecule."""
    from rdkit import Chem

    retained = set(atom_indices)
    editable = Chem.RWMol(molecule)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_correspondence_original_index", int(atom.GetIdx()))
    for atom_index in sorted(
        set(range(molecule.GetNumAtoms())) - retained,
        reverse=True,
    ):
        editable.RemoveAtom(int(atom_index))
    induced = editable.GetMol()
    Chem.GetSymmSSSR(induced)
    return induced


def _single_cut_fragment_correspondence_candidates(
    reactant: ReactionComponent,
    product: ReactionComponent,
    *,
    max_matches: int,
    max_fragments: int = 32,
) -> tuple[Tuple[Tuple[AtomPair, ...], ...], Tuple[str, ...]]:
    """Map the largest exact fragment exposed by one non-ring bond cut.

    This bounded fallback supports records where one reactant both loses a
    terminal fragment and participates in another event.  It proposes only
    atom-origin alternatives; edit minimization remains responsible for
    accepting a complete reaction correspondence.
    """
    from rdkit import Chem

    reactant_mol = parse_smiles(reactant.input_smiles)
    product_mol = parse_smiles(product.input_smiles)
    if reactant_mol is None or product_mol is None:
        return (), ("GLOBAL_CORRESPONDENCE_FRAGMENT_PARSE_FAILED",)
    candidates: set[Tuple[AtomPair, ...]] = set()
    largest_candidate_size = 0
    fragment_atom_sets: set[Tuple[int, ...]] = set()
    for bond in reactant_mol.GetBonds():
        if bond.IsInRing():
            continue
        editable = Chem.RWMol(reactant_mol)
        editable.RemoveBond(
            int(bond.GetBeginAtomIdx()),
            int(bond.GetEndAtomIdx()),
        )
        fragments = Chem.GetMolFrags(
            editable.GetMol(),
            asMols=False,
            sanitizeFrags=False,
        )
        if len(fragments) < 2:
            continue
        for fragment in fragments:
            atom_indices = tuple(
                sorted(
                    int(index)
                    for index in fragment
                    if reactant_mol.GetAtomWithIdx(int(index)).GetAtomicNum() > 1
                )
            )
            if len(atom_indices) < 2:
                continue
            fragment_atom_sets.add(atom_indices)
            if len(fragment_atom_sets) > max_fragments:
                return (), ("GLOBAL_CORRESPONDENCE_FRAGMENT_LIMIT",)

    # Only maximum-size fragment mappings are returned below.  Search larger
    # fragments first so a dominated, highly symmetric fragment cannot exhaust
    # the match limit before a more informative fragment is considered.  This
    # branch-and-bound ordering is independent of reactant bond serialization.
    ordered_fragments = tuple(
        sorted(fragment_atom_sets, key=lambda value: (-len(value), value))
    )
    for atom_indices in ordered_fragments:
        if len(atom_indices) < largest_candidate_size:
            break
        query = _induced_molecule(reactant_mol, atom_indices)
        matches = product_mol.GetSubstructMatches(
            query,
            uniquify=False,
            maxMatches=max_matches + 1,
        )
        if len(matches) > max_matches:
            return (), ("GLOBAL_CORRESPONDENCE_FRAGMENT_MATCH_LIMIT",)
        source_indices = tuple(
            int(
                query.GetAtomWithIdx(index).GetIntProp(
                    "_correspondence_original_index"
                )
            )
            for index in range(query.GetNumAtoms())
        )
        for match in matches:
            candidate = tuple(
                sorted(
                    (
                        reactant.component_index,
                        source_index,
                        product.component_index,
                        int(product_index),
                    )
                    for source_index, product_index in zip(
                        source_indices, match
                    )
                )
            )
            candidate_size = len(candidate)
            if candidate_size > largest_candidate_size:
                candidates.clear()
                largest_candidate_size = candidate_size
            candidates.add(candidate)
    if not candidates:
        return (), ()
    return tuple(sorted(candidates)), ()


def _global_correspondence_assignments(
    options: Mapping[int, Tuple[Tuple[AtomPair, ...], ...]],
    product_heavy_atoms: frozenset[int],
    *,
    max_combinations: int,
    max_candidates: int,
) -> tuple[set[Tuple[AtomPair, ...]], bool]:
    """Return non-overlapping component assignments covering the product."""
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

    results: set[Tuple[AtomPair, ...]] = set()
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
                results.add(
                    tuple(sorted(pair for mapping in selected for pair in mapping))
                )
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
    return results, limit_reached


def _fragment_correspondence_options(
    reactant_mol: Any,
    product_mol: Any,
    reactant_atoms: Tuple[int, ...],
    product_atoms: Tuple[int, ...],
    *,
    max_fragment_matches: int,
    max_singleton_product_atoms: int,
) -> tuple[str, Tuple[Tuple[Tuple[int, int], ...], ...]]:
    """Find the next conserved fragment between two residual atom sets."""
    from rdkit.Chem import rdFMCS

    product_fragment = _induced_molecule(product_mol, product_atoms)
    reactant_element_counts: Counter[str] = Counter(
        reactant_mol.GetAtomWithIdx(atom_index).GetSymbol()
        for atom_index in reactant_atoms
    )
    product_element_counts: Counter[str] = Counter(
        product_mol.GetAtomWithIdx(atom_index).GetSymbol()
        for atom_index in product_atoms
    )
    omission_candidates = tuple(
        atom_index
        for atom_index in reactant_atoms
        if reactant_element_counts[
            reactant_mol.GetAtomWithIdx(atom_index).GetSymbol()
        ]
        > product_element_counts.get(
            reactant_mol.GetAtomWithIdx(atom_index).GetSymbol(),
            0,
        )
    )
    reactant_variants = [reactant_atoms]
    if len(omission_candidates) <= 8:
        reactant_variants.extend(
            tuple(
                atom_index
                for atom_index in reactant_atoms
                if atom_index != omitted
            )
            for omitted in omission_candidates
        )
    options = set()
    found = False
    for reactant_variant in reactant_variants:
        reactant_fragment = _induced_molecule(
            reactant_mol,
            reactant_variant,
        )
        for ring_matches_ring_only in (False, True):
            result = rdFMCS.FindMCS(
                [reactant_fragment, product_fragment],
                atomCompare=rdFMCS.AtomCompare.CompareElements,
                bondCompare=rdFMCS.BondCompare.CompareAny,
                ringMatchesRingOnly=ring_matches_ring_only,
                completeRingsOnly=False,
                matchValences=False,
                timeout=2,
            )
            if result.canceled:
                return "FRAGMENT_CORRESPONDENCE_TIMEOUT", ()
            if result.numAtoms < 1:
                continue
            if (
                result.numAtoms == 1
                and len(product_atoms) > max_singleton_product_atoms
            ):
                continue
            query = compile_smarts(result.smartsString, validate=False)
            if query is None:
                return "FRAGMENT_CORRESPONDENCE_QUERY_FAILED", ()
            reactant_matches = reactant_fragment.GetSubstructMatches(
                query,
                uniquify=False,
                maxMatches=max_fragment_matches + 1,
            )
            product_matches = product_fragment.GetSubstructMatches(
                query,
                uniquify=False,
                maxMatches=max_fragment_matches + 1,
            )
            if (
                len(reactant_matches) > max_fragment_matches
                or len(product_matches) > max_fragment_matches
            ):
                return "FRAGMENT_CORRESPONDENCE_MATCH_LIMIT", ()
            for reactant_match in reactant_matches:
                original_reactant_atoms = tuple(
                    int(
                        reactant_fragment.GetAtomWithIdx(atom_index).GetIntProp(
                            "_correspondence_original_index"
                        )
                    )
                    for atom_index in reactant_match
                )
                for product_match in product_matches:
                    original_product_atoms = tuple(
                        int(
                            product_fragment.GetAtomWithIdx(atom_index).GetIntProp(
                                "_correspondence_original_index"
                            )
                        )
                        for atom_index in product_match
                    )
                    options.add(
                        tuple(
                            sorted(
                                zip(
                                    original_reactant_atoms,
                                    original_product_atoms,
                                )
                            )
                        )
                    )
                    found = True
    if not found and len(product_atoms) > max_singleton_product_atoms:
        return "FRAGMENT_CORRESPONDENCE_SINGLETON_LIMIT", ()
    if not options:
        return "FRAGMENT_CORRESPONDENCE_NOT_FOUND", ()
    return "ok", tuple(sorted(options))


def infer_fragmented_scaffold_correspondence_candidates(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_fragment_matches: int = 64,
    max_search_states: int = 5_000,
    max_candidates: int = 256,
    max_product_heavy_atoms: int = 60,
    max_unmapped_reactant_heavy_atoms: int = 4,
    max_singleton_product_atoms: int = 3,
) -> ScaffoldCorrespondenceCandidates:
    """Conservatively map a topology-changing single substrate.

    A connected MCS cannot span conserved regions separated by a deleted atom
    or a changed ring closure. This bounded fallback repeatedly removes the
    largest conserved fragment, enumerates every symmetry-equivalent match,
    and requires complete product-heavy-atom coverage. Edit-set consensus is
    applied by the caller after all minimum-edit mappings have been generated.
    """
    if any(component.atom_mapped for component in reactants + products):
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_SKIPPED_MAPPED_INPUT",), False
        )
    substantial_reactants = tuple(
        component for component in reactants if _heavy_atom_count(component) >= 2
    )
    substantial_products = tuple(
        component for component in products if _heavy_atom_count(component) >= 2
    )
    if len(substantial_reactants) != 1:
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_REQUIRES_ONE_SUBSTRATE",), False
        )
    if len(substantial_products) != 1 or any(
        _heavy_atom_count(component) > 0
        for component in products
        if component is not substantial_products[0]
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_REQUIRES_ONE_PRODUCT",), False
        )
    reactant = substantial_reactants[0]
    product = substantial_products[0]
    reactant_mol = parse_smiles(reactant.input_smiles)
    product_mol = parse_smiles(product.input_smiles)
    if reactant_mol is None or product_mol is None:
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_PARSE_FAILED",), False
        )
    reactant_atoms = tuple(
        int(atom.GetIdx())
        for atom in reactant_mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    )
    product_atoms = tuple(
        int(atom.GetIdx())
        for atom in product_mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    )
    if (
        len(product_atoms) > max_product_heavy_atoms
        or len(product_atoms) > len(reactant_atoms)
        or len(reactant_atoms) - len(product_atoms)
        > max_unmapped_reactant_heavy_atoms
        or len(product_atoms) / max(len(reactant_atoms), 1) < 0.6
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_COMPLEXITY_LIMIT",), False
        )
    reactant_elements = _element_counts((reactant,))
    product_elements = _element_counts((product,))
    if any(
        count > reactant_elements.get(element, 0)
        for element, count in product_elements.items()
    ):
        return ScaffoldCorrespondenceCandidates(
            (), ("FRAGMENT_CORRESPONDENCE_UNACCOUNTED_PRODUCT_ATOMS",), False
        )

    results = set()
    cache: Dict[
        Tuple[Tuple[int, ...], Tuple[int, ...]],
        tuple[str, Tuple[Tuple[Tuple[int, int], ...], ...]],
    ] = {}
    warnings = set()
    search_states = 0
    limit_reached = False

    def fragment_options(
        remaining_reactant: Tuple[int, ...],
        remaining_product: Tuple[int, ...],
    ) -> tuple[str, Tuple[Tuple[Tuple[int, int], ...], ...]]:
        key = (remaining_reactant, remaining_product)
        if key not in cache:
            cache[key] = _fragment_correspondence_options(
                reactant_mol,
                product_mol,
                remaining_reactant,
                remaining_product,
                max_fragment_matches=max_fragment_matches,
                max_singleton_product_atoms=max_singleton_product_atoms,
            )
        return cache[key]

    def search(
        remaining_reactant: Tuple[int, ...],
        remaining_product: Tuple[int, ...],
        selected: Tuple[Tuple[int, int], ...],
    ) -> None:
        nonlocal search_states, limit_reached
        if limit_reached:
            return
        search_states += 1
        if search_states > max_search_states:
            warnings.add("FRAGMENT_CORRESPONDENCE_SEARCH_LIMIT")
            limit_reached = True
            return
        if not remaining_product:
            results.add(tuple(sorted(selected)))
            if len(results) > max_candidates:
                warnings.add("FRAGMENT_CORRESPONDENCE_CANDIDATE_LIMIT")
                limit_reached = True
            return
        status, options = fragment_options(
            remaining_reactant,
            remaining_product,
        )
        if status != "ok":
            warnings.add(status)
            if status in {
                "FRAGMENT_CORRESPONDENCE_TIMEOUT",
                "FRAGMENT_CORRESPONDENCE_MATCH_LIMIT",
                "FRAGMENT_CORRESPONDENCE_QUERY_FAILED",
            }:
                limit_reached = True
            return
        for option in options:
            used_reactant = {left for left, _ in option}
            used_product = {right for _, right in option}
            search(
                tuple(
                    atom
                    for atom in remaining_reactant
                    if atom not in used_reactant
                ),
                tuple(
                    atom
                    for atom in remaining_product
                    if atom not in used_product
                ),
                selected + option,
            )
            if limit_reached:
                return

    search(reactant_atoms, product_atoms, ())
    if limit_reached:
        return ScaffoldCorrespondenceCandidates(
            (),
            tuple(sorted(warnings)),
            False,
        )
    if not results:
        return ScaffoldCorrespondenceCandidates(
            (),
            tuple(
                sorted(
                    warnings.union(
                        {"FRAGMENT_CORRESPONDENCE_NOT_FOUND"}
                    )
                )
            ),
            False,
        )
    candidates = tuple(
        sorted(
            tuple(
                (
                    reactant.component_index,
                    reactant_atom,
                    product.component_index,
                    product_atom,
                )
                for reactant_atom, product_atom in mapping
            )
            for mapping in results
        )
    )
    return ScaffoldCorrespondenceCandidates(candidates, (), True)


def infer_global_correspondence_candidates(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    *,
    max_component_matches: int = 32,
    max_combinations: int = 50_000,
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
        # Retain bounded one-atom-loss alternatives.  Maximum common
        # substructure matching otherwise forces a leaving atom to consume the
        # product atom that should be assigned to another partner (for example
        # choosing which alcohol oxygen survives an etherification).  These are
        # correspondence hypotheses only; minimum-edit consensus decides which
        # alternative, if any, is usable.
        if candidates:
            partial = []
            for mapping in candidates:
                if len(mapping) < 2:
                    continue
                for removed_index in range(len(mapping)):
                    partial.append(
                        mapping[:removed_index] + mapping[removed_index + 1 :]
                    )
                    if len(partial) >= max_component_matches:
                        break
                if len(partial) >= max_component_matches:
                    break
            candidates = tuple(sorted(set(candidates).union(partial)))
        # A small connected reagent may become disconnected in the product
        # (Br-Br addition is the simplest example).  MCS matching can map only
        # one resulting fragment.  Add a bounded, element-preserving atomwise
        # option for small components; edit-set minimization and consensus
        # remain responsible for accepting or rejecting the alternatives.
        reactant_mol = parse_smiles(reactant.input_smiles)
        if reactant_mol is not None and reactant_mol.GetNumHeavyAtoms() <= 4:
            from itertools import product as cartesian_product

            heavy_atoms = tuple(
                int(atom.GetIdx())
                for atom in reactant_mol.GetAtoms()
                if atom.GetAtomicNum() > 1
            )
            product_atoms_by_element: Dict[str, Tuple[int, ...]] = {}
            for atom in product_mol.GetAtoms():
                if atom.GetAtomicNum() > 1:
                    product_atoms_by_element.setdefault(
                        str(atom.GetSymbol()), ()
                    )
                    product_atoms_by_element[str(atom.GetSymbol())] += (
                        int(atom.GetIdx()),
                    )
            pools = tuple(
                product_atoms_by_element.get(
                    str(reactant_mol.GetAtomWithIdx(index).GetSymbol()),
                    (),
                )
                for index in heavy_atoms
            )
            atomwise = []
            if heavy_atoms and all(pools):
                for targets in cartesian_product(*pools):
                    if len(set(targets)) != len(targets):
                        continue
                    atomwise.append(
                        tuple(
                            (
                                reactant.component_index,
                                source,
                                product.component_index,
                                target,
                            )
                            for source, target in zip(heavy_atoms, targets)
                        )
                    )
                    if len(atomwise) >= max_component_matches:
                        break
            candidates = tuple(sorted(set(candidates).union(atomwise)))
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
    results, limit_reached = _global_correspondence_assignments(
        options,
        product_heavy_atoms,
        max_combinations=max_combinations,
        max_candidates=max_candidates,
    )
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
        fragment_options = dict(options)
        fragment_warnings = []
        for reactant in reactants:
            candidates, candidate_warnings = (
                _single_cut_fragment_correspondence_candidates(
                    reactant,
                    product,
                    max_matches=max_component_matches,
                )
            )
            fragment_warnings.extend(candidate_warnings)
            if candidates:
                fragment_options[reactant.component_index] = tuple(
                    sorted(
                        set(fragment_options.get(reactant.component_index, ())).union(
                            candidates
                        )
                    )
                )
        warnings.extend(fragment_warnings)
        results, limit_reached = _global_correspondence_assignments(
            fragment_options,
            product_heavy_atoms,
            max_combinations=max_combinations,
            max_candidates=max_candidates,
        )
        if results:
            warnings.append("INFERRED_SINGLE_CUT_FRAGMENT_CORRESPONDENCE")
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

    if any(
        _has_unsupported_mcs_bond(molecule)
        for molecule in (reactant_mol, product_mol)
    ):
        return ScaffoldCorrespondenceCandidates(
            (),
            ("SCAFFOLD_CORRESPONDENCE_UNSUPPORTED_DATIVE_BOND",),
            False,
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
    """Map a conserved scaffold while allowing unmatched connected branches.

    Unlike verified scaffold correspondence, this observation-only path does
    not require every product atom to have a reported reactant source. It is
    intentionally limited to one substantial substrate and one product; the
    caller must validate boundary topology and require one normalized
    transformation across all minimum-edit mappings.
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
    candidates, warnings = _component_correspondence_candidates(
        reactant,
        product,
        max_matches=min(max_candidates, 64),
    )
    if warnings:
        return ScaffoldCorrespondenceCandidates((), warnings, False)
    minimum_heavy_atoms = min(
        int(reactant_mol.GetNumHeavyAtoms()),
        int(product_mol.GetNumHeavyAtoms()),
    )
    conserved = tuple(
        mapping
        for mapping in candidates
        if len(mapping) >= 3
        and len(mapping) / max(minimum_heavy_atoms, 1) >= 0.5
    )
    if len(conserved) > max_candidates:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_CANDIDATE_LIMIT",), False
        )
    if not conserved:
        return ScaffoldCorrespondenceCandidates(
            (), ("PARTIAL_CORRESPONDENCE_NOT_FOUND",), False
        )
    return ScaffoldCorrespondenceCandidates(
        tuple(sorted(conserved)), (), True
    )


__all__ = [
    "AtomPair",
    "REACTION_CORRESPONDENCE_VERSION",
    "ScaffoldCorrespondenceCandidates",
    "infer_fragmented_scaffold_correspondence_candidates",
    "infer_global_correspondence_candidates",
    "infer_partial_scaffold_correspondence_candidates",
    "infer_scaffold_correspondence_candidates",
]
