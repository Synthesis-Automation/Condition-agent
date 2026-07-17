"""Conservative graph correspondence for unmapped scaffold transformations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

from .chemistry.rdkit_utils import parse_smiles
from .chemistry.smarts_cache import compile_smarts
from .reaction_models import ReactionComponent


AtomPair = Tuple[int, int, int, int]


@dataclass(frozen=True)
class ScaffoldCorrespondenceCandidates:
    """Candidate heavy-atom correspondences for one conserved scaffold."""

    candidates: Tuple[Tuple[AtomPair, ...], ...]
    warnings: Tuple[str, ...] = ()
    valid: bool = False


def _heavy_atom_count(component: ReactionComponent) -> int:
    mol = parse_smiles(component.input_smiles)
    return int(mol.GetNumHeavyAtoms()) if mol is not None else 0


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
        timeout=2,
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


__all__ = [
    "AtomPair",
    "ScaffoldCorrespondenceCandidates",
    "infer_scaffold_correspondence_candidates",
]
