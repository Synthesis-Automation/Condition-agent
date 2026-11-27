"""Compare deposited mapped and unmapped reaction structures without inferring maps."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache

from rdkit import Chem

from .rdkit_utils import parse_smiles
from .reaction_smiles import split_reaction_fields


SOURCE_MAPPING_SCREEN_VERSION = "source_mapping_screen.v1"
_SOURCE_MAPPING_CONFLICTS = frozenset({
    "source_mapping_unavailable_or_unreadable", "source_mapped_inputs_disagree",
    "source_mapped_products_disagree", "source_duplicate_atom_maps",
    "source_product_atom_origin_unverified",
})


def source_mapping_conflicts(reasons: Iterable[str]) -> tuple[str, ...]:
    """Return source contradictions and missing map evidence that preclude transfer."""
    return tuple(sorted(_SOURCE_MAPPING_CONFLICTS.intersection(reasons)))


@dataclass(frozen=True)
class SourceMappingEvidence:
    """Map-free source agreement and mapped product participation by input graph."""

    issues: tuple[str, ...]
    contributor_graphs: tuple[str, ...] = ()
    noncontributor_graphs: tuple[str, ...] = ()
    agent_graphs: tuple[str, ...] = ()
    definition_version: str = "source_mapping_evidence.v1"


@dataclass(frozen=True)
class _Fragment:
    graph: str
    atoms: tuple[tuple[int, int, int], ...]


def mapped_atom_identities(molecule: Chem.Mol) -> tuple[tuple[int, int, int], ...]:
    """Return map number, atomic number and isotope for each observed atom."""
    return tuple((a.GetAtomMapNum(), a.GetAtomicNum(), a.GetIsotope()) for a in molecule.GetAtoms())


def product_atom_origins_agree(
    reactants: Iterable[tuple[int, int, int]], products: Iterable[tuple[int, int, int]],
) -> bool:
    """Check observed nuclei for every product heavy atom and every explicitly mapped atom."""
    origins = {number: (element, isotope) for number, element, isotope in reactants if number}
    return all(number and origins.get(number) == (element, isotope)
               for number, element, isotope in products if element > 1 or number)


@lru_cache(maxsize=8192)
def _fragments(value: str) -> tuple[_Fragment, ...] | None:
    if not value:
        return ()
    molecule = parse_smiles(value)
    if molecule is None:
        return None
    result = []
    for fragment in Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=False):
        atoms = mapped_atom_identities(fragment)
        for atom in fragment.GetAtoms():
            atom.SetAtomMapNum(0)
        result.append(_Fragment(Chem.MolToSmiles(fragment), atoms))
    return tuple(result)


def source_material_graphs(value: str) -> tuple[str, ...]:
    """Return the full multiset of canonical constituent graphs for one claim."""
    return tuple(sorted(f.graph for f in _fragments(value) or ()))


def _reaction_sides(value: str) -> tuple | None:
    fields = split_reaction_fields(value)
    if len(fields) != 3 or not fields[0] or not fields[2]:
        return None
    sides = tuple(_fragments(side) for side in fields)
    return None if any(side is None for side in sides) else sides


def assess_source_mapping(unmapped: str, mapped: str) -> SourceMappingEvidence:
    """Check both graph multisets, unique maps and product heavy-atom origins."""
    original = _reaction_sides(unmapped)
    supplied = _reaction_sides(mapped)
    if not original or not supplied:
        return SourceMappingEvidence(("source_mapping_unavailable_or_unreadable",))
    raw_left, raw_agents, raw_right = original
    left, agents, right = supplied
    issues = []
    if Counter(f.graph for f in raw_left + raw_agents) != Counter(f.graph for f in left + agents):
        issues.append("source_mapped_inputs_disagree")
    if Counter(f.graph for f in raw_right) != Counter(f.graph for f in right):
        issues.append("source_mapped_products_disagree")
    left_atoms = [atom for f in left for atom in f.atoms]
    product_atoms = [atom for f in right for atom in f.atoms]
    right_maps = {number for number, _, _ in product_atoms if number}
    all_inputs = left_atoms + [atom for f in agents for atom in f.atoms]
    for atoms in (all_inputs, product_atoms):
        numbers = [n for n, _, _ in atoms if n]
        if len(numbers) != len(set(numbers)):
            issues.append("source_duplicate_atom_maps")
            break
    if not product_atom_origins_agree(left_atoms, product_atoms):
        issues.append("source_product_atom_origin_unverified")
    contributors = []
    noncontributors = []
    for f in left:
        bucket = contributors if any(n in right_maps for n, _, _ in f.atoms) else noncontributors
        bucket.append(f.graph)
    return SourceMappingEvidence(
        tuple(issues), tuple(sorted(contributors)),
        tuple(sorted(noncontributors)), tuple(sorted(f.graph for f in agents)),
    )
