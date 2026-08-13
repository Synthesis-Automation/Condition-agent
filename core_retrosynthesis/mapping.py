"""Materialize deterministic internal correspondence as mapped reaction SMILES."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

from rdkit import Chem

from reactive_taxonomy.reaction_edits import resolve_structural_evidence
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


_ATOM_MAP = re.compile(r":\d+\]")


@dataclass(frozen=True)
class MaterializedMapping:
    """Mapped reaction plus explicit correspondence provenance."""

    reaction_smiles: str
    evidence: str
    confidence: float


def _mapped_components(components: tuple, assignments: dict[tuple[int, int], int]) -> Optional[str]:
    values = []
    for component_index, component in enumerate(components):
        molecule = Chem.MolFromSmiles(component.input_smiles)
        if molecule is None:
            return None
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        for (assigned_component, atom_index), map_number in assignments.items():
            if assigned_component != component_index:
                continue
            if atom_index < 0 or atom_index >= molecule.GetNumAtoms():
                return None
            molecule.GetAtomWithIdx(atom_index).SetAtomMapNum(map_number)
        values.append(
            Chem.MolToSmiles(molecule, canonical=False, isomericSmiles=True)
        )
    return ".".join(values)


def materialize_atom_mapping(reaction_smiles: str) -> Optional[MaterializedMapping]:
    """Return supplied maps or re-run the deterministic correspondence provider."""

    if _ATOM_MAP.search(reaction_smiles):
        return MaterializedMapping(
            reaction_smiles=reaction_smiles,
            evidence="supplied_atom_mapping",
            confidence=1.0,
        )
    parsed = parse_reaction_smiles(
        reaction_smiles,
        include_molecular_interpretation=False,
    )
    if not parsed.valid or len(parsed.products) != 1:
        return None
    resolved = resolve_structural_evidence(parsed.reactants, parsed.products)
    if not resolved.valid or not resolved.atom_correspondence:
        return None
    reactant_assignments: dict[tuple[int, int], int] = {}
    product_assignments: dict[tuple[int, int], int] = {}
    for map_number, correspondence in enumerate(
        sorted(resolved.atom_correspondence),
        start=1,
    ):
        reactant_component, reactant_atom, product_component, product_atom = (
            correspondence
        )
        reactant_assignments[(reactant_component, reactant_atom)] = map_number
        product_assignments[(product_component, product_atom)] = map_number
    reactants = _mapped_components(parsed.reactants, reactant_assignments)
    products = _mapped_components(parsed.products, product_assignments)
    if reactants is None or products is None:
        return None
    product_molecule = Chem.MolFromSmiles(products)
    if product_molecule is None or any(
        atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() <= 0
        for atom in product_molecule.GetAtoms()
    ):
        return None
    return MaterializedMapping(
        reaction_smiles=f"{reactants}>>{products}",
        evidence=str(resolved.evidence),
        confidence=float(resolved.confidence),
    )


__all__ = ["MaterializedMapping", "materialize_atom_mapping"]
