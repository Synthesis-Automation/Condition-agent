"""Trace assessed bond formations through intermediates into the final target."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Mapping

from rdkit import Chem

from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields

from .chemistry import canonical_smiles, digest


def _marked_bonds(molecule: Chem.Mol) -> dict[int, str]:
    """Identify each bond by the isomeric graph with both endpoints marked.

    Symmetry-equivalent positions yield the same label. Marking the pair
    together preserves its relative position, including on a symmetric ring.
    """
    graph = Chem.Mol(molecule)
    for atom in graph.GetAtoms():
        atom.SetAtomMapNum(0)
    labels = {}
    for bond in graph.GetBonds():
        left, right = bond.GetBeginAtom(), bond.GetEndAtom()
        left.SetAtomMapNum(1)
        right.SetAtomMapNum(1)
        labels[bond.GetIdx()] = Chem.MolToSmiles(graph, isomericSmiles=True)
        left.SetAtomMapNum(0)
        right.SetAtomMapNum(0)
    return labels


def _atom_maps(bond: Chem.Bond) -> tuple[int, int]:
    return tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())))


@dataclass
class _StepBonds:
    product: str
    formed: set[str]
    # For each bond on a precursor, a product bond label or None if its atoms separate.
    retained: dict[tuple[str, str], set[str | None]]


def _assessed_bonds(assessment: Mapping[str, Any]) -> _StepBonds | None:
    """Read only correspondences and formed edits verified by the saved assessment."""
    passed = {g["gate_id"] for g in assessment.get("gates", ()) if g["status"] == "pass"}
    required = {
        "input_structure", "reaction_side_consistency", "atom_correspondence",
        "reaction_completeness", "reaction_signature", "stereochemistry",
    }
    signature = assessment.get("reaction_signature") or {}
    reaction = assessment.get("normalized_mapped_reaction_smiles") or ""
    fields = split_reaction_fields(reaction)
    if not required <= passed or "edits" not in signature or len(fields) != 3:
        return None
    left, _, right = fields
    product = assessment["canonical_target_smiles"]
    if canonical_smiles(left) != assessment["canonical_precursor_smiles"] or canonical_smiles(right) != product:
        return None
    product_mol = Chem.MolFromSmiles(right)
    reactants = Chem.MolFromSmiles(left)
    precursor_maps = {atom.GetAtomMapNum() for atom in reactants.GetAtoms()}
    before = {_atom_maps(bond) for bond in reactants.GetBonds()}
    labels = _marked_bonds(product_mol)
    after = {_atom_maps(bond): labels[bond.GetIdx()] for bond in product_mol.GetBonds()}
    formed = set()
    for edit in signature["edits"]:
        if edit["edit_type"] == "formed":
            numbers = [(edit.get(key) or {}).get("atom_map_number") for key in ("atom_1", "atom_2")]
            if any(type(n) is not int or n <= 0 for n in numbers):
                return None
            pair = tuple(sorted(numbers))
            if pair not in after or pair in before or not precursor_maps.issuperset(pair):
                return None
            formed.add(after[pair])

    retained = defaultdict(set)
    for precursor in Chem.GetMolFrags(reactants, asMols=True):
        smiles = canonical_smiles(Chem.MolToSmiles(precursor))
        labels = _marked_bonds(precursor)
        for bond in precursor.GetBonds():
            retained[smiles, labels[bond.GetIdx()]].add(after.get(_atom_maps(bond)))
    return _StepBonds(product, formed, dict(retained))


def construction_inventory(assessment: Mapping[str, Any]) -> dict[str, Any]:
    """List the proposed formed bonds that persist in the canonical target.

    Follow validated step mappings through all indistinguishable intermediate
    sites. Keep an unresolved step if the possibilities reach different target
    sites, include a lost bond, or encounter a gap in the recorded mapping.
    Indices are zero-based in the canonical target SMILES. Equal bond IDs
    identify one target site up to stereochemistry-preserving graph symmetry.
    """
    steps = assessment.get("step_assessments", ())
    step_ids = {item["external_step_id"] for item in steps}
    inventory = {
        "key_constructions": [], "unresolved_construction_step_ids": sorted(step_ids),
        "leaf_smiles": list(assessment.get("leaf_smiles", ())),
    }
    topology = {g["gate_id"] for g in assessment.get("topology_gates", ()) if g["status"] == "pass"}
    if not {"route_target", "step_identity", "product_ownership", "route_connectivity"} <= topology:
        return inventory
    target = assessment["canonical_target_smiles"]
    molecule = Chem.MolFromSmiles(target)
    if molecule is None:
        return inventory
    target_labels = _marked_bonds(molecule)
    sites = defaultdict(list)
    for bond in molecule.GetBonds():
        sites[target_labels[bond.GetIdx()]].append(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))

    by_step = {}
    parents = defaultdict(list)
    for item in steps:
        step = item["assessment"]
        mapped = _assessed_bonds(step)
        by_step[item["external_step_id"]] = mapped
        for precursor in set((step.get("canonical_precursor_smiles") or "").split(".")):
            parents[precursor].append(mapped)

    # None denotes a known lost bond. An unknown correspondence stays distinct.
    unknown = "unresolved"

    @lru_cache(maxsize=None)
    def trace(product: str, bond: str) -> frozenset[str | None]:
        if product == target:
            return frozenset((bond if bond in sites else unknown,))
        possibilities = set()
        for parent in parents[product]:
            if parent is None or (product, bond) not in parent.retained:
                possibilities.add(unknown)
                continue
            for downstream in parent.retained[product, bond]:
                if downstream is None:
                    possibilities.add(None)
                else:
                    possibilities.update(trace(parent.product, downstream))
        return frozenset(possibilities or {unknown})

    unresolved = set()
    for step_id, step in sorted(by_step.items()):
        if step is None:
            unresolved.add(step_id)
            continue
        projected = set()
        for bond in step.formed:
            destinations = trace(step.product, bond)
            if destinations == {None}:
                continue
            if len(destinations) != 1 or unknown in destinations:
                unresolved.add(step_id)
                continue
            projected.update(destinations)
        for label in sorted(projected):
            pairs = sorted(sites[label])
            i, j = pairs[0]
            inventory["key_constructions"].append({
                "step_id": step_id, "bond_id": digest("TARGETBOND1", target, label),
                "target_atom_indices": [i, j], "equivalent_target_atom_pairs": pairs,
                "elements": [molecule.GetAtomWithIdx(n).GetSymbol() for n in (i, j)],
                "bond_type": str(molecule.GetBondBetweenAtoms(i, j).GetBondType()),
            })
    inventory["unresolved_construction_step_ids"] = sorted(unresolved)
    return inventory
