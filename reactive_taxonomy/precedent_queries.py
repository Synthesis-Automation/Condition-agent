"""Deterministic nested scaffold queries for automatic precedent discovery.

Queries are graph hypotheses, not precursor molecules. They retain elements,
aromaticity, charge, isotope, bond order and selected stereo environments while
leaving nonstereo hydrogen counts and omitted substituents unconstrained.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

from rdkit import Chem

from .chemistry.smarts_cache import compile_smarts
from .fragment_search import compile_fragment_query, indexed_product, validate_fragment_target
from .search_fragments import _ring_systems, _neighborhood


_ATOM_QUERY = "[{isotope}#{element};{aromaticity};{charge:+d}]"


def discovery_policy() -> dict[str, Any]:
    """Load the versioned core-selection and bounded search policy."""
    value = json.loads((Path(__file__).parent / "definitions/precedent_discovery.v1.json").read_text("utf-8"))
    if (value.get("schema_version") != "precedent_discovery_policy.v1"
            or value.get("definition_version") != "precedent_discovery.v1@1.1"
            or value.get("stereo_hydrogen_policy") != "preserve_at_specified_centers"
            or value.get("levels") != ["core", "multiple_bond_context", "neighbor_context", "target_context"]
            or value.get("relaxations") != ["omitted_peripheral_atoms", "unconstrained_hydrogen_count", "additional_ring_fusion"]
            or value.get("core_priority") != ["heteroatoms", "atom_count", "ring_junctions"]
            or value.get("relationship_order") != ["constructed", "modified", "boundary_changed", "carried_through", "unresolved"]):
        raise ValueError("Invalid precedent discovery policy")
    for key in ("max_cores", "max_steps_per_core", "min_core_atoms", "max_target_atoms",
                "refine_above_products", "sufficient_construction_references",
                "default_timeout_seconds", "max_timeout_seconds"):
        if type(value.get(key)) is not int or value[key] < 1:
            raise ValueError(f"Invalid discovery policy: {key}")
    if value["default_timeout_seconds"] > value["max_timeout_seconds"] or value["max_steps_per_core"] > 4:
        raise ValueError("Invalid discovery policy bounds")
    return value


@dataclass(frozen=True)
class CoreQuery:
    """Query atom order is explicitly linked to the canonical target."""

    core_id: str
    level: str
    query: str
    target_atom_ids: tuple[int, ...]
    core_atom_ids: tuple[int, ...]
    query_format: str = "smarts"
    topology: str = "subgraph"


def _complete_context(mol: Any, selected: set[int], rings: tuple[frozenset[int], ...]) -> set[int]:
    """Keep the neighborhoods needed for specified stereo in selected ring subsets."""
    selected = set(selected)
    while True:
        before = set(selected)
        for i in tuple(selected):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                selected.update(n.GetIdx() for n in atom.GetNeighbors())
            for bond in atom.GetBonds():
                if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                    selected.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), *bond.GetStereoAtoms()))
        if selected == before:
            return selected


def _query(mol: Any, selected: set[int]) -> tuple[str, tuple[int, ...]]:
    editable = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        isotope = str(atom.GetIsotope()) if atom.GetIsotope() else ""
        aromaticity = "a" if atom.GetIsAromatic() else "A"
        expression = _ATOM_QUERY.format(isotope=isotope, element=atom.GetAtomicNum(),
                                        aromaticity=aromaticity, charge=atom.GetFormalCharge())
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            # A virtual H is part of tetrahedral neighbor ordering. Dropping it
            # can invert a ring-junction query on SMARTS serialization/parsing.
            expression = expression[:-1] + f";H{atom.GetTotalNumHs(includeNeighbors=True)}]"
        # Cached patterns are shared; never mutate their atoms' chiral tags.
        pattern = Chem.Mol(compile_smarts(expression, validate=True))
        replacement = pattern.GetAtomWithIdx(0)
        replacement.SetChiralTag(atom.GetChiralTag())
        editable.ReplaceAtom(atom.GetIdx(), replacement)
    for i in reversed(range(mol.GetNumAtoms())):
        if i not in selected:
            editable.RemoveAtom(i)
    query = editable.GetMol()
    Chem.GetSymmSSSR(query)
    expression = Chem.MolToSmarts(query, isomericSmiles=True)
    order = json.loads(query.GetProp("_smilesAtomOutputOrder"))
    original = sorted(selected)
    return expression, tuple(original[i] for i in order)


def plan_precedent_queries(target_smiles: str) -> dict[str, Any]:
    """Generate bounded core-first ladders; every query must match the target."""
    policy = discovery_policy()
    if not isinstance(target_smiles, str) or not target_smiles.strip() or len(target_smiles) > 5000:
        raise ValueError("Provide one target SMILES of at most 5000 characters")
    canonical, _ = indexed_product(target_smiles)
    mol = Chem.MolFromSmiles(canonical)
    if (len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() > policy["max_target_atoms"]
            or any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms())):
        raise ValueError("Provide one connected, closed-shell target without wildcard atoms")
    rings = _ring_systems(mol)
    individual = [set(r) for r in mol.GetRingInfo().AtomRings()]
    pairs = [a | b for i, a in enumerate(individual) for b in individual[i + 1:] if a & b]
    cores = [set(ring) for ring in (*pairs, *rings, *individual) if len(ring) >= policy["min_core_atoms"]]
    if not cores:
        cores = [_neighborhood(mol, {a.GetIdx()}, 1) for a in mol.GetAtoms()
                 if a.GetAtomicNum() not in (1, 6)]
        cores = [c for c in cores if len(c) >= min(policy["min_core_atoms"], mol.GetNumAtoms())]
    if not cores:
        cores = [set(range(mol.GetNumAtoms()))]

    def priority(core: set[int]) -> tuple[int, int, int, tuple[int, ...]]:
        return (-sum(mol.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6) for i in core),
                len(core), -sum(sum(b.IsInRing() for b in mol.GetAtomWithIdx(i).GetBonds()) > 2 for i in core),
                tuple(sorted(core)))

    unique = {tuple(sorted(_complete_context(mol, c, rings))) for c in cores}
    ordered = sorted((set(c) for c in unique), key=priority)
    ladders = []
    for number, core in enumerate(ordered[:policy["max_cores"]]):
        multiple = core | {b.GetOtherAtomIdx(i) for i in core for b in mol.GetAtomWithIdx(i).GetBonds()
                           if b.GetBondTypeAsDouble() > 1}
        neighbor = _neighborhood(mol, multiple, 1)
        levels = (core, multiple, neighbor, set(range(mol.GetNumAtoms())))
        steps, seen = [], set()
        for label, selected in zip(policy["levels"], levels):
            selected = _complete_context(mol, selected, rings)
            # A full target can include separate ring systems connected by a linker;
            # each ladder remains nested and connected as context grows.
            expression, order = _query(mol, selected)
            if expression in seen:
                continue
            compiled = compile_fragment_query(expression, "smarts", "subgraph")
            if not validate_fragment_target(compiled, canonical).matches_target:
                raise ValueError("Generated scaffold query lost target correspondence")
            steps.append(asdict(CoreQuery(f"core-{number + 1}", label, expression, order, tuple(sorted(core)))))
            seen.add(expression)
        ladders.append(steps[:policy["max_steps_per_core"]])
    return {"schema_version": "precedent_query_plan.v1", "definition_version": policy["definition_version"],
            "target_smiles": canonical, "ladders": ladders, "cores_truncated": len(ordered) > len(ladders),
            "relaxations": policy["relaxations"], "stereo_hydrogen_policy": policy["stereo_hydrogen_policy"],
            "target_atom_count": mol.GetNumAtoms()}
