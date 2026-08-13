"""Deterministic graph-derived strategic complexity for retrosynthesis.

The score measures how much a verified reverse step simplifies the target's
product-derived scaffold.  It is deliberately separate from feasibility,
conditions, reaction names, and learned synthetic-accessibility estimates.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any, Mapping, Tuple

from .chemistry.rdkit_utils import mol_to_canonical_smiles, parse_smiles


STRATEGIC_COMPLEXITY_DEFINITION_PATH = (
    Path(__file__).with_name("definitions") / "strategic_complexity.v1.json"
)
STRATEGIC_COMPLEXITY_DEFINITION_ID = "strategic_complexity.v1"
STRATEGIC_COMPLEXITY_SCHEMA_VERSION = "1.0"
_STRATEGIC_CLASSES = frozenset(
    {"scaffold_split", "ring_construction", "stereochemistry_installation"}
)
_CLASSIFICATIONS = frozenset(
    {
        *_STRATEGIC_CLASSES,
        "peripheral_attachment",
        "functional_group_interconversion",
        "protection_state_change",
        "complexity_increasing",
        "no_progress",
        "unresolved",
    }
)


@dataclass(frozen=True)
class MoleculeComplexityTrace:
    """Inspectable deterministic molecular graph-complexity components."""

    canonical_smiles: str
    heavy_atom_count: int
    heavy_bond_count: int
    cycle_rank: int
    ring_system_count: int
    largest_ring_system_atom_count: int
    fused_ring_bond_count: int
    branching_excess: int
    assigned_stereocenter_count: int
    atom_environment_count: int
    raw_complexity: float
    definition_id: str = STRATEGIC_COMPLEXITY_DEFINITION_ID
    schema_version: str = STRATEGIC_COMPLEXITY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible trace."""

        return asdict(self)


@dataclass(frozen=True)
class RetrosyntheticComplexityReduction:
    """Auditable strategic simplification produced by one reverse step."""

    score: float
    strategic_class: str
    is_strategic: bool
    evidence: str
    target: MoleculeComplexityTrace
    precursors: Tuple[MoleculeComplexityTrace, ...]
    product_derived_component_heavy_atom_counts: Tuple[int, ...]
    largest_retained_core_fraction: float
    core_fragmentation_score: float
    ring_topology_reduction_score: float
    graph_complexity_reduction_score: float
    graph_complexity_reduction_fraction: float
    stereochemical_simplification_score: float
    convergency_score: float
    tactical_penalty: float
    formed_product_bond_count: int
    formed_product_ring_bond_count: int
    product_cycle_rank_reduction: int
    extra_precursor_heavy_atom_count: int
    warnings: Tuple[str, ...]
    definition_id: str = STRATEGIC_COMPLEXITY_DEFINITION_ID
    schema_version: str = STRATEGIC_COMPLEXITY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible audit."""

        return asdict(self)


def _nonnegative_number(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field} must be a nonnegative number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise ValueError(f"{field} must be a nonnegative number")
    return normalized


def validate_strategic_complexity_definition(
    value: Mapping[str, Any],
) -> list[str]:
    """Return stable validation errors for the declarative policy."""

    errors = []
    if value.get("definition_id") != STRATEGIC_COMPLEXITY_DEFINITION_ID:
        errors.append("strategic_complexity:invalid_definition_id")
    if value.get("schema_version") != STRATEGIC_COMPLEXITY_SCHEMA_VERSION:
        errors.append("strategic_complexity:invalid_schema_version")
    molecule = value.get("molecule_complexity")
    reduction = value.get("retrosynthetic_reduction")
    if not isinstance(molecule, dict) or not isinstance(reduction, dict):
        return [*errors, "strategic_complexity:missing_sections"]
    required_molecule = {
        "heavy_atom",
        "heavy_bond",
        "cycle_rank",
        "ring_system",
        "fused_ring_bond",
        "branching",
        "assigned_stereocenter",
        "atom_environment",
    }
    if set(molecule) != required_molecule:
        errors.append("strategic_complexity:invalid_molecule_weights")
    weights = reduction.get("component_weights")
    required_weights = {
        "core_fragmentation",
        "ring_topology",
        "graph_complexity",
        "stereochemical_simplification",
        "convergency",
    }
    if not isinstance(weights, dict) or set(weights) != required_weights:
        errors.append("strategic_complexity:invalid_reduction_weights")
    else:
        try:
            normalized = [
                _nonnegative_number(raw, f"{name} weight")
                for name, raw in weights.items()
            ]
            if not math.isclose(sum(normalized), 1.0, abs_tol=1e-12):
                errors.append("strategic_complexity:weights_do_not_sum_to_one")
        except ValueError:
            errors.append("strategic_complexity:invalid_reduction_weight_value")
    penalties = reduction.get("tactical_penalties")
    if not isinstance(penalties, dict) or set(penalties) != _CLASSIFICATIONS:
        errors.append("strategic_complexity:invalid_tactical_penalties")
    for field in (
        "frontier_secondary_component_weight",
        "complexity_increase_fraction_threshold",
        "strategic_candidate_minimum_score",
    ):
        try:
            _nonnegative_number(reduction.get(field), field)
        except ValueError:
            errors.append(f"strategic_complexity:invalid_{field}")
    for field in (
        "minimum_substantive_fragment_heavy_atoms",
        "complex_target_minimum_heavy_atoms",
        "complex_target_minimum_cycle_rank",
    ):
        raw = reduction.get(field)
        if isinstance(raw, bool) or not isinstance(raw, int) or raw < 1:
            errors.append(f"strategic_complexity:invalid_{field}")
    return errors


@lru_cache(maxsize=1)
def load_strategic_complexity_definition() -> Mapping[str, Any]:
    """Load and validate the canonical strategic-complexity policy."""

    value = json.loads(
        STRATEGIC_COMPLEXITY_DEFINITION_PATH.read_text(encoding="utf-8")
    )
    errors = validate_strategic_complexity_definition(value)
    if errors:
        raise ValueError("; ".join(errors))
    return value


def _ring_systems(molecule: Any) -> tuple[frozenset[int], ...]:
    rings = [set(int(index) for index in ring) for ring in molecule.GetRingInfo().AtomRings()]
    systems: list[set[int]] = []
    for ring in rings:
        overlapping = [index for index, system in enumerate(systems) if ring & system]
        if not overlapping:
            systems.append(set(ring))
            continue
        merged = set(ring)
        for index in reversed(overlapping):
            merged.update(systems.pop(index))
        systems.append(merged)
    return tuple(sorted((frozenset(system) for system in systems), key=lambda item: tuple(sorted(item))))


def _canonical_without_maps(molecule: Any) -> str:
    from rdkit import Chem

    copied = Chem.Mol(molecule)
    for atom in copied.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol_to_canonical_smiles(copied) or ""


def _molecule_trace(molecule: Any) -> MoleculeComplexityTrace:
    from rdkit import Chem

    definition = load_strategic_complexity_definition()
    weights = definition["molecule_complexity"]
    heavy_atoms = tuple(atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1)
    heavy_indices = {int(atom.GetIdx()) for atom in heavy_atoms}
    heavy_bonds = tuple(
        bond
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in heavy_indices
        and int(bond.GetEndAtomIdx()) in heavy_indices
    )
    fragment_count = len(Chem.GetMolFrags(molecule, asMols=False))
    cycle_rank = max(0, len(heavy_bonds) - len(heavy_atoms) + fragment_count)
    systems = _ring_systems(molecule)
    fused_bonds = sum(
        molecule.GetRingInfo().NumBondRings(int(bond.GetIdx())) > 1
        for bond in heavy_bonds
    )
    branching = sum(
        max(
            0,
            sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors()) - 2,
        )
        for atom in heavy_atoms
    )
    stereocenters = len(
        Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=False,
            includeCIP=True,
        )
    )
    environments = {
        (
            int(atom.GetAtomicNum()),
            int(atom.GetFormalCharge()),
            bool(atom.GetIsAromatic()),
            sum(neighbor.GetAtomicNum() > 1 for neighbor in atom.GetNeighbors()),
            int(atom.GetTotalNumHs()),
        )
        for atom in heavy_atoms
    }
    components = {
        "heavy_atom": len(heavy_atoms),
        "heavy_bond": len(heavy_bonds),
        "cycle_rank": cycle_rank,
        "ring_system": len(systems),
        "fused_ring_bond": fused_bonds,
        "branching": branching,
        "assigned_stereocenter": stereocenters,
        "atom_environment": len(environments),
    }
    raw = sum(float(weights[name]) * value for name, value in components.items())
    return MoleculeComplexityTrace(
        canonical_smiles=_canonical_without_maps(molecule),
        heavy_atom_count=len(heavy_atoms),
        heavy_bond_count=len(heavy_bonds),
        cycle_rank=cycle_rank,
        ring_system_count=len(systems),
        largest_ring_system_atom_count=max((len(system) for system in systems), default=0),
        fused_ring_bond_count=fused_bonds,
        branching_excess=branching,
        assigned_stereocenter_count=stereocenters,
        atom_environment_count=len(environments),
        raw_complexity=round(raw, 8),
    )


def assess_molecule_complexity(smiles: str) -> MoleculeComplexityTrace:
    """Return deterministic molecular complexity for one parseable structure."""

    molecule = parse_smiles(smiles)
    if molecule is None:
        raise ValueError("invalid molecule for strategic-complexity assessment")
    return _molecule_trace(molecule)


def _mapped_bonds(molecules: tuple[Any, ...]) -> dict[tuple[int, int], tuple[float, bool, tuple[str, str]]]:
    values = {}
    for molecule in molecules:
        for bond in molecule.GetBonds():
            left = bond.GetBeginAtom()
            right = bond.GetEndAtom()
            left_map = int(left.GetAtomMapNum())
            right_map = int(right.GetAtomMapNum())
            if left_map <= 0 or right_map <= 0:
                continue
            key = tuple(sorted((left_map, right_map)))
            values[key] = (
                float(bond.GetBondTypeAsDouble()),
                bool(bond.IsInRing()),
                (left.GetSymbol(), right.GetSymbol()),
            )
    return values


def _induced_cycle_rank(molecule: Any, retained_maps: frozenset[int]) -> int:
    selected = {
        int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and int(atom.GetAtomMapNum()) in retained_maps
    }
    if not selected:
        return 0
    edges = [
        bond
        for bond in molecule.GetBonds()
        if int(bond.GetBeginAtomIdx()) in selected
        and int(bond.GetEndAtomIdx()) in selected
    ]
    adjacency = {index: set() for index in selected}
    for bond in edges:
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        adjacency[left].add(right)
        adjacency[right].add(left)
    components = 0
    remaining = set(selected)
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            current = stack.pop()
            neighbors = adjacency[current] & remaining
            remaining.difference_update(neighbors)
            stack.extend(neighbors)
    return max(0, len(edges) - len(selected) + components)


def _reaction_parts(reaction_smiles: str) -> tuple[tuple[Any, ...], tuple[Any, ...]]:
    sections = str(reaction_smiles or "").split(">")
    if len(sections) != 3 or not sections[0].strip() or not sections[2].strip():
        raise ValueError("strategic-complexity assessment requires a reaction SMILES")
    precursors = tuple(
        molecule
        for component in sections[0].split(".")
        if component.strip()
        for molecule in (parse_smiles(component),)
        if molecule is not None
    )
    products = tuple(
        molecule
        for component in sections[2].split(".")
        if component.strip()
        for molecule in (parse_smiles(component),)
        if molecule is not None
    )
    if not precursors or not products:
        raise ValueError("invalid reaction for strategic-complexity assessment")
    return precursors, products


def assess_retrosynthetic_complexity_reduction(
    reaction_smiles: str,
) -> RetrosyntheticComplexityReduction:
    """Assess strategic target simplification from a forward reaction graph."""

    definition = load_strategic_complexity_definition()
    policy = definition["retrosynthetic_reduction"]
    precursors, products = _reaction_parts(reaction_smiles)
    target_molecule = max(products, key=lambda molecule: molecule.GetNumHeavyAtoms())
    target = _molecule_trace(target_molecule)
    precursor_traces = tuple(
        sorted(
            (_molecule_trace(molecule) for molecule in precursors),
            key=lambda trace: (-trace.raw_complexity, trace.canonical_smiles),
        )
    )
    product_maps = {
        int(atom.GetAtomMapNum())
        for atom in target_molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and int(atom.GetAtomMapNum()) > 0
    }
    mapped = len(product_maps) == target.heavy_atom_count
    retained_counts = []
    precursor_product_cycle_rank = 0
    precursor_stereo_maps = set()
    extra_heavy_atoms = 0
    from rdkit import Chem

    for molecule in precursors:
        retained = frozenset(
            int(atom.GetAtomMapNum())
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() > 1
            and int(atom.GetAtomMapNum()) in product_maps
        )
        if retained:
            retained_counts.append(len(retained))
            precursor_product_cycle_rank += _induced_cycle_rank(molecule, retained)
        extra_heavy_atoms += sum(
            atom.GetAtomicNum() > 1
            and int(atom.GetAtomMapNum()) not in product_maps
            for atom in molecule.GetAtoms()
        )
        for atom_index, _ in Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=False,
            includeCIP=True,
        ):
            atom_map = int(molecule.GetAtomWithIdx(atom_index).GetAtomMapNum())
            if atom_map in product_maps:
                precursor_stereo_maps.add(atom_map)
    retained_counts = sorted(retained_counts, reverse=True) if mapped else []
    largest_fraction = (
        retained_counts[0] / target.heavy_atom_count
        if retained_counts and target.heavy_atom_count
        else 1.0
    )
    core_fragmentation = (
        min(1.0, 2.0 * (1.0 - largest_fraction))
        if len(retained_counts) >= 2
        else 0.0
    )
    minimum_fragment = int(policy["minimum_substantive_fragment_heavy_atoms"])
    convergency = (
        core_fragmentation
        if len(retained_counts) >= 2 and retained_counts[1] >= minimum_fragment
        else 0.0
    )
    precursor_bonds = _mapped_bonds(precursors)
    product_bonds = _mapped_bonds((target_molecule,))
    formed = {
        key: value
        for key, value in product_bonds.items()
        if key not in precursor_bonds
    }
    changed_orders = {
        key
        for key, value in product_bonds.items()
        if key in precursor_bonds
        and not math.isclose(value[0], precursor_bonds[key][0])
    }
    formed_ring_count = sum(value[1] for value in formed.values())
    cycle_reduction = max(0, target.cycle_rank - precursor_product_cycle_rank)
    ring_topology = (
        1.0
        if formed_ring_count
        else (
            min(1.0, cycle_reduction / max(1, target.cycle_rank))
            if mapped
            else 0.0
        )
    )
    target_stereo_maps = {
        int(target_molecule.GetAtomWithIdx(index).GetAtomMapNum())
        for index, _ in Chem.FindMolChiralCenters(
            target_molecule,
            includeUnassigned=False,
            includeCIP=True,
        )
    }
    stereo_gain = (
        len(target_stereo_maps - precursor_stereo_maps)
        / max(1, len(target_stereo_maps))
        if mapped and target_stereo_maps
        else 0.0
    )
    secondary_weight = float(policy["frontier_secondary_component_weight"])
    complexities = sorted(
        (trace.raw_complexity for trace in precursor_traces),
        reverse=True,
    )
    frontier = (
        complexities[0] + secondary_weight * sum(complexities[1:])
        if complexities
        else target.raw_complexity
    )
    complexity_fraction = (
        (target.raw_complexity - frontier) / max(target.raw_complexity, frontier)
        if max(target.raw_complexity, frontier) > 0.0
        else 0.0
    )
    graph_reduction = max(0.0, complexity_fraction)
    small_attachment = (
        len(retained_counts) >= 2 and retained_counts[1] < minimum_fragment
    )
    formed_has_heteroatom = any(
        any(element not in {"C", "H"} for element in value[2])
        for value in formed.values()
    )
    increase_threshold = float(policy["complexity_increase_fraction_threshold"])
    if ring_topology > 0.0 and formed_ring_count:
        strategic_class = "ring_construction"
    elif convergency > 0.0:
        strategic_class = "scaffold_split"
    elif stereo_gain > 0.0:
        strategic_class = "stereochemistry_installation"
    elif complexity_fraction < -increase_threshold:
        strategic_class = "complexity_increasing"
    elif not formed and extra_heavy_atoms > 0:
        strategic_class = "protection_state_change"
    elif not formed and not changed_orders:
        strategic_class = "no_progress"
    elif not formed or (small_attachment and formed_has_heteroatom):
        strategic_class = "functional_group_interconversion"
    elif small_attachment:
        strategic_class = "peripheral_attachment"
    else:
        strategic_class = "unresolved"
    weights = policy["component_weights"]
    components = {
        "core_fragmentation": core_fragmentation,
        "ring_topology": ring_topology,
        "graph_complexity": graph_reduction,
        "stereochemical_simplification": stereo_gain,
        "convergency": convergency,
    }
    tactical_penalty = float(policy["tactical_penalties"][strategic_class])
    raw_score = sum(float(weights[name]) * value for name, value in components.items())
    score = round(min(1.0, max(0.0, raw_score - tactical_penalty)), 8)
    strategic_threshold = float(policy["strategic_candidate_minimum_score"])
    warnings = []
    if not mapped:
        warnings.append("STRATEGIC_COMPLEXITY_MAPPING_UNAVAILABLE")
    if complexity_fraction < -increase_threshold:
        warnings.append("RETROSYNTHETIC_COMPLEXITY_INCREASE")
    if strategic_class in {
        "functional_group_interconversion",
        "protection_state_change",
        "peripheral_attachment",
        "no_progress",
    }:
        warnings.append("TACTICAL_NON_SCAFFOLD_TRANSFORMATION")
    return RetrosyntheticComplexityReduction(
        score=score,
        strategic_class=strategic_class,
        is_strategic=(
            mapped
            and strategic_class in _STRATEGIC_CLASSES
            and score >= strategic_threshold
        ),
        evidence="mapped_atom_correspondence" if mapped else "molecular_graph_only",
        target=target,
        precursors=precursor_traces,
        product_derived_component_heavy_atom_counts=tuple(retained_counts),
        largest_retained_core_fraction=round(largest_fraction, 8),
        core_fragmentation_score=round(core_fragmentation, 8),
        ring_topology_reduction_score=round(ring_topology, 8),
        graph_complexity_reduction_score=round(graph_reduction, 8),
        graph_complexity_reduction_fraction=round(complexity_fraction, 8),
        stereochemical_simplification_score=round(stereo_gain, 8),
        convergency_score=round(convergency, 8),
        tactical_penalty=round(tactical_penalty, 8),
        formed_product_bond_count=len(formed),
        formed_product_ring_bond_count=formed_ring_count,
        product_cycle_rank_reduction=cycle_reduction,
        extra_precursor_heavy_atom_count=extra_heavy_atoms,
        warnings=tuple(warnings),
    )


def complex_target_requires_strategic_candidate(
    assessment: RetrosyntheticComplexityReduction,
) -> bool:
    """Return whether a target is complex enough to report coverage absence."""

    policy = load_strategic_complexity_definition()["retrosynthetic_reduction"]
    return (
        assessment.target.heavy_atom_count
        >= int(policy["complex_target_minimum_heavy_atoms"])
        or assessment.target.cycle_rank
        >= int(policy["complex_target_minimum_cycle_rank"])
    )


__all__ = [
    "STRATEGIC_COMPLEXITY_DEFINITION_ID",
    "STRATEGIC_COMPLEXITY_DEFINITION_PATH",
    "STRATEGIC_COMPLEXITY_SCHEMA_VERSION",
    "MoleculeComplexityTrace",
    "RetrosyntheticComplexityReduction",
    "assess_molecule_complexity",
    "assess_retrosynthetic_complexity_reduction",
    "complex_target_requires_strategic_candidate",
    "load_strategic_complexity_definition",
    "validate_strategic_complexity_definition",
]
