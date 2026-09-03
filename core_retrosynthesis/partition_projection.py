"""Project mapped reactions and route frontiers into target partitions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from rdkit import Chem

from .chemistry import canonical_smiles, digest
from .mapping import materialize_atom_mapping
from .route_contract import MoleculeOccurrenceNode, ReactionRouteTree
from .synthetic_partition import (
    InterfaceHypothesis,
    LatentModuleState,
    NonTargetAtom,
    PartitionSearchDiagnostics,
    SyntheticPartition,
    SyntheticPartitionLandscape,
    analyze_partition_target,
    create_synthetic_partition,
)


PARTITION_PROJECTION_SCHEMA_VERSION = "1.0"
PARTITION_PROJECTION_ALGORITHM_VERSION = "route_partition_projection.v1"
MAXIMUM_ISOMORPHISMS = 33


@dataclass(frozen=True)
class ReactionPartitionProjection:
    """Target ownership induced by one mapped reverse transformation."""

    target_smiles: str
    mapped_reaction_smiles: str
    module_atom_sets: tuple[tuple[int, ...], ...]
    target_boundary_bonds: tuple[tuple[int, int, str], ...]
    precursor_mapped_smiles: tuple[str, ...]
    mapping_evidence: str
    mapping_confidence: float
    warnings: tuple[str, ...] = ()

    @property
    def k(self) -> int:
        """Return the number of target-contributing precursor components."""

        return len(self.module_atom_sets)


@dataclass(frozen=True)
class RouteFrontierPartition:
    """One depth-bounded route frontier and its target partition."""

    frontier_id: str
    depth: int
    component_occurrence_ids: tuple[str, ...]
    latent_states: tuple[LatentModuleState, ...]
    partition: SyntheticPartition
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible frontier projection."""

        return {
            "frontier_id": self.frontier_id,
            "depth": self.depth,
            "component_occurrence_ids": list(self.component_occurrence_ids),
            "latent_states": [state.to_dict() for state in self.latent_states],
            "partition": self.partition.to_dict(),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class RoutePartitionProjection:
    """All completely mapped level frontiers of one concrete route tree."""

    source_tree_id: str
    route_kind: str
    target_smiles: str
    frontiers: tuple[RouteFrontierPartition, ...]
    unresolved_occurrence_ids: tuple[str, ...]
    warnings: tuple[str, ...]
    algorithm_version: str = PARTITION_PROJECTION_ALGORITHM_VERSION
    schema_version: str = PARTITION_PROJECTION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible route projection."""

        return {
            "source_tree_id": self.source_tree_id,
            "route_kind": self.route_kind,
            "target_smiles": self.target_smiles,
            "frontiers": [frontier.to_dict() for frontier in self.frontiers],
            "unresolved_occurrence_ids": list(self.unresolved_occurrence_ids),
            "warnings": list(self.warnings),
            "algorithm_version": self.algorithm_version,
            "schema_version": self.schema_version,
        }

    def to_landscape(self) -> SyntheticPartitionLandscape:
        """Return a deduplicated role-neutral landscape from route frontiers."""

        canonical, target_id, target_atoms, _ = analyze_partition_target(
            self.target_smiles
        )
        by_id: dict[str, SyntheticPartition] = {}
        for frontier in self.frontiers:
            by_id.setdefault(frontier.partition.partition_id, frontier.partition)
        partitions = tuple(
            sorted(by_id.values(), key=lambda item: (item.k, item.partition_id))
        )
        return SyntheticPartitionLandscape(
            target_id=target_id,
            target_smiles=canonical,
            target_atoms=target_atoms,
            partitions=partitions,
            searched_k_values=tuple(sorted({item.k for item in partitions})),
            generation_coverage=("route_frontier_projection",),
            unresolved_motifs=(),
            abstained=False,
            abstention_reasons=(),
            diagnostics=PartitionSearchDiagnostics(
                returned_partition_count=len(partitions)
            ),
            warnings=self.warnings,
        )


@dataclass(frozen=True)
class _ComponentProjection:
    mapped_smiles: str
    canonical_smiles: str
    atom_ownership: tuple[tuple[int, int], ...]

    @property
    def target_atom_maps(self) -> tuple[int, ...]:
        return tuple(sorted({target_map for _, target_map in self.atom_ownership}))


@dataclass(frozen=True)
class _ProjectedOccurrence:
    node: MoleculeOccurrenceNode
    target_atom_maps: tuple[int, ...]
    atom_ownership: tuple[tuple[int, int], ...]
    children: tuple["_ProjectedOccurrence", ...]
    expansion_resolved: bool
    warnings: tuple[str, ...]


def _unmapped_copy(molecule: Chem.Mol) -> Chem.Mol:
    value = Chem.Mol(molecule)
    for atom in value.GetAtoms():
        atom.SetAtomMapNum(0)
    return value


def _canonical_unmapped(molecule: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        _unmapped_copy(molecule),
        canonical=True,
        isomericSmiles=True,
    )


def _full_isomorphisms(
    query: Chem.Mol,
    target: Chem.Mol,
) -> tuple[tuple[int, ...], ...]:
    if query.GetNumAtoms() != target.GetNumAtoms():
        return ()
    matches = target.GetSubstructMatches(
        _unmapped_copy(query),
        uniquify=False,
        useChirality=True,
        maxMatches=MAXIMUM_ISOMORPHISMS,
    )
    return tuple(sorted({tuple(int(index) for index in match) for match in matches}))


def _split_mapped_components(side: str) -> tuple[Chem.Mol, ...]:
    molecule = Chem.MolFromSmiles(side)
    if molecule is None:
        return ()
    return tuple(Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True))


def _propagate_reaction_ownership(
    product_smiles: str,
    reaction_smiles: str,
    product_atom_ownership: Mapping[int, int],
) -> tuple[
    tuple[_ComponentProjection, ...],
    str,
    float,
    tuple[str, ...],
]:
    materialized = materialize_atom_mapping(reaction_smiles)
    if materialized is None or materialized.reaction_smiles.count(">>") != 1:
        raise ValueError("reaction mapping could not be materialized")
    reactant_side, mapped_product_smiles = materialized.reaction_smiles.split(">>")
    product = Chem.MolFromSmiles(mapped_product_smiles)
    node_product = Chem.MolFromSmiles(product_smiles)
    if product is None or node_product is None:
        raise ValueError("reaction product could not be parsed")
    product_matches = _full_isomorphisms(node_product, product)
    if not product_matches:
        raise ValueError("reaction product does not match route product")
    warnings: list[str] = []
    if len(product_matches) > 1:
        warnings.append("PRODUCT_SYMMETRY_AMBIGUOUS_DETERMINISTIC_REPRESENTATIVE")
    match = product_matches[0]
    reaction_map_to_target: dict[int, int] = {}
    reaction_map_to_atomic_number: dict[int, int] = {}
    for node_index, reaction_index in enumerate(match):
        target_map = product_atom_ownership.get(node_index)
        if target_map is None:
            continue
        reaction_map = product.GetAtomWithIdx(reaction_index).GetAtomMapNum()
        if reaction_map <= 0:
            raise ValueError("mapped product omits a target-derived atom")
        if reaction_map in reaction_map_to_target:
            raise ValueError("mapped product reuses an atom-map number")
        reaction_map_to_target[int(reaction_map)] = int(target_map)
        reaction_map_to_atomic_number[int(reaction_map)] = int(
            product.GetAtomWithIdx(reaction_index).GetAtomicNum()
        )
    if set(reaction_map_to_target.values()) != set(product_atom_ownership.values()):
        raise ValueError("mapped product does not cover current target ownership")
    components: list[_ComponentProjection] = []
    distributed: list[int] = []
    for component in _split_mapped_components(reactant_side):
        mapped_smiles = Chem.MolToSmiles(
            component,
            canonical=True,
            isomericSmiles=True,
        )
        serialized_component = Chem.MolFromSmiles(mapped_smiles)
        if serialized_component is None:
            raise ValueError("mapped precursor component could not be serialized")
        ownership = []
        for atom in serialized_component.GetAtoms():
            reaction_map = int(atom.GetAtomMapNum())
            target_map = reaction_map_to_target.get(reaction_map)
            if target_map is not None:
                if (
                    atom.GetAtomicNum()
                    != reaction_map_to_atomic_number[reaction_map]
                ):
                    raise ValueError(
                        "reaction atom mapping changes target element"
                    )
                ownership.append((atom.GetIdx(), target_map))
                distributed.append(target_map)
        components.append(
            _ComponentProjection(
                mapped_smiles=mapped_smiles,
                canonical_smiles=_canonical_unmapped(serialized_component),
                atom_ownership=tuple(sorted(ownership)),
            )
        )
    if len(distributed) != len(set(distributed)):
        raise ValueError("reactant components duplicate target atom ownership")
    if set(distributed) != set(product_atom_ownership.values()):
        raise ValueError("reactant components do not cover current target ownership")
    return (
        tuple(components),
        materialized.evidence,
        materialized.confidence,
        tuple(warnings),
    )


def _boundary_bonds(
    target_smiles: str,
    module_atom_sets: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, int, str], ...]:
    canonical, _, _, index_to_map = analyze_partition_target(target_smiles)
    molecule = Chem.MolFromSmiles(canonical)
    if molecule is None:
        raise ValueError("target could not be parsed")
    owner = {
        atom_map: module_index
        for module_index, block in enumerate(module_atom_sets)
        for atom_map in block
    }
    values = []
    for bond in molecule.GetBonds():
        left = index_to_map.get(bond.GetBeginAtomIdx())
        right = index_to_map.get(bond.GetEndAtomIdx())
        if left is None or right is None or owner[left] == owner[right]:
            continue
        values.append((min(left, right), max(left, right), str(bond.GetBondType())))
    return tuple(sorted(values))


def project_reaction_to_target(
    target_smiles: str,
    reaction_smiles: str,
) -> ReactionPartitionProjection:
    """Project one reverse transformation into role-neutral target modules."""

    canonical, _, target_atoms, _ = analyze_partition_target(target_smiles)
    target = Chem.MolFromSmiles(canonical)
    if target is None:
        raise ValueError("target could not be parsed")
    map_by_index = {
        reference.atom_index: reference.atom_map for reference in target_atoms
    }
    components, evidence, confidence, warnings = _propagate_reaction_ownership(
        canonical,
        reaction_smiles,
        map_by_index,
    )
    contributing = tuple(
        sorted(
            (component for component in components if component.target_atom_maps),
            key=lambda item: (item.target_atom_maps, item.mapped_smiles),
        )
    )
    module_atom_sets = tuple(item.target_atom_maps for item in contributing)
    if not module_atom_sets:
        raise ValueError("reaction has no target-contributing precursor")
    materialized = materialize_atom_mapping(reaction_smiles)
    if materialized is None:
        raise ValueError("reaction mapping could not be materialized")
    return ReactionPartitionProjection(
        target_smiles=canonical,
        mapped_reaction_smiles=materialized.reaction_smiles,
        module_atom_sets=module_atom_sets,
        target_boundary_bonds=_boundary_bonds(canonical, module_atom_sets),
        precursor_mapped_smiles=tuple(item.mapped_smiles for item in contributing),
        mapping_evidence=evidence,
        mapping_confidence=confidence,
        warnings=warnings,
    )


def _match_children(
    children: tuple[MoleculeOccurrenceNode, ...],
    components: tuple[_ComponentProjection, ...],
) -> tuple[tuple[MoleculeOccurrenceNode, _ComponentProjection, tuple[int, ...]], ...]:
    available = set(range(len(components)))
    assignments = []
    for child in sorted(children, key=lambda item: item.occurrence_id):
        child_canonical = canonical_smiles(child.smiles)
        if child_canonical is None:
            raise ValueError("route child could not be parsed")
        candidates = tuple(
            index
            for index in sorted(available)
            if components[index].canonical_smiles == child_canonical
        )
        if not candidates:
            raise ValueError("route child does not match a reaction precursor")
        component_index = candidates[0]
        available.remove(component_index)
        component = components[component_index]
        child_molecule = Chem.MolFromSmiles(child.smiles)
        component_molecule = Chem.MolFromSmiles(component.mapped_smiles)
        if child_molecule is None or component_molecule is None:
            raise ValueError("route child matching could not be parsed")
        matches = _full_isomorphisms(child_molecule, component_molecule)
        if not matches:
            raise ValueError("route child graph differs from reaction precursor")
        assignments.append((child, component, matches[0]))
    omitted_target_maps = {
        target_map
        for index in available
        for target_map in components[index].target_atom_maps
    }
    if omitted_target_maps:
        raise ValueError("route tree omits a target-contributing precursor")
    return tuple(assignments)


def _project_occurrence(
    node: MoleculeOccurrenceNode,
    atom_ownership: Mapping[int, int],
    unresolved: set[str],
) -> _ProjectedOccurrence:
    target_maps = tuple(sorted(set(atom_ownership.values())))
    if node.reaction is None:
        return _ProjectedOccurrence(
            node=node,
            target_atom_maps=target_maps,
            atom_ownership=tuple(sorted(atom_ownership.items())),
            children=(),
            expansion_resolved=True,
            warnings=(),
        )
    try:
        components, _, _, reaction_warnings = _propagate_reaction_ownership(
            node.smiles,
            node.reaction.reaction_smiles,
            atom_ownership,
        )
        assignments = _match_children(node.reaction.children, components)
        children = []
        symmetry_warning = False
        for child, component, match in assignments:
            if len(
                _full_isomorphisms(
                    Chem.MolFromSmiles(component.canonical_smiles),
                    Chem.MolFromSmiles(component.mapped_smiles),
                )
            ) > 1:
                symmetry_warning = True
            component_ownership = dict(component.atom_ownership)
            child_ownership = {
                child_index: component_ownership[component_index]
                for child_index, component_index in enumerate(match)
                if component_index in component_ownership
            }
            children.append(
                _project_occurrence(child, child_ownership, unresolved)
            )
        warnings = list(reaction_warnings)
        if symmetry_warning:
            warnings.append(
                "PRECURSOR_SYMMETRY_AMBIGUOUS_DETERMINISTIC_REPRESENTATIVE"
            )
        return _ProjectedOccurrence(
            node=node,
            target_atom_maps=target_maps,
            atom_ownership=tuple(sorted(atom_ownership.items())),
            children=tuple(children),
            expansion_resolved=True,
            warnings=tuple(sorted(set(warnings))),
        )
    except ValueError:
        unresolved.add(node.occurrence_id)
        return _ProjectedOccurrence(
            node=node,
            target_atom_maps=target_maps,
            atom_ownership=tuple(sorted(atom_ownership.items())),
            children=(),
            expansion_resolved=False,
            warnings=("ROUTE_ATOM_PROJECTION_UNRESOLVED",),
        )


def _frontier_at_depth(
    projected: _ProjectedOccurrence,
    depth: int,
) -> tuple[tuple[_ProjectedOccurrence, ...], bool]:
    if projected.node.depth >= depth or projected.node.reaction is None:
        return (projected,), True
    if not projected.expansion_resolved:
        return (), False
    values = []
    complete = True
    for child in projected.children:
        child_values, child_complete = _frontier_at_depth(child, depth)
        values.extend(child_values)
        complete = complete and child_complete
    return tuple(values), complete


def _validate_projected_elements(
    projected: _ProjectedOccurrence,
    target_elements: Mapping[int, str],
) -> None:
    """Reject ownership paths that change a target atom's element."""

    molecule = Chem.MolFromSmiles(projected.node.smiles)
    if molecule is None:
        raise ValueError("projected route molecule could not be parsed")
    for atom_index, target_map in projected.atom_ownership:
        if atom_index >= molecule.GetNumAtoms():
            raise ValueError("projected atom ownership index is out of range")
        actual = molecule.GetAtomWithIdx(atom_index).GetSymbol()
        expected = target_elements[target_map]
        if actual != expected:
            raise ValueError(
                "route projection changes target element identity: "
                f"{projected.node.occurrence_id}:atom-{atom_index}:"
                f"map-{target_map}:{actual}!={expected}"
            )
    for child in projected.children:
        _validate_projected_elements(child, target_elements)


def _latent_state(
    tree_id: str,
    projected: _ProjectedOccurrence,
    partition: SyntheticPartition,
) -> LatentModuleState:
    target_maps = projected.target_atom_maps
    module_ids = tuple(
        module.module_id
        for module in partition.modules
        if set(module.target_atom_maps).issubset(target_maps)
    )
    molecule = Chem.MolFromSmiles(projected.node.smiles)
    owned_indices = {index for index, _ in projected.atom_ownership}
    non_target = []
    if molecule is not None:
        non_target = [
            NonTargetAtom(
                element=atom.GetSymbol(),
                atom_map=(atom.GetAtomMapNum() or None),
                classification="UNKNOWN",
            )
            for atom in molecule.GetAtoms()
            if atom.GetIdx() not in owned_indices
        ]
    return LatentModuleState(
        latent_state_id=digest(
            "SPLATENT1",
            tree_id,
            projected.node.occurrence_id,
            ",".join(map(str, target_maps)),
        ),
        module_ids=module_ids,
        mapped_smiles=projected.node.smiles,
        target_atom_maps=target_maps,
        non_target_atoms=tuple(non_target),
        evidence_status="route_projected",
        source_occurrence_id=projected.node.occurrence_id,
    )


def project_route_partitions(tree: ReactionRouteTree) -> RoutePartitionProjection:
    """Project every complete level frontier of an observed or planned route."""

    canonical, _, target_atoms, _ = analyze_partition_target(tree.target_smiles)
    canonical_root = Chem.MolFromSmiles(canonical)
    route_root = Chem.MolFromSmiles(tree.root.smiles)
    if canonical_root is None or route_root is None:
        raise ValueError("route target could not be parsed")
    root_matches = _full_isomorphisms(canonical_root, route_root)
    if not root_matches:
        raise ValueError("route root does not match its declared target")
    canonical_map_by_index = {
        reference.atom_index: reference.atom_map for reference in target_atoms
    }
    root_ownership = {
        route_index: canonical_map_by_index[canonical_index]
        for canonical_index, route_index in enumerate(root_matches[0])
        if canonical_index in canonical_map_by_index
    }
    unresolved: set[str] = set()
    projected_root = _project_occurrence(tree.root, root_ownership, unresolved)
    _validate_projected_elements(
        projected_root,
        {reference.atom_map: reference.element for reference in target_atoms},
    )
    evidence_level = "E4" if tree.route_kind == "observed" else "E3"
    source_kind = (
        "observed_route_frontier"
        if tree.route_kind == "observed"
        else "planned_route_frontier"
    )
    frontiers = []
    projection_warnings: set[str] = set()
    for depth in range(tree.maximum_depth + 1):
        components, complete = _frontier_at_depth(projected_root, depth)
        if not complete:
            projection_warnings.add("ROUTE_FRONTIER_PROJECTION_PARTIAL")
            continue
        target_components = tuple(
            component for component in components if component.target_atom_maps
        )
        module_sets = tuple(
            sorted(
                (component.target_atom_maps for component in target_components),
                key=lambda item: (item[0], len(item), item),
            )
        )
        if not module_sets:
            continue
        boundary_bonds = _boundary_bonds(canonical, module_sets)
        hypothesis = InterfaceHypothesis(
            interface_kind=(
                "route_frontier_boundary" if boundary_bonds else "unary_state"
            ),
            target_bonds=boundary_bonds,
            evidence_level=evidence_level,
        )
        component_warnings = tuple(
            sorted(
                {
                    warning
                    for component in target_components
                    for warning in component.warnings
                }
            )
        )
        partition = create_synthetic_partition(
            canonical,
            module_sets,
            source_kind=source_kind,
            evidence_level=evidence_level,
            interface_hypotheses=(hypothesis,),
            generation_evidence=(tree.tree_id, f"frontier_depth:{depth}"),
            warnings=component_warnings,
        )
        occurrence_ids = tuple(
            sorted(component.node.occurrence_id for component in target_components)
        )
        latent_states = tuple(
            sorted(
                (
                    _latent_state(tree.tree_id, component, partition)
                    for component in target_components
                ),
                key=lambda state: state.latent_state_id,
            )
        )
        frontiers.append(
            RouteFrontierPartition(
                frontier_id=digest(
                    "SPFRONT1",
                    tree.tree_id,
                    str(depth),
                    *occurrence_ids,
                ),
                depth=depth,
                component_occurrence_ids=occurrence_ids,
                latent_states=latent_states,
                partition=partition,
                warnings=component_warnings,
            )
        )
    if unresolved:
        projection_warnings.add("ROUTE_ATOM_PROJECTION_HAS_UNRESOLVED_STEPS")
    return RoutePartitionProjection(
        source_tree_id=tree.tree_id,
        route_kind=tree.route_kind,
        target_smiles=canonical,
        frontiers=tuple(frontiers),
        unresolved_occurrence_ids=tuple(sorted(unresolved)),
        warnings=tuple(sorted(projection_warnings)),
    )


__all__ = [
    "PARTITION_PROJECTION_ALGORITHM_VERSION",
    "PARTITION_PROJECTION_SCHEMA_VERSION",
    "ReactionPartitionProjection",
    "RouteFrontierPartition",
    "RoutePartitionProjection",
    "project_reaction_to_target",
    "project_route_partitions",
]
