"""General component and ring topology for reaction transformations."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_edits import EditNormalizationResult
from .reaction_models import (
    ReactionCandidate,
    ReactionComponent,
    ReactionEdit,
    ReactionSiteReference,
    ReactionTopology,
)


def _cycle_rank(components: Tuple[ReactionComponent, ...]) -> Optional[int]:
    """Return total graph cycle rank, or None when any component is unparseable."""
    from rdkit import Chem

    total = 0
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            return None
        fragment_count = len(Chem.GetMolFrags(molecule, asMols=False))
        total += molecule.GetNumBonds() - molecule.GetNumAtoms() + fragment_count
    return total


def _same_component_role_groups(
    selected: Optional[ReactionCandidate],
) -> Tuple[Tuple[str, ...], ...]:
    if selected is None:
        return ()
    by_component: dict[int, list[str]] = defaultdict(list)
    for role, site in selected.role_assignments.items():
        by_component[int(site.component_index)].append(str(role))
    return tuple(
        sorted(
            tuple(sorted(roles))
            for roles in by_component.values()
            if len(roles) > 1
        )
    )


def _formed_edits(edits: Iterable[ReactionEdit]) -> Tuple[ReactionEdit, ...]:
    return tuple(
        edit
        for edit in edits
        if edit.edit_type == "formed" and edit.atom_2 is not None
    )


def build_reaction_topology(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    selected: Optional[ReactionCandidate],
    edit_result: EditNormalizationResult,
) -> Optional[ReactionTopology]:
    """Derive topology from atom-provenanced edits, independent of family names."""
    if not edit_result.edits:
        return None

    components_by_index = {
        component.component_index: component for component in reactants
    }
    participating = {
        atom.component_index
        for edit in edit_result.edits
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None and atom.side == "reactant"
    }
    formed_scopes = []
    tether_distances = []
    formed_ring_sizes = []
    from rdkit import Chem

    for edit in _formed_edits(edit_result.edits):
        assert edit.atom_2 is not None
        if edit.atom_1.component_index != edit.atom_2.component_index:
            formed_scopes.append("intermolecular")
            continue
        formed_scopes.append("intramolecular")
        component = components_by_index.get(edit.atom_1.component_index)
        molecule = parse_smiles(component.input_smiles) if component else None
        if molecule is None:
            continue
        left = int(edit.atom_1.atom_index)
        right = int(edit.atom_2.atom_index)
        if left == right or molecule.GetBondBetweenAtoms(left, right) is not None:
            continue
        try:
            path = Chem.GetShortestPath(molecule, left, right)
        except Exception:
            continue
        if path:
            distance = len(path) - 1
            tether_distances.append(distance)
            formed_ring_sizes.append(distance + 1)

    scope_set = set(formed_scopes)
    if scope_set == {"intramolecular"}:
        reaction_scope = "intramolecular"
    elif scope_set == {"intermolecular"}:
        reaction_scope = "intermolecular"
    elif len(scope_set) > 1:
        reaction_scope = "mixed"
    elif len(participating) == 1:
        reaction_scope = "unimolecular"
    else:
        reaction_scope = "unresolved"

    reactant_rank = _cycle_rank(reactants)
    product_rank = _cycle_rank(products)
    ring_count_delta = (
        product_rank - reactant_rank
        if reactant_rank is not None and product_rank is not None
        else None
    )
    role_components = (
        {
            str(role): int(site.component_index)
            for role, site in selected.role_assignments.items()
        }
        if selected is not None
        else {}
    )
    return ReactionTopology(
        reaction_scope=reaction_scope,
        participating_component_indices=tuple(sorted(participating)),
        role_component_indices=role_components,
        same_component_role_groups=_same_component_role_groups(selected),
        formed_bond_scopes=tuple(sorted(formed_scopes)),
        reactant_tether_distances=tuple(sorted(tether_distances)),
        formed_ring_sizes=tuple(sorted(formed_ring_sizes)),
        ring_count_delta=ring_count_delta,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
    )


def topology_label_prefix(topology: Optional[ReactionTopology]) -> str:
    """Return a concise display-only prefix for intramolecular transformations."""
    if topology is None or topology.reaction_scope != "intramolecular":
        return ""
    if len(topology.formed_ring_sizes) == 1:
        return f"intramolecular ({topology.formed_ring_sizes[0]}-membered ring) "
    return "intramolecular "


def assignment_component_scope(
    assignment: dict[str, ReactionSiteReference],
) -> str:
    """Describe whether assigned grammar roles share reactant components."""
    component_indices = [site.component_index for site in assignment.values()]
    if len(component_indices) < 2:
        return "unresolved"
    unique_count = len(set(component_indices))
    if unique_count == 1:
        return "intramolecular"
    if unique_count == len(component_indices):
        return "intermolecular"
    return "mixed"


__all__ = [
    "assignment_component_scope",
    "build_reaction_topology",
    "topology_label_prefix",
]
