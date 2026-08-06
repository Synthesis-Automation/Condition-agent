"""Observation-only product correspondence for incomplete reaction records."""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from typing import Iterable, Optional, Tuple

from .chemistry.rdkit_utils import (
    parse_smiles,
    prepare_fragment_serialization_copy,
)
from .reaction_correspondence import (
    AtomPair,
    infer_partial_scaffold_correspondence_candidates,
)
from .reaction_edits import reaction_atom_reference
from .reaction_models import (
    PartialProductTransformation,
    ProductAtomProvenance,
    ProductFragmentAttachment,
    ProductFragmentSourceCandidate,
    ProductOriginGap,
    ReactionCompletenessAssessment,
    ReactionComponent,
)

PRODUCT_ORIGIN_GAP_VERSION = "1.0"


@dataclass(frozen=True)
class _UnmappedFragment:
    atom_indices: Tuple[int, ...]
    boundaries: Tuple[Tuple[int, int, str], ...]


def _component(
    components: Iterable[ReactionComponent], component_index: int
) -> Optional[ReactionComponent]:
    return next(
        (
            component
            for component in components
            if component.component_index == component_index
        ),
        None,
    )


def _exclude_unchanged_components(
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
) -> tuple[Tuple[ReactionComponent, ...], Tuple[ReactionComponent, ...], int]:
    """Remove exact side-preserved components from correspondence inputs."""
    reactants_by_identity: dict[str, list[ReactionComponent]] = defaultdict(list)
    for component in reactants:
        reactants_by_identity[component.canonical_smiles].append(component)
    for values in reactants_by_identity.values():
        values.sort(key=lambda component: component.component_index)

    excluded_reactant_indices = set()
    excluded_product_indices = set()
    for product in sorted(products, key=lambda component: component.component_index):
        matches = reactants_by_identity.get(product.canonical_smiles, [])
        if not matches:
            continue
        reactant = matches.pop(0)
        excluded_reactant_indices.add(reactant.component_index)
        excluded_product_indices.add(product.component_index)
    return (
        tuple(
            component
            for component in reactants
            if component.component_index not in excluded_reactant_indices
        ),
        tuple(
            component
            for component in products
            if component.component_index not in excluded_product_indices
        ),
        len(excluded_product_indices),
    )


def _bond_order(molecule: object, left: int, right: int) -> Optional[str]:
    bond = molecule.GetBondBetweenAtoms(int(left), int(right))
    return str(bond.GetBondType()).upper() if bond is not None else None


def _mapped_bonds_are_conserved(
    mapping: Tuple[AtomPair, ...],
    reactant_molecule: object,
    product_molecule: object,
) -> bool:
    for index, (_, reactant_left, _, product_left) in enumerate(mapping):
        for _, reactant_right, _, product_right in mapping[index + 1 :]:
            if _bond_order(
                reactant_molecule, reactant_left, reactant_right
            ) != _bond_order(product_molecule, product_left, product_right):
                return False
    return True


def _unmapped_fragments(
    molecule: object,
    mapped_indices: set[int],
) -> Tuple[_UnmappedFragment, ...]:
    unmapped = {
        int(atom.GetIdx())
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and int(atom.GetIdx()) not in mapped_indices
    }
    fragments = []
    while unmapped:
        pending = [min(unmapped)]
        atoms = set()
        while pending:
            atom_index = pending.pop()
            if atom_index in atoms:
                continue
            atoms.add(atom_index)
            atom = molecule.GetAtomWithIdx(int(atom_index))
            pending.extend(
                int(neighbor.GetIdx())
                for neighbor in atom.GetNeighbors()
                if neighbor.GetAtomicNum() > 1
                and int(neighbor.GetIdx()) in unmapped
                and int(neighbor.GetIdx()) not in atoms
            )
        unmapped.difference_update(atoms)
        boundaries = []
        for atom_index in atoms:
            atom = molecule.GetAtomWithIdx(int(atom_index))
            for neighbor in atom.GetNeighbors():
                neighbor_index = int(neighbor.GetIdx())
                if neighbor.GetAtomicNum() <= 1 or neighbor_index not in mapped_indices:
                    continue
                order = _bond_order(molecule, atom_index, neighbor_index)
                if order is not None:
                    boundaries.append((neighbor_index, atom_index, order))
        fragments.append(
            _UnmappedFragment(
                atom_indices=tuple(sorted(atoms)),
                boundaries=tuple(sorted(boundaries)),
            )
        )
    return tuple(
        sorted(
            fragments,
            key=lambda fragment: (
                fragment.boundaries,
                fragment.atom_indices,
            ),
        )
    )


def _fragment_smiles(
    molecule: object,
    atom_indices: Tuple[int, ...],
) -> str:
    from rdkit import Chem

    copied = prepare_fragment_serialization_copy(molecule, atom_indices)
    return str(
        Chem.MolFragmentToSmiles(
            copied,
            atomsToUse=list(atom_indices),
            canonical=True,
            isomericSmiles=True,
        )
    )


def _rooted_fragment_smiles(
    molecule: object,
    atom_indices: Tuple[int, ...],
    root_index: int,
    bond_order: str,
) -> str:
    from rdkit import Chem

    editable = Chem.RWMol()
    index_map = {}
    for atom_index in atom_indices:
        atom = Chem.Atom(molecule.GetAtomWithIdx(int(atom_index)))
        atom.SetAtomMapNum(0)
        index_map[atom_index] = int(editable.AddAtom(atom))
    selected = set(atom_indices)
    for bond in molecule.GetBonds():
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        if left in selected and right in selected:
            editable.AddBond(
                index_map[left],
                index_map[right],
                bond.GetBondType(),
            )
    dummy = Chem.Atom(0)
    dummy.SetNoImplicit(True)
    dummy_index = int(editable.AddAtom(dummy))
    bond_types = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
    }
    editable.AddBond(
        dummy_index,
        index_map[root_index],
        bond_types.get(bond_order, Chem.BondType.SINGLE),
    )
    fragment = editable.GetMol()
    try:
        Chem.SanitizeMol(fragment)
    except Exception:
        pass
    return str(
        Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True)
    )


def _fragment_key(rooted_smiles: str) -> str:
    digest = hashlib.sha256(
        f"{PRODUCT_ORIGIN_GAP_VERSION}:{rooted_smiles}".encode("utf-8")
    ).hexdigest()[:24]
    return f"PFG1:{digest}"


def _fragment_query(
    molecule: object,
    atom_indices: Tuple[int, ...],
) -> object:
    """Build an element/bond query while allowing donor charge redistribution."""
    from rdkit import Chem
    from rdkit.Chem import rdqueries

    editable = Chem.RWMol()
    index_map = {}
    for atom_index in atom_indices:
        atom = molecule.GetAtomWithIdx(int(atom_index))
        index_map[atom_index] = int(
            editable.AddAtom(
                rdqueries.AtomNumEqualsQueryAtom(int(atom.GetAtomicNum()))
            )
        )
    selected = set(atom_indices)
    for bond in molecule.GetBonds():
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        if left in selected and right in selected:
            editable.AddBond(
                index_map[left],
                index_map[right],
                bond.GetBondType(),
            )
    return editable.GetMol()


def _source_candidates(
    *,
    product_molecule: object,
    fragment: _UnmappedFragment,
    agents: Tuple[ReactionComponent, ...],
) -> Tuple[ProductFragmentSourceCandidate, ...]:
    query = _fragment_query(product_molecule, fragment.atom_indices)
    candidates = []
    for component in agents:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None or not molecule.HasSubstructMatch(query):
            continue
        candidates.append(
            ProductFragmentSourceCandidate(
                side="agent",
                component_index=component.component_index,
                canonical_smiles=component.canonical_smiles,
                evidence="agent_fragment_subgraph_support",
                confidence=0.7,
            )
        )
    return tuple(
        sorted(
            candidates,
            key=lambda candidate: (
                candidate.side,
                candidate.component_index,
                candidate.canonical_smiles,
            ),
        )
    )


def _internal_bond_types(
    molecule: object,
    atom_indices: Tuple[int, ...],
) -> Tuple[str, ...]:
    selected = set(atom_indices)
    values = []
    for bond in molecule.GetBonds():
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        if left not in selected or right not in selected:
            continue
        elements = sorted(
            (
                str(molecule.GetAtomWithIdx(left).GetSymbol()),
                str(molecule.GetAtomWithIdx(right).GetSymbol()),
            )
        )
        values.append(
            f"{elements[0]}-{elements[1]}:{str(bond.GetBondType()).upper()}"
        )
    return tuple(sorted(values))


def _is_acyl_center(molecule: object, center_index: int) -> bool:
    atom = molecule.GetAtomWithIdx(int(center_index))
    if atom.GetSymbol() != "C":
        return False
    return any(
        neighbor.GetSymbol() == "O"
        and _bond_order(molecule, center_index, int(neighbor.GetIdx()))
        == "DOUBLE"
        for neighbor in atom.GetNeighbors()
    )


def _observation_from_mapping(
    mapping: Tuple[AtomPair, ...],
    reactants: Tuple[ReactionComponent, ...],
    agents: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
) -> Optional[PartialProductTransformation]:
    if not mapping:
        return None
    reactant_component = _component(reactants, mapping[0][0])
    product_component = _component(products, mapping[0][2])
    if reactant_component is None or product_component is None:
        return None
    reactant_molecule = parse_smiles(reactant_component.input_smiles)
    product_molecule = parse_smiles(product_component.input_smiles)
    if reactant_molecule is None or product_molecule is None:
        return None
    if not _mapped_bonds_are_conserved(
        mapping, reactant_molecule, product_molecule
    ):
        return None

    mapped_reactant = {reactant_atom for _, reactant_atom, _, _ in mapping}
    mapped_product = {product_atom for _, _, _, product_atom in mapping}
    removed_fragments = _unmapped_fragments(
        reactant_molecule, mapped_reactant
    )
    added_fragments = _unmapped_fragments(product_molecule, mapped_product)
    if (
        len(removed_fragments) != 1
        or len(added_fragments) != 1
        or len(removed_fragments[0].boundaries) != 1
        or len(added_fragments[0].boundaries) != 1
    ):
        return None
    removed_fragment = removed_fragments[0]
    added_fragment = added_fragments[0]
    reactant_center, removed, old_order = removed_fragment.boundaries[0]
    product_center, added, new_order = added_fragment.boundaries[0]
    mapped_centers = {
        reactant_atom: product_atom
        for _, reactant_atom, _, product_atom in mapping
    }
    if mapped_centers.get(reactant_center) != product_center:
        return None

    transformation_class = (
        "acyl_heteroatom_substitution"
        if _is_acyl_center(reactant_molecule, reactant_center)
        and _is_acyl_center(product_molecule, product_center)
        else "attachment_replacement"
    )
    missing_counts = Counter(
        str(product_molecule.GetAtomWithIdx(index).GetSymbol())
        for index in added_fragment.atom_indices
    )
    missing_elements = tuple(
        element
        for element, count in sorted(missing_counts.items())
        for _ in range(int(count))
    )
    source_candidates = _source_candidates(
        product_molecule=product_molecule,
        fragment=added_fragment,
        agents=agents,
    )
    source_status = (
        "agent_supported"
        if len(source_candidates) == 1
        else "ambiguous"
        if source_candidates
        else "unresolved"
    )
    canonical_fragment_smiles = _fragment_smiles(
        product_molecule, added_fragment.atom_indices
    )
    rooted_fragment_smiles = _rooted_fragment_smiles(
        product_molecule,
        added_fragment.atom_indices,
        added,
        new_order,
    )
    product_atom_references = tuple(
        reaction_atom_reference(product_component, atom_index)
        for atom_index in added_fragment.atom_indices
    )
    installed_fragment = ProductOriginGap(
        product_component_index=product_component.component_index,
        atom_references=product_atom_references,
        internal_bond_types=_internal_bond_types(
            product_molecule, added_fragment.atom_indices
        ),
        attachments=(
            ProductFragmentAttachment(
                scaffold_atom=reaction_atom_reference(
                    product_component, product_center
                ),
                fragment_atom=reaction_atom_reference(
                    product_component, added
                ),
                bond_order=new_order,
            ),
        ),
        canonical_fragment_smiles=canonical_fragment_smiles,
        rooted_fragment_smiles=rooted_fragment_smiles,
        fragment_key=_fragment_key(rooted_fragment_smiles),
        element_counts=dict(sorted(missing_counts.items())),
        formal_charge=sum(
            int(product_molecule.GetAtomWithIdx(index).GetFormalCharge())
            for index in added_fragment.atom_indices
        ),
        source_status=source_status,
        source_candidates=source_candidates,
        evidence=(
            "agent_fragment_subgraph_support"
            if source_candidates
            else "product_only_fragment"
        ),
        confidence=0.8,
    )
    warnings = ["PARTIAL_PRODUCT_CORRESPONDENCE"]
    if source_status == "unresolved":
        warnings.append("PRODUCT_FRAGMENT_SOURCE_UNRESOLVED")
        warnings.extend(
            f"PRODUCT_ATOM_SOURCE_UNRESOLVED:{element}"
            for element in sorted(set(missing_elements))
        )
    elif source_status == "ambiguous":
        warnings.append("AMBIGUOUS_PRODUCT_FRAGMENT_SOURCES")
    else:
        warnings.append("PRODUCT_FRAGMENT_SOURCE_AGENT_SUPPORTED")
    removed_fragment_smiles = _fragment_smiles(
        reactant_molecule, removed_fragment.atom_indices
    )
    removed_rooted_smiles = _rooted_fragment_smiles(
        reactant_molecule,
        removed_fragment.atom_indices,
        removed,
        old_order,
    )
    product_atom_provenance = [
        ProductAtomProvenance(
            product_atom=reaction_atom_reference(
                product_component, product_atom
            ),
            source_kind="reactant_correspondence",
            source_atom=reaction_atom_reference(
                reactant_component, reactant_atom
            ),
            evidence="partial_scaffold_correspondence",
            confidence=0.8,
        )
        for _, reactant_atom, _, product_atom in mapping
    ]
    product_atom_provenance.extend(
        ProductAtomProvenance(
            product_atom=reference,
            source_kind=(
                "agent_supported"
                if source_status == "agent_supported"
                else "unresolved"
            ),
            source_atom=None,
            evidence=installed_fragment.evidence,
            confidence=(
                0.7 if source_status == "agent_supported" else 0.0
            ),
        )
        for reference in product_atom_references
    )
    return PartialProductTransformation(
        transformation_type=(
            "attachment_replacement"
            if len(added_fragment.atom_indices) == 1
            and len(removed_fragment.atom_indices) == 1
            else "attachment_fragment_replacement"
        ),
        transformation_class=transformation_class,
        reactant_center=reaction_atom_reference(
            reactant_component, reactant_center
        ),
        product_center=reaction_atom_reference(
            product_component, product_center
        ),
        removed_attachment=reaction_atom_reference(
            reactant_component, removed
        ),
        added_attachment=reaction_atom_reference(product_component, added),
        removed_fragment_atom_indices=removed_fragment.atom_indices,
        removed_fragment_smiles=removed_fragment_smiles,
        removed_fragment_key=_fragment_key(removed_rooted_smiles),
        installed_fragment=installed_fragment,
        product_atom_provenance=tuple(
            sorted(
                product_atom_provenance,
                key=lambda provenance: (
                    provenance.product_atom.component_index,
                    provenance.product_atom.atom_index,
                ),
            )
        ),
        old_order=old_order,
        new_order=new_order,
        conserved_atom_count=len(mapping),
        product_heavy_atom_coverage=round(
            len(mapping) / max(int(product_molecule.GetNumHeavyAtoms()), 1),
            6,
        ),
        missing_product_atom_elements=missing_elements,
        evidence="partial_product_correspondence",
        confidence=0.8,
        warnings=tuple(sorted(warnings)),
    )


def _identity(value: PartialProductTransformation) -> tuple[object, ...]:
    return (
        value.transformation_type,
        value.transformation_class,
        value.reactant_center.local_environment_id,
        value.removed_fragment_key,
        value.installed_fragment.fragment_key,
        value.old_order,
        value.new_order,
    )


def infer_partial_product_transformation(
    *,
    reactants: Tuple[ReactionComponent, ...],
    agents: Tuple[ReactionComponent, ...] = (),
    products: Tuple[ReactionComponent, ...],
    completeness: ReactionCompletenessAssessment,
) -> Optional[PartialProductTransformation]:
    """Return one consensus fragment replacement or preserve ambiguity."""
    if (
        completeness.status != "incomplete"
        or not completeness.product_element_excess
    ):
        return None
    effective_reactants, effective_products, excluded_count = (
        _exclude_unchanged_components(reactants, products)
    )
    correspondence = infer_partial_scaffold_correspondence_candidates(
        effective_reactants, effective_products
    )
    if not correspondence.valid:
        return None
    observations = tuple(
        observation
        for mapping in correspondence.candidates
        if (
            observation := _observation_from_mapping(
                mapping,
                effective_reactants,
                agents,
                effective_products,
            )
        )
        is not None
    )
    identities = {_identity(observation) for observation in observations}
    if len(identities) != 1:
        return None
    selected = min(
        observations,
        key=lambda observation: (
            observation.reactant_center.component_index,
            observation.reactant_center.atom_index,
            observation.removed_attachment.atom_index,
            observation.added_attachment.atom_index,
        ),
    )
    if not excluded_count:
        return selected
    return replace(
        selected,
        warnings=tuple(
            sorted(
                set(selected.warnings).union(
                    {
                        "UNCHANGED_PRODUCT_COMPONENTS_EXCLUDED_FROM_"
                        "PARTIAL_CORRESPONDENCE"
                    }
                )
            )
        ),
    )


__all__ = [
    "infer_partial_product_transformation",
]
