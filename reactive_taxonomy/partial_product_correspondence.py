"""Observation-only product correspondence for incomplete reaction records."""

from __future__ import annotations

from typing import Iterable, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_correspondence import (
    AtomPair,
    infer_partial_scaffold_correspondence_candidates,
)
from .reaction_edits import reaction_atom_reference
from .reaction_models import (
    PartialProductTransformation,
    ReactionCompletenessAssessment,
    ReactionComponent,
)


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


def _terminal_unmapped_neighbors(
    molecule: object,
    center_index: int,
    mapped_indices: set[int],
) -> Tuple[int, ...]:
    values = []
    center = molecule.GetAtomWithIdx(int(center_index))
    for neighbor in center.GetNeighbors():
        neighbor_index = int(neighbor.GetIdx())
        if neighbor.GetAtomicNum() <= 1 or neighbor_index in mapped_indices:
            continue
        heavy_degree = sum(
            adjacent.GetAtomicNum() > 1 for adjacent in neighbor.GetNeighbors()
        )
        if heavy_degree == 1:
            values.append(neighbor_index)
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
    products: Tuple[ReactionComponent, ...],
    completeness: ReactionCompletenessAssessment,
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
    replacements = []
    for _, reactant_center, _, product_center in mapping:
        removed = _terminal_unmapped_neighbors(
            reactant_molecule, reactant_center, mapped_reactant
        )
        added = _terminal_unmapped_neighbors(
            product_molecule, product_center, mapped_product
        )
        if removed or added:
            if len(removed) != 1 or len(added) != 1:
                return None
            replacements.append(
                (reactant_center, product_center, removed[0], added[0])
            )
    if len(replacements) != 1:
        return None
    reactant_center, product_center, removed, added = replacements[0]
    removed_atom = reactant_molecule.GetAtomWithIdx(int(removed))
    added_atom = product_molecule.GetAtomWithIdx(int(added))
    removed_element = str(removed_atom.GetSymbol())
    added_element = str(added_atom.GetSymbol())
    if completeness.reactant_element_excess.get(removed_element, 0) < 1:
        return None
    if completeness.product_element_excess.get(added_element, 0) < 1:
        return None
    old_order = _bond_order(reactant_molecule, reactant_center, removed)
    new_order = _bond_order(product_molecule, product_center, added)
    if old_order is None or new_order is None:
        return None

    transformation_class = (
        "acyl_heteroatom_substitution"
        if _is_acyl_center(reactant_molecule, reactant_center)
        and _is_acyl_center(product_molecule, product_center)
        else "attachment_replacement"
    )
    missing_elements = tuple(
        element
        for element, count in sorted(
            completeness.product_element_excess.items()
        )
        for _ in range(int(count))
    )
    warnings = (
        "PARTIAL_PRODUCT_CORRESPONDENCE",
        *(f"PRODUCT_ATOM_SOURCE_UNRESOLVED:{element}" for element in missing_elements),
    )
    return PartialProductTransformation(
        transformation_type="attachment_replacement",
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


def _identity(
    value: PartialProductTransformation,
) -> tuple[str, str, str, str, str]:
    return (
        value.transformation_type,
        value.transformation_class,
        value.reactant_center.element,
        value.removed_attachment.element,
        value.added_attachment.element,
    )


def infer_partial_product_transformation(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    completeness: ReactionCompletenessAssessment,
) -> Optional[PartialProductTransformation]:
    """Return one consensus attachment replacement or preserve ambiguity."""
    if (
        completeness.status != "incomplete"
        or not completeness.product_element_excess
        or not completeness.reactant_element_excess
    ):
        return None
    correspondence = infer_partial_scaffold_correspondence_candidates(
        reactants, products
    )
    if not correspondence.valid:
        return None
    observations = tuple(
        observation
        for mapping in correspondence.candidates
        if (
            observation := _observation_from_mapping(
                mapping, reactants, products, completeness
            )
        )
        is not None
    )
    identities = {_identity(observation) for observation in observations}
    if len(identities) != 1:
        return None
    return min(
        observations,
        key=lambda observation: (
            observation.reactant_center.component_index,
            observation.reactant_center.atom_index,
            observation.removed_attachment.atom_index,
            observation.added_attachment.atom_index,
        ),
    )


def _attachment_label(
    component: ReactionComponent, atom_index: int
) -> str:
    molecule = parse_smiles(component.input_smiles)
    if molecule is None:
        return "?"
    atom = molecule.GetAtomWithIdx(int(atom_index))
    symbol = str(atom.GetSymbol())
    if symbol in {"O", "N", "S"} and int(atom.GetTotalNumHs()) > 0:
        return symbol + ("H" if int(atom.GetTotalNumHs()) == 1 else f"H{atom.GetTotalNumHs()}")
    return symbol


def render_partial_product_transformation(
    transformation: PartialProductTransformation,
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    style: str,
) -> tuple[str, str]:
    """Render a mechanism-neutral display label from the typed observation."""
    bond = "–" if style == "unicode" else "-"
    arrow = "→" if style == "unicode" else "->"
    reactant_component = _component(
        reactants, transformation.removed_attachment.component_index
    )
    product_component = _component(
        products, transformation.added_attachment.component_index
    )
    old_attachment = (
        _attachment_label(
            reactant_component, transformation.removed_attachment.atom_index
        )
        if reactant_component
        else transformation.removed_attachment.element
    )
    new_attachment = (
        _attachment_label(
            product_component, transformation.added_attachment.atom_index
        )
        if product_component
        else transformation.added_attachment.element
    )
    if transformation.transformation_class == "acyl_heteroatom_substitution":
        before = f"R{bond}C(=O){bond}{old_attachment}"
        after = f"R{bond}C(=O){bond}{new_attachment}"
    else:
        center = transformation.reactant_center.element
        before = f"{center}{bond}{old_attachment}"
        after = f"{center}{bond}{new_attachment}"
    missing = ", ".join(transformation.missing_product_atom_elements)
    note = f"[{missing} source missing]"
    concise = f"{before} {arrow} {after} {note}"
    detailed = (
        f"Partial conserved-scaffold observation: {before} {arrow} {after}; "
        f"the reported reactants do not account for {missing} in the product."
    )
    return concise, detailed


__all__ = [
    "infer_partial_product_transformation",
    "render_partial_product_transformation",
]
