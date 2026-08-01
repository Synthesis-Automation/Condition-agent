"""Conservative product-atom accounting for reaction observations."""

from __future__ import annotations

from collections import Counter
from typing import Iterable, Mapping, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_edits import EditNormalizationResult
from .reaction_models import (
    ReactionCompletenessAssessment,
    ReactionComponent,
)


def _side_statistics(
    components: Iterable[ReactionComponent],
) -> tuple[Counter[str], int, int, dict[int, str]]:
    elements: Counter[str] = Counter()
    mapped_heavy_atom_count = 0
    mapped_elements: dict[int, str] = {}
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            element = str(atom.GetSymbol())
            elements[element] += 1
            map_number = int(atom.GetAtomMapNum())
            if map_number:
                mapped_heavy_atom_count += 1
                mapped_elements[map_number] = element
    return elements, sum(elements.values()), mapped_heavy_atom_count, mapped_elements


def _positive_difference(
    left: Mapping[str, int], right: Mapping[str, int]
) -> dict[str, int]:
    return {
        element: int(count - right.get(element, 0))
        for element, count in sorted(left.items())
        if count > right.get(element, 0)
    }


def build_reaction_completeness(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
) -> ReactionCompletenessAssessment:
    """Assess whether every reported product heavy atom has reactant provenance.

    Reactant atoms absent from the reported main product are retained as an
    observation rather than treated as an error because reaction records
    commonly omit leaving-group products, salts, and other byproducts.
    """
    (
        reactant_elements,
        reactant_heavy_atom_count,
        reactant_mapped_heavy_atom_count,
        reactant_maps,
    ) = _side_statistics(reactants)
    (
        product_elements,
        product_heavy_atom_count,
        product_mapped_heavy_atom_count,
        product_maps,
    ) = _side_statistics(products)
    product_element_excess = _positive_difference(
        product_elements, reactant_elements
    )
    reactant_element_excess = _positive_difference(
        reactant_elements, product_elements
    )
    shared_maps = {
        map_number
        for map_number, element in product_maps.items()
        if reactant_maps.get(map_number) == element
    }
    product_maps_missing = set(product_maps) - set(reactant_maps)
    reactant_maps_missing = set(reactant_maps) - set(product_maps)
    warnings = []
    any_mapping = bool(reactant_maps or product_maps)
    if any_mapping and (
        reactant_mapped_heavy_atom_count < reactant_heavy_atom_count
        or product_mapped_heavy_atom_count < product_heavy_atom_count
    ):
        warnings.append("PARTIAL_ATOM_MAPPING")
    if product_maps_missing:
        warnings.append("PRODUCT_MAPS_MISSING_FROM_REACTANTS")
    if reactant_maps_missing:
        warnings.append("REACTANT_MAPS_MISSING_FROM_PRODUCTS")

    insufficient_multiplicity = False
    suspected_missing_reactant = bool(product_element_excess)
    if product_element_excess:
        warnings.append("UNACCOUNTED_PRODUCT_HEAVY_ATOMS")
        warnings.append(
            "INSUFFICIENT_REACTANT_MULTIPLICITY"
            if insufficient_multiplicity
            else "MISSING_REACTANT_SUSPECTED"
        )

    correspondence_verified = edit_result.evidence in {
        "fragmented_scaffold_correspondence",
        "unique_scaffold_correspondence",
        "global_atom_correspondence",
        "conflicting_stereochemical_evidence",
    }
    mapping_verified = edit_result.evidence.startswith("validated") and bool(
        product_maps
    ) and (
        product_mapped_heavy_atom_count == product_heavy_atom_count
        and not product_maps_missing
    )
    if product_element_excess:
        status = "incomplete"
        evidence = "product_element_excess"
    elif correspondence_verified:
        status = "verified"
        evidence = edit_result.evidence
    elif mapping_verified:
        status = "verified"
        evidence = "complete_product_atom_mapping"
    else:
        status = "unresolved"
        evidence = "insufficient_atom_provenance"
        warnings.append("REACTION_COMPLETENESS_UNRESOLVED")

    if product_element_excess and product_heavy_atom_count:
        product_heavy_atom_coverage = round(
            (
                product_heavy_atom_count
                - sum(product_element_excess.values())
            )
            / product_heavy_atom_count,
            6,
        )
    elif status == "verified":
        product_heavy_atom_coverage = 1.0
    else:
        product_heavy_atom_coverage = None
    return ReactionCompletenessAssessment(
        status=status,
        reactant_heavy_atom_count=reactant_heavy_atom_count,
        product_heavy_atom_count=product_heavy_atom_count,
        reactant_element_counts=dict(sorted(reactant_elements.items())),
        product_element_counts=dict(sorted(product_elements.items())),
        product_element_excess=product_element_excess,
        reactant_element_excess=reactant_element_excess,
        reactant_mapped_heavy_atom_count=reactant_mapped_heavy_atom_count,
        product_mapped_heavy_atom_count=product_mapped_heavy_atom_count,
        shared_mapped_heavy_atom_count=len(shared_maps),
        reactant_mapping_coverage=round(
            reactant_mapped_heavy_atom_count
            / max(reactant_heavy_atom_count, 1),
            6,
        ),
        product_mapping_coverage=round(
            product_mapped_heavy_atom_count
            / max(product_heavy_atom_count, 1),
            6,
        ),
        product_heavy_atom_coverage=product_heavy_atom_coverage,
        suspected_missing_reactant=suspected_missing_reactant,
        suspected_insufficient_reactant_multiplicity=insufficient_multiplicity,
        evidence=evidence,
        warnings=tuple(sorted(set(warnings))),
    )


__all__ = ["build_reaction_completeness"]
