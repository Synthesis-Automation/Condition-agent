"""Structural transfer cautions for a proposed step and a recorded reaction.

This is an inventory comparison, not reaction-center mapping, a mechanism
classifier, or a feasibility score. It does not change admission or ranking.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any

from .api import analyze_molecule
from .descriptors.registry import descriptor_definition_versions
from .fragment_search import indexed_product


@dataclass(frozen=True)
class SiteInventoryDifference:
    """Difference in registered graph-derived site signatures and multiplicity."""

    site_type: str
    signature: str
    proposed_count: int
    recorded_count: int


@dataclass(frozen=True)
class PrecedentTransferReview:
    """Exact structure scope and unresolved transfer questions, without a score."""

    structure_scope: str
    same_product: bool | None
    same_precursors: bool | None
    proposed_only_components: tuple[str, ...] = ()
    recorded_only_components: tuple[str, ...] = ()
    precursor_site_differences: tuple[SiteInventoryDifference, ...] = ()
    warnings: tuple[str, ...] = ()
    definition_versions: tuple[tuple[str, str], ...] = ()
    experimental_transfer: str = "not_established"
    coverage: str = "registered_site_inventory_not_reaction_center_alignment"
    limitations: tuple[str, ...] = (
        "Exact structures do not verify a source procedure, yield or transferred conditions.",
        "Site differences may concern spectators or unused handles; inspect the actual reaction center.",
        "Equal inventories do not establish positional, electronic, stereochemical or condition compatibility.",
        "No atom correspondence or named reaction family is inferred by this comparison.",
    )
    schema_version: str = "precedent_transfer_review.v1"

    def to_dict(self) -> dict[str, Any]:
        """Return serializable structural evidence and its limitations."""
        return asdict(self)


def review_precedent_transfer(
    proposed_precursors: str, proposed_product: str,
    recorded_precursors: str, recorded_product: str,
) -> PrecedentTransferReview:
    """Compare molecular identities and detected precursor handles conservatively.

    Atom maps and serialization order are ignored; stereo, isotope, charge,
    component multiplicity and chemical form are retained. Invalid or missing
    evidence remains unresolved, never an assertion of absent differences.
    """
    versions = descriptor_definition_versions()
    try:
        p_reactants, p_product, r_reactants, r_product = (
            indexed_product(value)[0] for value in (
                proposed_precursors, proposed_product, recorded_precursors, recorded_product
            )
        )
    except (TypeError, ValueError):
        return PrecedentTransferReview(
            "unresolved", None, None, warnings=("INVALID_OR_MISSING_STRUCTURE",),
            definition_versions=versions,
        )
    same_product, same_precursors = p_product == r_product, p_reactants == r_reactants
    proposed = analyze_molecule(p_reactants, include_context_features=False)
    recorded = analyze_molecule(r_reactants, include_context_features=False)
    inventories = [Counter(
        (site.site_type, site.canonical_signature)
        for site in analysis.reactive_site_hypotheses
    ) for analysis in (proposed, recorded)]
    differences = tuple(
        SiteInventoryDifference(kind, signature, inventories[0][kind, signature], inventories[1][kind, signature])
        for kind, signature in sorted(inventories[0].keys() | inventories[1].keys())
        if inventories[0][kind, signature] != inventories[1][kind, signature]
    )
    p_components, r_components = Counter(p_reactants.split(".")), Counter(r_reactants.split("."))
    warnings = []
    if not same_precursors:
        warnings.append("PRECURSOR_STRUCTURES_DIFFER")
    if any(item.site_type == "transfer_group" for item in differences):
        warnings.append("TRANSFER_GROUP_INVENTORY_DIFFERS_CHECK_DONOR_CHEMISTRY")
    if any(item.site_type == "leaving_group" for item in differences):
        warnings.append("LEAVING_GROUP_INVENTORY_DIFFERS")
    if not same_product:
        warnings.append("PRODUCT_STRUCTURES_DIFFER_INSPECT_LOCAL_ENVIRONMENT_AND_STEREO")
    return PrecedentTransferReview(
        structure_scope=("exact_structures" if same_product and same_precursors else
                         "same_product_different_precursors" if same_product else
                         "different_product_requires_review"),
        same_product=same_product, same_precursors=same_precursors,
        proposed_only_components=tuple(sorted((p_components - r_components).elements())),
        recorded_only_components=tuple(sorted((r_components - p_components).elements())),
        precursor_site_differences=differences, warnings=tuple(warnings),
        definition_versions=versions,
    )
