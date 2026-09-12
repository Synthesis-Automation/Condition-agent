"""Observation-aligned molecular annotations for recommendation, outside identity."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping, Tuple

from reactive_taxonomy.descriptors.models import SiteReactivityProfile
from reactive_taxonomy.descriptors.registry import descriptor_definition_versions


MOLECULAR_FEATURE_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ObservedSiteFeatures:
    """A molecular hypothesis attached to explicitly observed reactant edits."""

    hypothesis_id: str
    component_index: int
    site_type: str
    active_atom_indices: Tuple[int, ...]
    edit_indices: Tuple[int, ...]
    reactivity_profile: SiteReactivityProfile
    nearby_groups: Tuple[dict[str, Any], ...]
    role: None = None


@dataclass(frozen=True)
class ReactionMolecularFeatures:
    """Versioned scoring projection; never used to construct reaction identity."""

    partners: Tuple[ObservedSiteFeatures, ...]
    status: str
    definition_versions: Tuple[Tuple[str, str], ...]
    schema_version: str = MOLECULAR_FEATURE_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Serialize the projection independently from the reaction signature."""
        return {**asdict(self), "definition_versions": dict(self.definition_versions)}


def build_reaction_molecular_features(analysis: Any) -> ReactionMolecularFeatures:
    """Select existing molecular hypotheses only where observed edits touch them.

    Unresolved profiles remain visible with their provenance; similarity excludes
    their missing descriptors. Edit hypotheses never substitute for observed edits.
    """
    observation = analysis.observation
    edits = tuple(observation.edits) if observation is not None else ()
    active: dict[tuple[int, int], set[int]] = {}
    for number, edit in enumerate(edits):
        for atom in (edit.atom_1, edit.atom_2):
            if atom is not None and atom.side == "reactant":
                active.setdefault((atom.component_index, atom.atom_index), set()).add(
                    number
                )
    partners = []
    for component in analysis.reactants:
        molecule = component.molecule_analysis
        sites = {site.hypothesis_id: site for site in molecule.reactive_site_hypotheses}
        for environment in molecule.reactive_site_environments:
            site = sites[environment.hypothesis_id]
            locus = (
                site.atom_indices
                if site.topology == "bond"
                else (environment.center_atom_index,)
            )
            observed_atoms = tuple(
                sorted(
                    index
                    for index in locus
                    if (component.component_index, index) in active
                )
            )
            if not observed_atoms:
                continue
            edit_indices = tuple(
                sorted(
                    {
                        number
                        for atom in observed_atoms
                        for number in active[(component.component_index, atom)]
                    }
                )
            )
            partners.append(
                ObservedSiteFeatures(
                    hypothesis_id=site.hypothesis_id,
                    component_index=component.component_index,
                    site_type=site.site_type,
                    active_atom_indices=observed_atoms,
                    edit_indices=edit_indices,
                    reactivity_profile=environment.reactivity_profile,
                    nearby_groups=environment.nearby_motifs,
                )
            )
    return ReactionMolecularFeatures(
        partners=tuple(partners),
        status="observed_edits" if edits else "unavailable",
        definition_versions=descriptor_definition_versions(),
    )


def validate_molecular_features(
    features: Mapping[str, Any],
    signature: Mapping[str, Any] | None = None,
) -> None:
    """Reject stale profiles and unsupported edit provenance before indexing."""
    if not features:
        return  # Pre-projection records have no molecular similarity evidence.
    if features.get("schema_version") != MOLECULAR_FEATURE_SCHEMA_VERSION:
        raise ValueError(
            "Incompatible molecular feature schema; regenerate converted records"
        )
    if features.get("definition_versions") != dict(descriptor_definition_versions()):
        raise ValueError(
            "Incompatible molecular feature definitions; regenerate converted records"
        )
    if features.get("status") not in {"observed_edits", "unavailable"}:
        raise ValueError("Invalid molecular feature observation status")
    for partner in features.get("partners") or ():
        profile = partner.get("reactivity_profile") or {}
        if profile.get("schema_version") != "1.1" or dict(
            profile.get("definition_versions") or ()
        ) != dict(descriptor_definition_versions()):
            raise ValueError(
                "Incompatible molecular profile provenance; regenerate converted records"
            )
        if features.get("status") != "observed_edits" or not partner.get(
            "edit_indices"
        ):
            raise ValueError("Molecular feature requires observed edit provenance")
        if signature:
            edits = signature.get("edits") or ()
            for number in partner["edit_indices"]:
                if not isinstance(number, int) or not 0 <= number < len(edits):
                    raise ValueError("Molecular feature references an unobserved edit")
                atoms = (edits[number].get("atom_1"), edits[number].get("atom_2"))
                if not any(
                    atom
                    and atom.get("side") == "reactant"
                    and atom.get("component_index") == partner.get("component_index")
                    and atom.get("atom_index") in partner.get("active_atom_indices", ())
                    for atom in atoms
                ):
                    raise ValueError(
                        "Molecular feature locus conflicts with observed edit"
                    )
