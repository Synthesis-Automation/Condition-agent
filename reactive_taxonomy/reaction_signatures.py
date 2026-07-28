"""Construct deterministic, type-agnostic reaction signatures."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_edits import EditNormalizationResult
from .reaction_events import build_reaction_events
from .reaction_correspondence import REACTION_CORRESPONDENCE_VERSION
from .reaction_models import (
    ProductConnection,
    ProductTransformation,
    REACTION_SIGNATURE_SCHEMA_VERSION,
    ReactionAtomReference,
    ReactionCandidate,
    ReactionCompletenessAssessment,
    ReactionComponent,
    ReactionFamilyEnvironment,
    ReactionPartner,
    ReactionSignature,
    ReactionSpectatorGroup,
    ReactionStereoChange,
    ReactionTopology,
)

_DEFINITIONS = Path(__file__).with_name("definitions")
_SIGNATURE_RULES_PATH = _DEFINITIONS / "signature_features.v1.json"


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def _digest(prefix: str, value: Any, *, length: int = 24) -> str:
    encoded = _canonical_json(value).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()[:length]


@lru_cache(maxsize=1)
def _definition_versions() -> Dict[str, str]:
    versions: Dict[str, str] = {}
    for filename in (
        "handles.v1.json",
        "reaction_grammars.v1.json",
        "descriptor_rules.v1.json",
        "signature_features.v1.json",
    ):
        path = _DEFINITIONS / filename
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        versions[filename] = str(payload.get("schema_version") or "unknown")
    versions["reaction_correspondence"] = REACTION_CORRESPONDENCE_VERSION
    return dict(sorted(versions.items()))


def reaction_signature_definition_versions() -> Dict[str, str]:
    """Return the definition versions that participate in signature identity."""
    return dict(_definition_versions())


def _unmapped_canonical(component: ReactionComponent) -> str:
    from rdkit import Chem

    mol = parse_smiles(component.input_smiles)
    if mol is None:
        return component.canonical_smiles
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _atom_token(atom: ReactionAtomReference) -> Tuple[Any, ...]:
    return (
        atom.element,
        atom.formal_charge,
        atom.aromatic,
        atom.hybridization,
        atom.local_environment_id,
    )


def _edit_token(edit: Any, *, include_environment: bool) -> Tuple[Any, ...]:
    def endpoint(atom: Optional[ReactionAtomReference]) -> Tuple[Any, ...]:
        if atom is None:
            return ("H",)
        if include_environment:
            return _atom_token(atom)
        return (atom.element, atom.formal_charge, atom.aromatic, atom.hybridization)

    endpoints = tuple(sorted((endpoint(edit.atom_1), endpoint(edit.atom_2))))
    return (
        edit.edit_type,
        endpoints,
        edit.old_order or "NONE",
        edit.new_order or "NONE",
    )


def _stereo_token(
    change: ReactionStereoChange,
    *,
    include_environment: bool,
) -> Tuple[Any, ...]:
    def endpoint(atom: Optional[ReactionAtomReference]) -> Tuple[Any, ...]:
        if atom is None:
            return ()
        if include_environment:
            return _atom_token(atom)
        return (
            atom.element,
            atom.formal_charge,
            atom.aromatic,
            atom.hybridization,
        )

    endpoints = tuple(
        sorted(
            endpoint(atom)
            for atom in (change.atom_1, change.atom_2)
            if atom is not None
        )
    )
    return (
        change.stereo_type,
        endpoints,
        change.old_descriptor or "NONE",
        change.new_descriptor or "NONE",
        change.change_type,
    )


def _bond_type(edit: Any, order: Optional[str]) -> str:
    if edit.atom_2 is None:
        return f"{edit.atom_1.element}-H:{order or 'NONE'}"
    elements = sorted((edit.atom_1.element, edit.atom_2.element))
    return f"{elements[0]}-{elements[1]}:{order or 'NONE'}"


def _site_environment(component: ReactionComponent, site_id: str) -> Any:
    return next(
        (
            environment
            for environment in component.compound_analysis.site_environments
            if environment.site_id == site_id
        ),
        None,
    )


def _component(
    components: Tuple[ReactionComponent, ...], component_index: int
) -> Optional[ReactionComponent]:
    return next(
        (
            component
            for component in components
            if component.component_index == component_index
        ),
        None,
    )


def _family_partner(environment: Optional[ReactionFamilyEnvironment], role: str) -> Any:
    if environment is None:
        return None
    return next(
        (partner for partner in environment.partners if partner.role == role), None
    )


def _selected_partners(
    components: Tuple[ReactionComponent, ...],
    selected: ReactionCandidate,
    family_environment: Optional[ReactionFamilyEnvironment],
    spectators: Tuple[ReactionSpectatorGroup, ...],
) -> Tuple[ReactionPartner, ...]:
    partners = []
    for role, site in selected.role_assignments.items():
        component = _component(components, site.component_index)
        if component is None:
            continue
        overlay = _family_partner(family_environment, role)
        environment = _site_environment(component, site.site_id)
        handles = tuple(
            sorted(
                {
                    str(value)
                    for value in (
                        site.details.get("handle_token"),
                        site.details.get("center_token"),
                    )
                    if value
                }
            )
        )
        contexts = tuple(
            sorted(
                {
                    str(value)
                    for value in (
                        [site.details.get("anchor_context")]
                        + list(site.details.get("contexts") or ())
                    )
                    if value
                }
            )
        )
        identity = {
            "component": _unmapped_canonical(component),
            "site_signature": site.canonical_signature,
            "handles": handles,
            "contexts": contexts,
        }
        role_spectators = tuple(
            sorted(
                group.group_id
                for group in spectators
                if group.component_index == site.component_index
            )
        )
        partners.append(
            ReactionPartner(
                partner_id=_digest("RP1", identity),
                component_index=site.component_index,
                role=role,
                role_confidence=1.0,
                reactive_site_ids=(site.site_id,),
                handle_tokens=handles,
                anchor_contexts=contexts,
                chemist_label=site.chemist_label,
                steric=dict(overlay.steric)
                if overlay
                else (dict(environment.steric) if environment else {}),
                electronic=dict(overlay.electronic)
                if overlay
                else (dict(environment.electronic) if environment else {}),
                nearby_groups=overlay.nearby_groups
                if overlay
                else (environment.nearby_groups if environment else ()),
                spectator_group_ids=role_spectators,
                flags=tuple(overlay.flags) if overlay else (),
            )
        )
    return tuple(sorted(partners, key=lambda partner: partner.partner_id))


def _mapped_partners(
    components: Tuple[ReactionComponent, ...],
    edit_result: EditNormalizationResult,
    spectators: Tuple[ReactionSpectatorGroup, ...],
) -> Tuple[ReactionPartner, ...]:
    atom_indices: Dict[int, set[int]] = {}
    for edit in edit_result.edits:
        for atom in (edit.atom_1, edit.atom_2):
            if atom is not None and atom.side == "reactant":
                atom_indices.setdefault(atom.component_index, set()).add(
                    atom.atom_index
                )
    partners = []
    for component_index, indices in atom_indices.items():
        component = _component(components, component_index)
        if component is None:
            continue
        sites = tuple(
            sorted(
                (
                    site
                    for site in component.compound_analysis.sites
                    if indices.intersection(site.atom_indices)
                ),
                key=lambda site: site.site_id,
            )
        )
        handles = tuple(
            sorted(
                {
                    str(token)
                    for site in sites
                    for token in (
                        site.details.get("handle_token"),
                        site.details.get("center_token"),
                    )
                    if token
                }
            )
        )
        contexts = tuple(
            sorted(
                {
                    str(context)
                    for site in sites
                    for context in (
                        [site.details.get("anchor_context")]
                        + list(site.details.get("contexts") or ())
                    )
                    if context
                }
            )
        )
        environment = _site_environment(component, sites[0].site_id) if sites else None
        identity = {
            "component": _unmapped_canonical(component),
            "atoms": sorted(
                atom.local_environment_id
                for edit in edit_result.edits
                for atom in (edit.atom_1, edit.atom_2)
                if atom is not None and atom.component_index == component_index
            ),
            "handles": handles,
            "contexts": contexts,
        }
        partners.append(
            ReactionPartner(
                partner_id=_digest("RP1", identity),
                component_index=component_index,
                role=None,
                role_confidence=0.0,
                reactive_site_ids=tuple(site.site_id for site in sites),
                handle_tokens=handles,
                anchor_contexts=contexts,
                chemist_label=" + ".join(site.chemist_label for site in sites)
                or _unmapped_canonical(component),
                steric=dict(environment.steric) if environment else {},
                electronic=dict(environment.electronic) if environment else {},
                nearby_groups=environment.nearby_groups if environment else (),
                spectator_group_ids=tuple(
                    sorted(
                        group.group_id
                        for group in spectators
                        if group.component_index == component_index
                    )
                ),
            )
        )
    return tuple(sorted(partners, key=lambda partner: partner.partner_id))


def _partner_token(partner: ReactionPartner, *, environment: bool) -> Any:
    base = {
        "handles": partner.handle_tokens,
        "contexts": partner.anchor_contexts,
    }
    if environment:
        base.update(
            {
                "steric": partner.steric,
                "electronic": partner.electronic,
                "nearby_groups": partner.nearby_groups,
                "flags": partner.flags,
            }
        )
    return base


def build_reaction_signature(
    *,
    reactants: Tuple[ReactionComponent, ...],
    selected: Optional[ReactionCandidate],
    selected_events: Tuple[ReactionCandidate, ...],
    edit_result: EditNormalizationResult,
    family_environment: Optional[ReactionFamilyEnvironment],
    product_connection: Optional[ProductConnection],
    spectators: Tuple[ReactionSpectatorGroup, ...],
    named_family: Optional[str],
    compatible_named_families: Tuple[str, ...],
    topology: ReactionTopology,
    completeness: ReactionCompletenessAssessment,
    contextual_product_label: Optional[str] = None,
    warnings: Iterable[str] = (),
) -> Optional[ReactionSignature]:
    """Build a versioned signature when verified edit evidence is available."""
    if not edit_result.edits:
        return None
    partners = (
        _selected_partners(reactants, selected, family_environment, spectators)
        if selected is not None
        else _mapped_partners(reactants, edit_result, spectators)
    )
    edit_tokens = tuple(
        sorted(
            _edit_token(edit, include_environment=False) for edit in edit_result.edits
        )
    )
    environment_edit_tokens = tuple(
        sorted(
            _edit_token(edit, include_environment=True) for edit in edit_result.edits
        )
    )
    formed = tuple(
        sorted(
            _bond_type(edit, edit.new_order)
            for edit in edit_result.edits
            if edit.edit_type == "formed"
        )
    )
    broken = tuple(
        sorted(
            _bond_type(edit, edit.old_order)
            for edit in edit_result.edits
            if edit.edit_type == "broken"
        )
    )
    order_changes = tuple(
        sorted(
            f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
            for edit in edit_result.edits
            if edit.edit_type == "order_changed"
        )
    )
    stereo_tokens = tuple(
        sorted(
            _stereo_token(change, include_environment=False)
            for change in edit_result.stereo_changes
        )
    )
    environment_stereo_tokens = tuple(
        sorted(
            _stereo_token(change, include_environment=True)
            for change in edit_result.stereo_changes
        )
    )
    hydrogen_changes = tuple(
        sorted(
            f"{_bond_type(edit, edit.old_order)}>{edit.new_order or 'NONE'}"
            for edit in edit_result.edits
            if edit.edit_type == "hydrogen_change"
        )
    )
    partner_handles = tuple(
        sorted(
            (_partner_token(partner, environment=False) for partner in partners),
            key=_canonical_json,
        )
    )
    partner_environments = tuple(
        sorted(
            (_partner_token(partner, environment=True) for partner in partners),
            key=_canonical_json,
        )
    )
    spectator_tokens = tuple(
        sorted(
            (
                (group.group_id, group.graph_distance, group.tags)
                for group in spectators
            ),
            key=_canonical_json,
        )
    )
    topology_token = {
        "reaction_scope": topology.reaction_scope,
        "same_component_role_groups": topology.same_component_role_groups,
        "formed_bond_scopes": topology.formed_bond_scopes,
        "reactant_tether_distances": topology.reactant_tether_distances,
        "formed_ring_sizes": topology.formed_ring_sizes,
        "ring_count_delta": topology.ring_count_delta,
    }
    events, event_relations = build_reaction_events(
        reactants=reactants,
        edits=edit_result.edits,
        partners=partners,
        selected=selected,
        selected_events=selected_events,
        named_family=named_family,
        compatible_named_families=compatible_named_families,
        evidence=edit_result.evidence,
        confidence=edit_result.confidence,
    )
    if selected is not None:
        transformation_class = selected.transformation_class
    else:
        event_classes = tuple(
            sorted(
                {
                    event.transformation_class
                    for event in events
                    if event.transformation_class is not None
                }
            )
        )
        if len(event_classes) == 1:
            transformation_class = event_classes[0]
        elif events:
            transformation_class = (
                "generic_graph_transformation"
                if len(events) == 1
                else "generic_multi_event_graph_transformation"
            )
        else:
            transformation_class = None
    event_archetypes = tuple(
        sorted({event.edit_archetype for event in events})
    )
    edit_archetype = (
        event_archetypes[0]
        if len(event_archetypes) == 1
        else "composite"
        if event_archetypes
        else "unresolved"
    )
    exact_event_tokens = tuple(event.event_signature_key for event in events)
    handle_event_tokens = tuple(
        tuple(
            sorted(_edit_token(edit, include_environment=False) for edit in event.edits)
        )
        for event in events
    )
    transformation_event_tokens = tuple(
        (
            event.formed_bond_types,
            event.broken_bond_types,
            event.order_changes,
            event.hydrogen_changes,
            event.edit_archetype,
            event.transformation_class,
            event.topology.reaction_scope,
        )
        for event in events
    )
    exact_key = _digest(
        "L0",
        {
            "edits": environment_edit_tokens,
            "stereo": environment_stereo_tokens,
            "events": exact_event_tokens,
            "partners": partner_environments,
            "topology": topology_token,
        },
    )
    handle_key = _digest(
        "L1",
        {
            "edits": edit_tokens,
            "stereo": stereo_tokens,
            "events": handle_event_tokens,
            "partners": partner_handles,
            "topology": topology_token,
        },
    )
    transformation_key = _digest(
        "L2",
        {
            "formed": formed,
            "broken": broken,
            "order_changes": order_changes,
            "hydrogen_changes": hydrogen_changes,
            "transformation_class": transformation_class,
            "events": transformation_event_tokens,
            "topology": topology_token,
        },
    )
    bond_key = _digest(
        "L3",
        {
            "formed": formed,
            "broken": broken,
            "order_changes": order_changes,
            "hydrogen_changes": hydrogen_changes,
        },
    )
    environment_key = _digest(
        "L4",
        {"partners": partner_environments, "spectators": spectator_tokens},
    )
    definition_versions = _definition_versions()
    signature_id = _digest(
        "RS1",
        {
            "keys": (
                exact_key,
                handle_key,
                transformation_key,
                bond_key,
                environment_key,
            ),
            "definitions": definition_versions,
            "schema_version": REACTION_SIGNATURE_SCHEMA_VERSION,
        },
        length=64,
    )
    product_transformation = ProductTransformation(
        edits=edit_result.edits,
        stereo_changes=edit_result.stereo_changes,
        formed_connection_labels=(
            (product_connection.concise_label,)
            if product_connection
            else tuple(
                label
                for event in events
                for label in event.formed_connection_labels
            )
        ),
        concise_label=(
            product_connection.concise_label
            if product_connection
            else contextual_product_label
        ),
        exact_product_verified=bool(
            (selected and selected.verification == "exact_product_reconstruction")
            or edit_result.evidence.startswith("validated_atom_mapping")
            or edit_result.evidence.startswith("validated_mapping")
            or edit_result.evidence == "exact_multi_event_reconstruction"
        ),
        evidence=edit_result.evidence,
    )
    return ReactionSignature(
        signature_id=signature_id,
        exact_signature_key=exact_key,
        handle_signature_key=handle_key,
        transformation_signature_key=transformation_key,
        bond_edit_signature_key=bond_key,
        environment_signature_key=environment_key,
        edit_signature=_canonical_json(edit_tokens),
        formed_bond_types=formed,
        broken_bond_types=broken,
        order_changes=order_changes,
        hydrogen_changes=hydrogen_changes,
        stereo_changes=edit_result.stereo_changes,
        edits=edit_result.edits,
        events=events,
        event_count=len(events),
        event_scope=(
            "single_event"
            if len(events) == 1
            else "multi_event"
            if len(events) > 1
            else "unresolved"
        ),
        event_relations=event_relations,
        partners=partners,
        product_transformation=product_transformation,
        topology=topology,
        edit_archetype=edit_archetype,
        transformation_class=transformation_class,
        transformation_confidence=edit_result.confidence,
        named_family=named_family,
        family_confidence=1.0 if named_family else 0.0,
        compatible_named_families=tuple(
            sorted(
                set(compatible_named_families).union(
                    family
                    for event in events
                    for family in event.compatible_named_families
                )
            )
        ),
        spectator_groups=spectators,
        completeness=completeness,
        global_descriptors={
            "edit_count": len(edit_result.edits),
            "event_count": len(events),
            "partner_count": len(partners),
            "spectator_count": len(spectators),
        },
        warnings=tuple(
            sorted(
                set(warnings)
                .union(edit_result.warnings)
                .union(
                    warning
                    for event in events
                    for warning in event.warnings
                )
            )
        ),
        evidence_quality=edit_result.evidence,
        definition_versions=definition_versions,
    )


__all__ = ["build_reaction_signature", "reaction_signature_definition_versions"]
