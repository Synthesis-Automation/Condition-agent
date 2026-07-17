"""Construct deterministic, type-agnostic reaction signatures."""

from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from .chemistry.rdkit_utils import parse_smiles
from .reaction_edits import EditNormalizationResult
from .reaction_models import (
    ProductConnection,
    ProductTransformation,
    ReactionAtomReference,
    ReactionCandidate,
    ReactionComponent,
    ReactionFamilyEnvironment,
    ReactionPartner,
    ReactionSignature,
    ReactionSpectatorGroup,
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
    return versions


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


def _family_partner(
    environment: Optional[ReactionFamilyEnvironment], role: str
) -> Any:
    if environment is None:
        return None
    return next((partner for partner in environment.partners if partner.role == role), None)


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
                steric=dict(overlay.steric) if overlay else (
                    dict(environment.steric) if environment else {}
                ),
                electronic=dict(overlay.electronic) if overlay else (
                    dict(environment.electronic) if environment else {}
                ),
                nearby_groups=overlay.nearby_groups if overlay else (
                    environment.nearby_groups if environment else ()
                ),
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
                atom_indices.setdefault(atom.component_index, set()).add(atom.atom_index)
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
    edit_result: EditNormalizationResult,
    family_environment: Optional[ReactionFamilyEnvironment],
    product_connection: Optional[ProductConnection],
    spectators: Tuple[ReactionSpectatorGroup, ...],
    named_family: Optional[str],
    compatible_named_families: Tuple[str, ...],
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
        sorted(_edit_token(edit, include_environment=False) for edit in edit_result.edits)
    )
    environment_edit_tokens = tuple(
        sorted(_edit_token(edit, include_environment=True) for edit in edit_result.edits)
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
    transformation_class = selected.transformation_class if selected else None
    partner_handles = tuple(sorted(
        (_partner_token(partner, environment=False) for partner in partners),
        key=_canonical_json,
    ))
    partner_environments = tuple(sorted(
        (_partner_token(partner, environment=True) for partner in partners),
        key=_canonical_json,
    ))
    spectator_tokens = tuple(
        sorted(
            (group.group_id, group.graph_distance, group.tags)
            for group in spectators
        )
    )
    exact_key = _digest(
        "L0",
        {
            "edits": environment_edit_tokens,
            "partners": partner_environments,
        },
    )
    handle_key = _digest(
        "L1", {"edits": edit_tokens, "partners": partner_handles}
    )
    transformation_key = _digest(
        "L2",
        {
            "formed": formed,
            "broken": broken,
            "order_changes": order_changes,
            "transformation_class": transformation_class,
        },
    )
    bond_key = _digest(
        "L3", {"formed": formed, "broken": broken, "order_changes": order_changes}
    )
    environment_key = _digest(
        "L4",
        {"partners": partner_environments, "spectators": spectator_tokens},
    )
    definition_versions = _definition_versions()
    signature_id = _digest(
        "RS1",
        {
            "keys": (exact_key, handle_key, transformation_key, bond_key, environment_key),
            "definitions": definition_versions,
            "schema_version": "1.0",
        },
        length=64,
    )
    product_transformation = ProductTransformation(
        edits=edit_result.edits,
        formed_connection_labels=(product_connection.concise_label,)
        if product_connection
        else (),
        concise_label=(
            product_connection.concise_label
            if product_connection
            else contextual_product_label
        ),
        exact_product_verified=bool(
            selected and selected.verification == "exact_product_reconstruction"
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
        edits=edit_result.edits,
        partners=partners,
        product_transformation=product_transformation,
        transformation_class=transformation_class,
        transformation_confidence=edit_result.confidence,
        named_family=named_family,
        family_confidence=1.0 if named_family else 0.0,
        compatible_named_families=tuple(sorted(set(compatible_named_families))),
        spectator_groups=spectators,
        global_descriptors={
            "edit_count": len(edit_result.edits),
            "partner_count": len(partners),
            "spectator_count": len(spectators),
        },
        warnings=tuple(sorted(set(warnings).union(edit_result.warnings))),
        evidence_quality=edit_result.evidence,
        definition_versions=definition_versions,
    )


__all__ = ["build_reaction_signature"]
