"""Build conservative reaction-retrieval descriptors without claiming edits."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Optional, Tuple

from rdkit import Chem

from .chemistry.rdkit_utils import parse_smiles
from .reaction_models import (
    REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION,
    PartialProductTransformation,
    ReactionCandidate,
    ReactionCompletenessAssessment,
    ReactionComponent,
    ReactionFallbackDescriptor,
    ReactionSignature,
    ReactionSiteReference,
)
from .reaction_signatures import reaction_signature_definition_versions


_DEFINITION_ID = "reaction_fallback_descriptor.v1"
_DEFINITION_VERSION = "1.1"
_DEFINITIONS = Path(__file__).with_name("definitions")
_FALLBACK_FRAGMENT_RULES_PATH = _DEFINITIONS / "fallback_fragments.v1.json"


def reaction_fallback_definition_versions() -> dict[str, str]:
    """Return chemistry definitions participating in fallback descriptor identity."""
    raw = _FALLBACK_FRAGMENT_RULES_PATH.read_bytes()
    with _FALLBACK_FRAGMENT_RULES_PATH.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    rule_version = (
        f"{payload.get('schema_version', 'unknown')}@sha256:"
        f"{hashlib.sha256(raw).hexdigest()[:16]}"
    )
    return {
        **reaction_signature_definition_versions(),
        _DEFINITION_ID: _DEFINITION_VERSION,
        _FALLBACK_FRAGMENT_RULES_PATH.name: rule_version,
    }


@lru_cache(maxsize=1)
def _load_fallback_fragment_rules() -> tuple[dict[str, object], ...]:
    with _FALLBACK_FRAGMENT_RULES_PATH.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if str(payload.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported fallback fragment definition schema")
    if str(payload.get("definition_id") or "") != "fallback_fragments.v1":
        raise ValueError("unexpected fallback fragment definition ID")
    rules = tuple(dict(rule) for rule in payload.get("rules") or ())
    if not rules:
        raise ValueError("fallback fragment definition requires rules")
    required_lists = (
        "transformation_types",
        "center_elements",
        "removed_attachment_elements",
        "installed_fragment_elements",
        "old_orders",
        "new_orders",
    )
    rule_ids = tuple(str(rule.get("rule_id") or "") for rule in rules)
    if any(not rule_id for rule_id in rule_ids) or len(set(rule_ids)) != len(rule_ids):
        raise ValueError("fallback fragment rules require unique non-empty IDs")
    for rule in rules:
        if any(not tuple(rule.get(field) or ()) for field in required_lists):
            raise ValueError("fallback fragment rule has an empty chemistry constraint")
        if int(rule.get("maximum_missing_atom_count") or 0) < 1:
            raise ValueError("fallback fragment rule has an invalid atom-count limit")
        for field in ("minimum_confidence", "minimum_product_heavy_atom_coverage"):
            value = float(rule.get(field) or 0.0)
            if not 0.0 < value <= 1.0:
                raise ValueError(f"fallback fragment rule has an invalid {field}")
        if not str(rule.get("condition_source_requirement_id") or ""):
            raise ValueError(
                "fallback fragment rule needs a condition-source requirement"
            )
    return rules


def _canonical_without_maps(smiles: str) -> str:
    molecule = parse_smiles(smiles)
    if molecule is None:
        return ""
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _component_tokens(
    components: Tuple[ReactionComponent, ...],
) -> tuple[str, ...]:
    return tuple(
        sorted(
            token
            for component in components
            if (token := _canonical_without_maps(component.input_smiles))
        )
    )


def _site_tokens(
    components: Tuple[ReactionComponent, ...],
) -> tuple[str, ...]:
    tokens = []
    for component in components:
        for site in component.compound_analysis.sites:
            tokens.extend(
                (
                    f"type:{site.site_type}",
                    f"signature:{site.canonical_signature}",
                    f"availability:{site.site_type}:{site.availability}",
                )
            )
    return tuple(sorted(tokens))


def _group_tokens(
    components: Tuple[ReactionComponent, ...],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    groups = []
    tags = []
    for component in components:
        for group in component.compound_analysis.functional_groups:
            groups.append(f"group:{group.group_id}")
            tags.extend(str(tag) for tag in group.tags)
    return tuple(sorted(groups)), tuple(sorted(tags))


def _context_tokens(
    components: Tuple[ReactionComponent, ...],
) -> tuple[str, ...]:
    tokens = []
    for component in components:
        for environment in component.compound_analysis.site_environments:
            profile = environment.reactivity_profile
            tokens.extend(f"shell:{value}" for value in environment.first_shell)
            tokens.extend(
                (
                    f"center:{profile.reactive_center.element}:"
                    f"{profile.reactive_center.hybridization}:"
                    f"{profile.reactive_center.formal_charge}:"
                    f"{int(profile.reactive_center.aromatic)}",
                    f"steric:{profile.steric.accessibility_class}:"
                    f"{profile.steric.approach_burden_class}",
                    f"electronic:{profile.electronic.activation_axis}:"
                    f"{profile.electronic.activation_class}",
                    f"context:{profile.context_kind}",
                )
            )
            tokens.extend(
                f"nearby:{group.get('group_id')}:{group.get('distance')}"
                for group in environment.nearby_groups
                if group.get("group_id")
            )
    return tuple(sorted(tokens))


def _atom_for_role(
    role: str,
    assignments: Mapping[str, ReactionSiteReference],
    components: Mapping[int, ReactionComponent],
) -> str:
    participant, separator, atom_role = role.partition(".")
    site = assignments.get(participant)
    if site is None:
        return f"role:{role}"
    indices = site.atom_roles.get(atom_role if separator else role) or ()
    component = components.get(int(site.component_index))
    molecule = parse_smiles(component.input_smiles) if component else None
    if not indices or molecule is None:
        return f"site:{site.site_type}:{atom_role or role}"
    atom = molecule.GetAtomWithIdx(int(indices[0]))
    return (
        f"{atom.GetSymbol()}:{atom.GetFormalCharge()}:"
        f"{int(atom.GetIsAromatic())}:{str(atom.GetHybridization())}"
    )


def _candidate_features(
    candidates: Tuple[ReactionCandidate, ...],
    reactants: Tuple[ReactionComponent, ...],
) -> tuple[
    tuple[str, ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[str, ...],
    tuple[str, ...],
]:
    components = {component.component_index: component for component in reactants}
    grammars = set()
    transformations = set()
    handles = set()
    edits = set()
    hypotheses = set()
    for candidate in candidates:
        grammars.add(str(candidate.grammar_id))
        transformations.add(str(candidate.transformation_class))
        role_tokens = tuple(
            sorted(
                f"{role}:{site.site_type}:{site.canonical_signature}"
                for role, site in candidate.role_assignments.items()
            )
        )
        handles.update(
            f"{site.site_type}:{site.canonical_signature}"
            for site in candidate.role_assignments.values()
        )
        candidate_edits = []
        for edit in candidate.predicted_bond_changes:
            atom_1 = _atom_for_role(
                edit.atom_1_role,
                candidate.role_assignments,
                components,
            )
            atom_2 = (
                _atom_for_role(
                    edit.atom_2_role,
                    candidate.role_assignments,
                    components,
                )
                if edit.atom_2_role
                else "H"
            )
            endpoints = tuple(sorted((atom_1, atom_2)))
            token = (
                f"{edit.change_type}:{endpoints[0]}~{endpoints[1]}:"
                f"{edit.old_order or 'NONE'}>{edit.new_order or 'NONE'}"
            )
            edits.add(token)
            candidate_edits.append(token)
        hypotheses.add(
            json.dumps(
                {
                    "grammar": candidate.grammar_id,
                    "transformation": candidate.transformation_class,
                    "archetype": candidate.edit_archetype,
                    "roles": role_tokens,
                    "edits": tuple(sorted(candidate_edits)),
                },
                ensure_ascii=True,
                sort_keys=True,
                separators=(",", ":"),
            )
        )
    return (
        tuple(sorted(grammars)),
        tuple(sorted(transformations)),
        tuple(sorted(handles)),
        tuple(sorted(edits)),
        tuple(sorted(hypotheses)),
    )


def _bond_inventory(
    components: Tuple[ReactionComponent, ...],
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        for bond in molecule.GetBonds():
            endpoints = sorted(
                (
                    bond.GetBeginAtom().GetSymbol(),
                    bond.GetEndAtom().GetSymbol(),
                )
            )
            counts[
                f"{endpoints[0]}-{endpoints[1]}:{str(bond.GetBondType()).upper()}"
            ] += 1
    return counts


def _element_inventory(
    components: Tuple[ReactionComponent, ...],
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in components:
        molecule = parse_smiles(component.input_smiles)
        if molecule is None:
            continue
        counts.update(
            atom.GetSymbol() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1
        )
    return counts


def _delta_tokens(
    reactant_counts: Counter[str],
    product_counts: Counter[str],
) -> tuple[str, ...]:
    tokens = []
    for token in sorted(set(reactant_counts).union(product_counts)):
        delta = product_counts[token] - reactant_counts[token]
        direction = "gained" if delta > 0 else "lost"
        tokens.extend(f"{direction}:{token}" for _ in range(abs(delta)))
    return tuple(tokens)


def _verified_edit_tokens(
    signature: Optional[ReactionSignature],
    partial: Optional[PartialProductTransformation],
) -> tuple[str, ...]:
    if signature is not None:
        return tuple(
            sorted(
                f"{edit.edit_type}:"
                f"{'-'.join(sorted((edit.atom_1.element, edit.atom_2.element if edit.atom_2 else 'H')))}:"
                f"{edit.old_order or 'NONE'}>{edit.new_order or 'NONE'}"
                for edit in signature.edits
            )
        )
    if partial is not None:
        center = partial.reactant_center.element
        removed = partial.removed_attachment.element
        added = partial.added_attachment.element
        return tuple(
            sorted(
                (
                    f"broken:{'-'.join(sorted((center, removed)))}:"
                    f"{partial.old_order}>NONE",
                    f"formed:{'-'.join(sorted((center, added)))}:"
                    f"NONE>{partial.new_order}",
                )
            )
        )
    return ()


def _condition_source_requirement(
    partial: Optional[PartialProductTransformation],
    completeness: Optional[ReactionCompletenessAssessment],
) -> tuple[Optional[str], tuple[str, ...]]:
    """Recognize bounded product fragments that conditions may legitimately supply."""
    if partial is None or completeness is None:
        return None, ()
    missing_elements = tuple(sorted(partial.missing_product_atom_elements))
    installed_elements = tuple(
        sorted(
            element
            for element, count in partial.installed_fragment.element_counts.items()
            for _ in range(count)
        )
    )
    missing_count = len(installed_elements)
    if not missing_elements or missing_elements != installed_elements:
        return None, ()
    for rule in _load_fallback_fragment_rules():
        if partial.transformation_type not in tuple(
            str(value) for value in rule.get("transformation_types") or ()
        ):
            continue
        if partial.reactant_center.element not in tuple(
            str(value) for value in rule.get("center_elements") or ()
        ):
            continue
        if partial.removed_attachment.element not in tuple(
            str(value) for value in rule.get("removed_attachment_elements") or ()
        ):
            continue
        if not set(missing_elements) <= set(
            str(value) for value in rule.get("installed_fragment_elements") or ()
        ):
            continue
        if missing_count > int(rule.get("maximum_missing_atom_count") or 0):
            continue
        if partial.confidence < float(rule.get("minimum_confidence") or 0.0):
            continue
        coverage = completeness.product_heavy_atom_coverage
        if coverage is None or coverage < float(
            rule.get("minimum_product_heavy_atom_coverage") or 0.0
        ):
            continue
        if partial.old_order not in tuple(
            str(value) for value in rule.get("old_orders") or ()
        ):
            continue
        if partial.new_order not in tuple(
            str(value) for value in rule.get("new_orders") or ()
        ):
            continue
        return (
            str(rule.get("condition_source_requirement_id") or ""),
            missing_elements,
        )
    return None, ()


def _digest(payload: Mapping[str, object]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "RFD1:" + hashlib.sha256(encoded).hexdigest()


def build_reaction_fallback_descriptor(
    *,
    reactants: Tuple[ReactionComponent, ...],
    products: Tuple[ReactionComponent, ...],
    candidates: Tuple[ReactionCandidate, ...],
    signature: Optional[ReactionSignature],
    partial_transformation: Optional[PartialProductTransformation],
    completeness: Optional[ReactionCompletenessAssessment],
    evidence_quality: str,
    warnings: Iterable[str] = (),
) -> ReactionFallbackDescriptor:
    """Build a deterministic descriptor from typed observations and hypotheses."""
    warning_values = tuple(sorted(set(str(value) for value in warnings)))
    reactant_groups, reactant_tags = _group_tokens(reactants)
    product_groups, product_tags = _group_tokens(products)
    candidate_features = _candidate_features(candidates, reactants)
    bond_delta = _delta_tokens(
        _bond_inventory(reactants),
        _bond_inventory(products),
    )
    element_delta = _delta_tokens(
        _element_inventory(reactants),
        _element_inventory(products),
    )
    verified_edits = _verified_edit_tokens(signature, partial_transformation)
    (
        condition_source_requirement_id,
        required_condition_source_elements,
    ) = _condition_source_requirement(partial_transformation, completeness)
    condition_source_required = bool(
        condition_source_requirement_id and required_condition_source_elements
    )
    if signature is not None:
        evidence_mode = "verified_signature"
        confidence = 1.0
    elif partial_transformation is not None:
        evidence_mode = "partial_product_correspondence"
        confidence = float(partial_transformation.confidence)
    elif candidates:
        evidence_mode = "candidate_hypotheses"
        confidence = 0.4
    else:
        evidence_mode = "structure_inventory_only"
        confidence = 0.25

    ineligibility_reasons = []
    completeness_status = completeness.status if completeness is not None else ""
    reactant_components = Counter(_component_tokens(reactants))
    product_components = Counter(_component_tokens(products))
    if not reactants or not products:
        ineligibility_reasons.append("missing_reaction_side")
    if completeness_status == "incomplete" and not condition_source_required:
        ineligibility_reasons.append("incomplete_product_atom_provenance")
    if evidence_quality in {
        "invalid_atom_mapping",
        "conflicting_edit_evidence",
        "conflicting_stereochemical_evidence",
    }:
        ineligibility_reasons.append(evidence_quality)
    blocking_warnings = {
        "PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES",
        "PRODUCT_MAPS_MISSING_FROM_REACTANTS",
        "UNACCOUNTED_PRODUCT_HEAVY_ATOMS",
        "MISSING_REACTANT_SUSPECTED",
        "INSUFFICIENT_REACTANT_MULTIPLICITY",
    }
    if blocking_warnings.intersection(warning_values) and not condition_source_required:
        ineligibility_reasons.append("contradicted_or_incomplete_structure")
    if any(
        warning.startswith(
            (
                "DUPLICATE_ATOM_MAP:",
                "CONTRADICTORY_MAPPED_BOND:",
                "ATOM_MAP_ELEMENT_MISMATCH:",
            )
        )
        for warning in warning_values
    ):
        ineligibility_reasons.append("invalid_atom_mapping")
    if (
        signature is None
        and not candidate_features[3]
        and not bond_delta
        and not element_delta
    ):
        ineligibility_reasons.append("no_discriminating_transformation_features")
    if (
        signature is None
        and product_components
        and all(
            count <= reactant_components.get(token, 0)
            for token, count in product_components.items()
        )
    ):
        ineligibility_reasons.append("no_product_structural_change")

    versions = reaction_fallback_definition_versions()
    payload = {
        "evidence_mode": evidence_mode,
        "reactant_components": tuple(sorted(reactant_components.elements())),
        "product_components": tuple(sorted(product_components.elements())),
        "reactant_sites": _site_tokens(reactants),
        "product_sites": _site_tokens(products),
        "reactant_groups": reactant_groups,
        "product_groups": product_groups,
        "contexts": _context_tokens(reactants),
        "candidate_grammars": candidate_features[0],
        "candidate_transformations": candidate_features[1],
        "candidate_handles": candidate_features[2],
        "candidate_edits": candidate_features[3],
        "candidate_hypotheses": candidate_features[4],
        "verified_edits": verified_edits,
        "bond_inventory_delta": bond_delta,
        "element_delta": element_delta,
        "required_condition_source_elements": required_condition_source_elements,
        "condition_source_requirement_id": condition_source_requirement_id,
        "definitions": versions,
        "schema_version": REACTION_FALLBACK_DESCRIPTOR_SCHEMA_VERSION,
    }
    return ReactionFallbackDescriptor(
        descriptor_id=_digest(payload),
        evidence_mode=evidence_mode,
        confidence=confidence,
        retrieval_eligible=not ineligibility_reasons,
        ineligibility_reasons=tuple(sorted(set(ineligibility_reasons))),
        reactant_component_tokens=payload["reactant_components"],
        product_component_tokens=payload["product_components"],
        reactant_site_tokens=payload["reactant_sites"],
        product_site_tokens=payload["product_sites"],
        reactant_group_tokens=reactant_groups,
        product_group_tokens=product_groups,
        context_tokens=payload["contexts"],
        candidate_grammar_tokens=candidate_features[0],
        candidate_transformation_tokens=candidate_features[1],
        candidate_handle_tokens=candidate_features[2],
        candidate_edit_tokens=candidate_features[3],
        candidate_hypothesis_tokens=candidate_features[4],
        verified_edit_tokens=verified_edits,
        bond_inventory_delta_tokens=bond_delta,
        element_delta_tokens=element_delta,
        compatibility_tags=tuple(sorted(set(reactant_tags).union(product_tags))),
        required_condition_source_elements=required_condition_source_elements,
        condition_source_requirement_id=condition_source_requirement_id,
        definition_versions=versions,
        warnings=tuple(
            sorted(
                set(warning_values).union(
                    {
                        "CONDITION_SOURCE_REQUIRED:"
                        f"{condition_source_requirement_id}:"
                        f"{','.join(required_condition_source_elements)}"
                    }
                    if condition_source_required
                    else set()
                )
            )
        ),
    )


__all__ = [
    "build_reaction_fallback_descriptor",
    "reaction_fallback_definition_versions",
]
