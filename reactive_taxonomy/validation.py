"""Validation for the isolated reactive-taxonomy data bundle."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from .chemistry.smarts_cache import compile_smarts
from .reaction_patterns import load_reaction_pattern_definitions
from .reaction_core.substituents import load_substituent_profile_definition


DEFINITIONS_DIR = Path(__file__).with_name("definitions")


def load_taxonomy_data() -> Dict[str, Any]:
    """Load all versioned taxonomy documents."""
    payload: Dict[str, Any] = {}
    for path in sorted(DEFINITIONS_DIR.glob("*.json")):
        with path.open("r", encoding="utf-8-sig") as handle:
            payload[path.stem] = json.load(handle)
    return payload


def validate_taxonomy() -> List[str]:
    """Return validation errors; an empty list means the bundle is valid."""
    errors: List[str] = []
    try:
        payload = load_taxonomy_data()
    except (OSError, json.JSONDecodeError) as exc:
        return [f"taxonomy_load_failed:{exc}"]
    expected = {
        "context_facets.v2",
        "chemist_notation.v1",
        "descriptor_rules.v1",
        "molecular_motifs.v1",
        "site_interfaces.v2",
        "site_patterns.v2",
        "taxonomy_manifest.v4",
        "rendering.v1",
        "reactivity_descriptor_rules.v1",
        "aromatic_systems.v1",
        "reactivity_rendering.v1",
        "signature_features.v3",
        "substituent_profiles.v1",
    }
    missing = expected - set(payload)
    if missing:
        errors.append(f"missing_taxonomy_files:{','.join(sorted(missing))}")
        return errors
    descriptor_profile_rules = payload["reactivity_descriptor_rules.v1"]
    if descriptor_profile_rules.get("profile_schema_version") != "1.0":
        errors.append("invalid_reactivity_profile_schema_version")
    if not descriptor_profile_rules.get("context_kinds"):
        errors.append("missing_reactivity_context_kinds")
    if not descriptor_profile_rules.get("descriptor_statuses"):
        errors.append("missing_reactivity_descriptor_statuses")
    descriptor_steric = descriptor_profile_rules.get("steric") or {}
    descriptor_electronic = descriptor_profile_rules.get("electronic") or {}
    if int(descriptor_steric.get("radius", 0)) < 1:
        errors.append("invalid_reactivity_steric_radius")
    if int(descriptor_electronic.get("radius", 0)) < 1:
        errors.append("invalid_reactivity_electronic_radius")
    if set(descriptor_profile_rules.get("context_kinds") or ()) != set(
        (descriptor_electronic.get("activation_axes") or {}).keys()
    ):
        errors.append("reactivity_activation_axis_context_mismatch")
    aromatic_rules = payload["aromatic_systems.v1"]
    if aromatic_rules.get("classification_method") != "ring_graph_v1":
        errors.append("invalid_aromatic_system_classification_method")
    if not aromatic_rules.get("monocyclic_families"):
        errors.append("missing_aromatic_system_families")
    rendering_rules = payload["reactivity_rendering.v1"]
    if not rendering_rules.get("field_order"):
        errors.append("missing_reactivity_rendering_order")
    contexts = payload["context_facets.v2"]
    context_records = contexts.get("contexts") or []
    tokens = [str(record.get("id") or "") for record in context_records]
    precedence = contexts.get("precedence") or []
    if len(tokens) != len(set(tokens)):
        errors.append("duplicate_context_tokens")
    if set(tokens) != set(precedence):
        errors.append("context_precedence_mismatch")
    allowed_methods = {
        "mapped_smarts",
        "aromatic_ring_system",
        "atom_property",
        "element",
        "fallback",
    }
    for record in context_records:
        context_id = str(record.get("id") or "<missing>")
        if record.get("facet") not in {
            "activation",
            "scaffold",
            "element",
            "fallback",
        }:
            errors.append(f"invalid_context_facet:{context_id}")
        if not record.get("semantic_id") or not record.get("notation_id"):
            errors.append(f"incomplete_context_identity:{context_id}")
        method = str(record.get("classification_method") or "")
        if method not in allowed_methods:
            errors.append(f"invalid_context_method:{context_id}")
        if method == "mapped_smarts":
            compiled = compile_smarts(str(record.get("smarts") or ""), validate=False)
            if compiled is None:
                errors.append(f"invalid_context_smarts:{context_id}")
                continue
            available_maps = {
                int(atom.GetAtomMapNum())
                for atom in compiled.GetAtoms()
                if atom.GetAtomMapNum()
            }
            roles = record.get("atom_roles") or {}
            if "context_anchor" not in roles:
                errors.append(f"missing_context_anchor:{context_id}")
            for role, raw_maps in roles.items():
                role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
                if {int(value) for value in role_maps} - available_maps:
                    errors.append(f"unknown_context_atom_map:{context_id}:{role}")
    families = payload["site_patterns.v2"].get("site_families") or {}
    required = {
        "leaving_group",
        "pronucleophile_XH",
        "nucleophile_anion",
        "transfer_group",
        "addition_donor",
        "eliminable_pair",
        "electrophilic_center",
        "aromatic_CH",
        "unsaturated_bond",
        "dipolar_group",
        "heteroatom_bond",
    }
    if set(families) != required:
        errors.append("site_family_mismatch")
    interface_payload = payload["site_interfaces.v2"]
    if interface_payload.get("interface_schema_version") != "2.0":
        errors.append("invalid_site_interface_schema_version")
    adapters = interface_payload.get("adapters") or []
    adapter_ids = [str(adapter.get("id") or "") for adapter in adapters]
    if not adapters or len(adapter_ids) != len(set(adapter_ids)):
        errors.append("invalid_site_interface_adapter_ids")
    allowed_interfaces = {
        "reactive_link",
        "bond_capacity",
        "connection_endpoint",
        "connection_endpoints",
    }
    for adapter in adapters:
        adapter_id = str(adapter.get("id") or "<missing>")
        if adapter.get("annotation_type") not in required:
            errors.append(f"invalid_interface_annotation_type:{adapter_id}")
        if not set(adapter.get("emits") or ()) <= allowed_interfaces:
            errors.append(f"invalid_emitted_interface:{adapter_id}")
        if not adapter.get("requires_roles"):
            errors.append(f"missing_interface_roles:{adapter_id}")
    manifest = payload["taxonomy_manifest.v4"]
    identity_files = set(manifest.get("identity_definitions") or ())
    annotation_files = set(manifest.get("annotation_definitions") or ())
    manifest_files = identity_files | annotation_files
    missing_manifest_files = {
        filename
        for filename in manifest_files
        if not (DEFINITIONS_DIR / filename).is_file()
    }
    if manifest.get("taxonomy_version") != "4.0":
        errors.append("invalid_taxonomy_manifest_version")
    if missing_manifest_files:
        errors.append(
            "missing_manifest_definitions:"
            + ",".join(sorted(missing_manifest_files))
        )
    try:
        load_reaction_pattern_definitions()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"invalid_reaction_patterns:{exc}")
    try:
        load_substituent_profile_definition()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"invalid_substituent_profiles:{exc}")
    for required_identity in {"signature_features.v3.json"}:
        if required_identity not in identity_files:
            errors.append(f"missing_identity_definition:{required_identity}")
    descriptor_rules = payload["descriptor_rules.v1"].get("site_environment") or {}
    if int(descriptor_rules.get("local_group_radius", 0)) < 1:
        errors.append("invalid_local_group_radius")
    family_rules = payload["descriptor_rules.v1"].get("reaction_families") or {}
    suzuki_rules = family_rules.get("suzuki_miyaura") or {}
    if set(suzuki_rules.get("roles") or []) != {"electrophile", "transfer_partner"}:
        errors.append("invalid_suzuki_environment_roles")
    if not suzuki_rules.get("competing_site_types"):
        errors.append("missing_suzuki_competing_site_types")
    cn_rules = family_rules.get("c_n_coupling") or {}
    if set(cn_rules.get("roles") or []) != {"electrophile", "nucleophile"}:
        errors.append("invalid_cn_environment_roles")
    for family_id in ("c_o_coupling", "c_s_coupling"):
        if set((family_rules.get(family_id) or {}).get("roles") or []) != {
            "electrophile",
            "nucleophile",
        }:
            errors.append(f"invalid_heteroatom_environment_roles:{family_id}")
    group_records = payload["molecular_motifs.v1"].get("groups") or []
    group_ids = [str(record.get("id") or "") for record in group_records]
    if not group_records:
        errors.append("missing_molecular_motifs")
    if len(group_ids) != len(set(group_ids)):
        errors.append("duplicate_molecular_motif_ids")
    for record in group_records:
        group_id = str(record.get("id") or "<missing>")
        if not record.get("label"):
            errors.append(f"missing_molecular_motif_label:{group_id}")
        raw_group_patterns = record.get("smarts") or ""
        group_patterns = (
            raw_group_patterns
            if isinstance(raw_group_patterns, list)
            else [raw_group_patterns]
        )
        if not group_patterns or any(
            compile_smarts(str(smarts), validate=False) is None
            for smarts in group_patterns
        ):
            errors.append(f"invalid_molecular_motif_smarts:{group_id}")
        unknown_suppressed = set(record.get("suppresses_on_overlap") or []) - set(
            group_ids
        )
        if unknown_suppressed:
            errors.append(f"unknown_suppressed_molecular_motif:{group_id}")
    patterns = payload["site_patterns.v2"].get("patterns") or []
    pattern_ids = [str(pattern.get("id") or "") for pattern in patterns]
    if not patterns:
        errors.append("missing_handle_patterns")
    if len(pattern_ids) != len(set(pattern_ids)):
        errors.append("duplicate_handle_pattern_ids")
    for pattern in patterns:
        pattern_id = str(pattern.get("id") or "<missing>")
        if pattern.get("site_type") not in required:
            errors.append(f"invalid_pattern_site_type:{pattern_id}")
        smarts = str(pattern.get("smarts") or "")
        compiled = compile_smarts(smarts, validate=False)
        if compiled is None:
            errors.append(f"invalid_pattern_smarts:{pattern_id}")
            continue
        available_maps = {
            int(atom.GetAtomMapNum())
            for atom in compiled.GetAtoms()
            if atom.GetAtomMapNum()
        }
        atom_roles = pattern.get("atom_roles") or {}
        required_roles = {
            "leaving_group": {"anchor", "handle"},
            "pronucleophile_XH": {"center"},
            "nucleophile_anion": {"center"},
            "transfer_group": {"anchor", "center"},
            "addition_donor": {"addend_a"},
            "eliminable_pair": {"endpoint_a", "endpoint_b", "departing_a"},
            "electrophilic_center": {"center"},
            "aromatic_CH": {"center"},
            "unsaturated_bond": {"endpoint_a", "endpoint_b"},
            "dipolar_group": {
                "attachment",
                "proximal_nitrogen",
                "central_nitrogen",
                "terminal_nitrogen",
            },
            "heteroatom_bond": {
                "attachment_a",
                "endpoint_a",
                "endpoint_b",
                "attachment_b",
            },
        }.get(str(pattern.get("site_type")), set())
        if not required_roles <= set(atom_roles):
            errors.append(f"missing_required_atom_role:{pattern_id}")
        if "Csp3" in (pattern.get("tokens") or []):
            if (
                pattern.get("site_type") != "pronucleophile_XH"
                or not pattern.get("activation_token")
                or pattern.get("activation_relationship") != "alpha_to"
                or "activation_anchor" not in atom_roles
            ):
                errors.append(f"invalid_activated_c_h_pattern:{pattern_id}")
        for role, raw_maps in atom_roles.items():
            role_maps = raw_maps if isinstance(raw_maps, list) else [raw_maps]
            unknown_maps = {int(value) for value in role_maps} - available_maps
            if unknown_maps:
                errors.append(f"unknown_atom_map:{pattern_id}:{role}")
        for rule in pattern.get("suppresses") or []:
            if rule.get("site_type") not in required:
                errors.append(f"invalid_suppression_site_type:{pattern_id}")
            if rule.get("owned_role") not in (pattern.get("atom_roles") or {}):
                errors.append(f"invalid_suppression_role:{pattern_id}")
    rendering = payload["rendering.v1"]
    notation = payload["chemist_notation.v1"]
    styles = notation.get("styles") or {}
    if notation.get("default_style") not in styles:
        errors.append("invalid_default_notation_style")
    for style_id, style in styles.items():
        if not isinstance(style, dict) or not style.get("single"):
            errors.append(f"invalid_notation_style:{style_id}")
        elif not style.get("double") or not style.get("triple"):
            errors.append(f"missing_unsaturated_bond_style:{style_id}")
        elif not style.get("negative_charge") or not style.get(
            "positive_charge"
        ):
            errors.append(f"missing_charge_rendering_style:{style_id}")
    required_handle_templates = {
        "carbonyl_formaldehyde",
        "carbonyl_aldehyde",
        "carbonyl_ketone",
        "aromatic_ch",
        "nitrile",
        "organic_azide",
        "azo_bond",
        "disulfide_bond",
        "peroxide_bond",
        "addition_donor",
        "eliminable_pair",
    }
    if not required_handle_templates <= set(
        rendering.get("named_handle_templates") or {}
    ):
        errors.append("missing_named_handle_templates")
    unsaturated_templates = rendering.get("unsaturated_bond_templates") or {}
    if set(unsaturated_templates) != {"alkene", "alkyne"}:
        errors.append("missing_unsaturated_bond_templates")
    else:
        alkene_templates = unsaturated_templates["alkene"]
        if set(alkene_templates.get("left_endpoints") or {}) != {"0", "1", "2"}:
            errors.append("invalid_alkene_left_endpoint_templates")
        if set(alkene_templates.get("right_endpoints") or {}) != {"0", "1", "2"}:
            errors.append("invalid_alkene_right_endpoint_templates")
        if set(unsaturated_templates["alkyne"]) != {
            "acetylene",
            "terminal",
            "internal",
        }:
            errors.append("invalid_alkyne_templates")
    fragment_records = notation.get("fragment_notations") or []
    fragment_ids = [str(record.get("id") or "") for record in fragment_records]
    fragment_symbols = [str(record.get("symbol") or "") for record in fragment_records]
    if not fragment_ids or len(fragment_ids) != len(set(fragment_ids)):
        errors.append("invalid_fragment_notation_ids")
    if not all(fragment_symbols) or len(fragment_symbols) != len(set(fragment_symbols)):
        errors.append("invalid_fragment_notation_symbols")
    context_notations = notation.get("context_notations") or {}
    if set(context_notations) != set(tokens):
        errors.append("context_notation_mismatch")
    if any(
        str(record.get("notation_id") or "") not in context_notations
        for record in context_records
    ):
        errors.append("unknown_context_notation")
    expected_remote_classes = {
        "aryl", "heteroaryl", "alkyl", "alkenyl", "alkynyl", "acyl",
        "ring_aliphatic", "heteroatom", "generic_R",
    }
    if set(notation.get("remote_class_notations") or {}) != expected_remote_classes:
        errors.append("remote_class_notation_mismatch")
    if "HeteroAr" in json.dumps(payload, sort_keys=True):
        errors.append("removed_heteroar_notation_present")
    rendering_rules = rendering.get("xh_rules") or []
    rendering_rule_ids = [str(rule.get("id") or "") for rule in rendering_rules]
    if len(rendering_rule_ids) != len(set(rendering_rule_ids)):
        errors.append("duplicate_rendering_rule_ids")
    if any(not rule.get("template") for rule in rendering_rules):
        errors.append("missing_rendering_template")
    signature_features = payload["signature_features.v3"]
    if signature_features.get("signature_schema_version") != "3.4":
        errors.append("invalid_signature_schema_version")
    if (
        signature_features.get("environment_feature_contract")
        != "graph_local_environment.v1"
    ):
        errors.append("invalid_signature_environment_contract")
    signature_levels = signature_features.get("levels") or {}
    if any(
        "reaction_topology" not in (signature_levels.get(level) or [])
        for level in ("L0", "L1", "L2")
    ):
        errors.append("missing_signature_reaction_topology")
    return errors


__all__ = ["load_taxonomy_data", "validate_taxonomy"]
