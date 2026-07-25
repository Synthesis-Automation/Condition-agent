"""Validation for the isolated reactive-taxonomy data bundle."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from .chemistry.smarts_cache import compile_smarts


DEFINITIONS_DIR = Path(__file__).with_name("definitions")


def load_taxonomy_data() -> Dict[str, Any]:
    """Load all v1 taxonomy documents."""
    payload: Dict[str, Any] = {}
    for path in sorted(DEFINITIONS_DIR.glob("*.v1.json")):
        with path.open("r", encoding="utf-8") as handle:
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
        "contexts.v1",
        "descriptor_rules.v1",
        "functional_groups.v1",
        "handles.v1",
        "rendering.v1",
        "reaction_grammars.v1",
        "reaction_label_patterns.v1",
        "reaction_label_rendering.v1",
        "reaction_rendering.v1",
        "signature_features.v1",
    }
    missing = expected - set(payload)
    if missing:
        errors.append(f"missing_taxonomy_files:{','.join(sorted(missing))}")
        return errors
    contexts = payload["contexts.v1"]
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
    families = payload["handles.v1"].get("site_families") or {}
    required = {
        "leaving_group",
        "pronucleophile_XH",
        "transfer_group",
        "electrophilic_center",
        "aromatic_CH",
        "unsaturated_bond",
        "dipolar_group",
        "heteroatom_bond",
    }
    if set(families) != required:
        errors.append("site_family_mismatch")
    descriptor_rules = payload["descriptor_rules.v1"].get("site_environment") or {}
    if int(descriptor_rules.get("local_group_radius", 0)) < 1:
        errors.append("invalid_local_group_radius")
    if int(descriptor_rules.get("steric_radius", 0)) < 1:
        errors.append("invalid_steric_radius")
    if not isinstance(descriptor_rules.get("electronic_tag_weights"), dict):
        errors.append("invalid_electronic_tag_weights")
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
    group_records = payload["functional_groups.v1"].get("groups") or []
    group_ids = [str(record.get("id") or "") for record in group_records]
    if not group_records:
        errors.append("missing_functional_groups")
    if len(group_ids) != len(set(group_ids)):
        errors.append("duplicate_functional_group_ids")
    for record in group_records:
        group_id = str(record.get("id") or "<missing>")
        if not record.get("label"):
            errors.append(f"missing_functional_group_label:{group_id}")
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
            errors.append(f"invalid_functional_group_smarts:{group_id}")
        unknown_suppressed = set(record.get("suppresses_on_overlap") or []) - set(
            group_ids
        )
        if unknown_suppressed:
            errors.append(f"unknown_suppressed_functional_group:{group_id}")
    patterns = payload["handles.v1"].get("patterns") or []
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
            "transfer_group": {"anchor", "center"},
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
    styles = rendering.get("styles") or {}
    if rendering.get("default_style") not in styles:
        errors.append("invalid_default_rendering_style")
    for style_id, style in styles.items():
        if not isinstance(style, dict) or not style.get("bond"):
            errors.append(f"invalid_rendering_style:{style_id}")
        elif not style.get("double") or not style.get("triple"):
            errors.append(f"missing_unsaturated_bond_style:{style_id}")
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
    rendering_contexts = set((rendering.get("context_labels") or {}).keys())
    if not rendering_contexts <= set(tokens):
        errors.append("unknown_rendering_context")
    rendering_rules = rendering.get("xh_rules") or []
    rendering_rule_ids = [str(rule.get("id") or "") for rule in rendering_rules]
    if len(rendering_rule_ids) != len(set(rendering_rule_ids)):
        errors.append("duplicate_rendering_rule_ids")
    if any(not rule.get("template") for rule in rendering_rules):
        errors.append("missing_rendering_template")
    grammar_ids: List[str] = []
    allowed_operators = {
        "join_two_anchors",
        "replace_handle_with_center",
        "replace_handle_with_alkene_endpoint",
        "replace_carbonyl_oxygen_with_amine",
    }
    known_roles: Dict[str, set[str]] = {site_type: set() for site_type in required}
    for pattern in patterns:
        known_roles.get(str(pattern.get("site_type")), set()).update(
            (pattern.get("atom_roles") or {}).keys()
        )
    for grammar in payload["reaction_grammars.v1"].get("grammars") or []:
        grammar_id = str(grammar.get("id") or "<missing>")
        grammar_ids.append(grammar_id)
        roles = grammar.get("roles") or {}
        if not roles:
            errors.append(f"missing_reaction_roles:{grammar_id}")
        for role_name, constraint in roles.items():
            if constraint.get("site_type") not in required:
                errors.append(f"invalid_reaction_site_type:{grammar_id}:{role_name}")
        operator = grammar.get("operator") or {}
        if operator.get("id") not in allowed_operators:
            errors.append(f"invalid_reaction_operator:{grammar_id}")
        if grammar.get("role_relationships") and grammar.get("distinct_components"):
            errors.append(f"conflicting_component_relationship_rules:{grammar_id}")
        for relationship in grammar.get("role_relationships") or []:
            relationship_roles = relationship.get("roles") or []
            if (
                len(relationship_roles) < 2
                or any(role not in roles for role in relationship_roles)
                or relationship.get("component_relation")
                not in {"same", "different", "same_or_different"}
            ):
                errors.append(f"invalid_role_relationship:{grammar_id}")
        for pair in grammar.get("distinct_components") or []:
            if len(pair) != 2 or any(role not in roles for role in pair):
                errors.append(f"invalid_distinct_component_rule:{grammar_id}")
    if len(grammar_ids) != len(set(grammar_ids)):
        errors.append("duplicate_reaction_grammar_ids")
    reaction_rendering = payload["reaction_rendering.v1"].get("rules") or {}
    product_precedence = (
        payload["reaction_rendering.v1"].get("product_context_precedence") or []
    )
    if (
        len(product_precedence) != len(set(product_precedence))
        or not product_precedence
    ):
        errors.append("invalid_product_context_precedence")
    if "Other" not in product_precedence:
        errors.append("missing_product_context_fallback")
    fragment_indexing = payload["reaction_rendering.v1"].get("fragment_indexing") or {}
    context_symbols = fragment_indexing.get("context_symbols") or {}
    if not isinstance(context_symbols, dict) or not context_symbols:
        errors.append("invalid_reaction_fragment_context_symbols")
    alias_template = str(fragment_indexing.get("alias_template") or "")
    if "{symbol}" not in alias_template or "{index}" not in alias_template:
        errors.append("invalid_reaction_fragment_alias_template")
    if set(reaction_rendering) != set(grammar_ids):
        errors.append("reaction_rendering_coverage_mismatch")
    signature_features = payload["signature_features.v1"]
    if signature_features.get("signature_schema_version") != "1.2":
        errors.append("invalid_signature_schema_version")
    signature_levels = signature_features.get("levels") or {}
    if any(
        "reaction_topology" not in (signature_levels.get(level) or [])
        for level in ("L0", "L1", "L2")
    ):
        errors.append("missing_signature_reaction_topology")
    allowed_product_kinds = {
        "join_contexts",
        "nitrogen_substitution",
        "heteroatom_substitution",
        "terminal_alkyne",
        "activated_carbon_substitution",
        "heck_alkene",
        "chan_lam",
        "reductive_amination",
        "amide",
        "acyl_heteroatom",
        "aryl_acylation",
        "sulfonamide",
    }
    for grammar_id, rule in reaction_rendering.items():
        if rule.get("product_kind") not in allowed_product_kinds:
            errors.append(f"invalid_reaction_product_renderer:{grammar_id}")
    label_rendering = payload["reaction_label_rendering.v1"]
    label_styles = label_rendering.get("styles") or {}
    if label_rendering.get("default_style") not in label_styles:
        errors.append("invalid_default_reaction_label_style")
    if set(label_styles) != set(styles):
        errors.append("reaction_label_style_mismatch")
    edit_types = {"formed", "broken", "order_changed", "hydrogen_change"}
    clause_order = label_rendering.get("clause_order") or []
    if set(clause_order) != edit_types or len(clause_order) != len(edit_types):
        errors.append("invalid_reaction_label_clause_order")
    required_label_templates = {
        "formed",
        "broken",
        "order_changed",
        "hydrogen_gain",
        "hydrogen_loss",
        "counted_clause",
        "mapped_atom",
        "conflict",
        "exact_detail",
        "contextual_detail",
        "event_detail",
    }
    if not required_label_templates <= set(label_rendering.get("templates") or {}):
        errors.append("missing_reaction_label_templates")
    label_patterns = payload["reaction_label_patterns.v1"].get("patterns") or []
    pattern_ids = [str(pattern.get("id") or "") for pattern in label_patterns]
    allowed_pattern_matchers = {
        "substitution",
        "hydrogenation",
        "complete_alkyne_hydrogenation",
        "partial_alkyne_hydrogenation",
        "heteroatom_bond_reduction",
        "dehydrogenation",
        "heteroatom_bond_oxidation",
        "reductive_bond_cleavage",
        "intramolecular_bond_formation",
    }
    if not pattern_ids or any(not pattern_id for pattern_id in pattern_ids):
        errors.append("missing_reaction_label_pattern_id")
    if len(pattern_ids) != len(set(pattern_ids)):
        errors.append("duplicate_reaction_label_pattern_ids")
    for pattern in label_patterns:
        pattern_id = str(pattern.get("id") or "<missing>")
        if pattern.get("matcher") not in allowed_pattern_matchers:
            errors.append(f"invalid_reaction_label_pattern_matcher:{pattern_id}")
        if set(pattern.get("templates") or {}) != set(label_styles):
            errors.append(f"reaction_label_pattern_style_mismatch:{pattern_id}")
    return errors


__all__ = ["load_taxonomy_data", "validate_taxonomy"]
