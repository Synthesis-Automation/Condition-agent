import copy
import json

from condition_recommender.rules import (
    load_condition_rule_set,
    validate_condition_rule_payload,
)
from condition_recommender.rules.loader import RULE_SET_PATH


def _payload() -> dict:
    return json.loads(RULE_SET_PATH.read_text(encoding="utf-8"))


def test_clean_pd_sp2_cn_rules_load_with_explicit_fallback() -> None:
    rule_set = load_condition_rule_set()

    assert rule_set.definition_id == "pd_sp2_cn_condition_rules.v1"
    assert rule_set.selection_mode == "first_nonempty_tier"
    assert rule_set.tier_order == ("specific", "fallback")
    assert len(rule_set.rules) == 5
    assert sum(rule.selection.tier == "specific" for rule in rule_set.rules) == 4
    assert sum(rule.selection.tier == "fallback" for rule in rule_set.rules) == 1
    assert all(rule.status == "draft" for rule in rule_set.rules)


def test_rules_use_taxonomy_facts_instead_of_names_or_legacy_labels() -> None:
    payload = _payload()
    serialized = json.dumps(payload, sort_keys=True)

    assert "sp2_c_n_substitution" in serialized
    assert "pronucleophile_XH" in serialized
    for forbidden in (
        "named_family",
        "reaction_type",
        "reactant_1",
        "Alkyl-NH2",
        "@sp2_electrophiles",
        "BH_CN_",
        "yield",
        "z_score",
        "datasets/rules",
    ):
        assert forbidden not in serialized


def test_every_rule_references_a_distinct_registry_template() -> None:
    rules = load_condition_rule_set().rules
    template_ids = [
        recommendation.recipe_template_id
        for rule in rules
        for recommendation in rule.recommendations
    ]

    assert len(template_ids) == len(set(template_ids)) == 5


def test_rule_validation_rejects_display_label_predicates() -> None:
    payload = _payload()
    payload["rules"][0]["match"]["partner_constraints"][1][
        "chemist_label"
    ] = "Alkyl-NH2"

    assert any(
        "unknown_key:chemist_label" in error
        for error in validate_condition_rule_payload(payload)
    )


def test_rule_validation_rejects_unknown_taxonomy_vocabulary() -> None:
    payload = _payload()
    payload["rules"][0]["scope"]["transformation_classes_any"] = [
        "named_reaction_guess"
    ]

    assert any(
        "unknown_transformation_class:named_reaction_guess" in error
        for error in validate_condition_rule_payload(payload)
    )


def test_rule_validation_reports_malformed_h_count_without_raising() -> None:
    payload = _payload()
    payload["rules"][0]["match"]["partner_constraints"][1][
        "h_count_min"
    ] = "two"

    assert any(
        "h_count_min_must_be_integer" in error
        for error in validate_condition_rule_payload(payload)
    )


def test_active_rule_cannot_reference_a_draft_recipe_template() -> None:
    payload = copy.deepcopy(_payload())
    payload["rules"][0]["status"] = "active"

    assert any(
        "active_rule_references_draft_template" in error
        for error in validate_condition_rule_payload(payload)
    )
