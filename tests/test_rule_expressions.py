import pytest

from chemtools.util.boolean_expr import evaluate
from chemtools.rule.database import RuleDatabase
from chemtools.rule.models import RuleSpec
from chemtools.rule.analyzer import FeatureAnalyzer
from chemtools.util.rdkit_helpers import rdkit_available


def test_boolean_expr_supports_motif_ids_with_hyphens():
    values = {"Ar-Br": True, "Ar-OTf": False, "Ar-F": False}
    assert evaluate("Ar-Br OR Ar-OTf", values) is True
    assert evaluate("(Ar-Br OR Ar-OTf) AND NOT Ar-F", values) is True


def test_boolean_expr_does_not_split_token_parentheses():
    values = {"Ar-B(OH)2": True}
    assert evaluate("Ar-B(OH)2", values) is True
    assert evaluate("NOT Ar-B(OH)2", values) is False


def test_rule_database_applies_if_expr_combines_all_any():
    db = RuleDatabase(
        applies_if={
            "expr": "Ar-Br",
            "all": ["must_have"],
            "any": ["option_a", "option_b"],
        }
    )
    assert db.check_applies({"Ar-Br": True, "must_have": True, "option_b": True}) is True
    assert db.check_applies({"Ar-Br": True, "must_have": False, "option_b": True}) is False
    assert db.check_applies({"Ar-Br": False, "must_have": True, "option_b": True}) is False
    assert db.check_applies({"Ar-Br": True, "must_have": True, "option_a": False, "option_b": False}) is False


def test_rule_spec_matches_expr():
    rule = RuleSpec(
        name="test",
        conditions={},
        reactant_features={"expr": "Ar-Br AND Any-NH2"},
    )
    ok, matched = rule.matches({"Ar-Br": True, "Any-NH2": True})
    assert ok is True
    assert matched


def test_rule_spec_from_dict_preserves_v2_fields():
    rule = RuleSpec.from_dict(
        {
            "id": "BR_test",
            "name": "Example Rule",
            "description": "Hello",
            "reactant_features": {"expr": "Ar-Br"},
            "conditions": {"base": "K2CO3"},
        }
    )
    assert rule.rule_id == "BR_test"
    assert rule.description == "Hello"
    dumped = rule.to_dict()
    assert dumped["id"] == "BR_test"
    assert dumped["description"] == "Hello"


def test_rule_database_v2_default_rule_fallback_uses_conditions_only():
    db = RuleDatabase.from_dict(
        {
            "schema_version": "2.0",
            "source_type": "rule",
            "metadata": {"id": "demo", "name": "Demo DB", "version": "v2"},
            "reaction": {"family": "suzuki_miyaura"},
            "applies_if": {"expr": "Ar-Br"},
            "default_rule": {
                "id": "DEF_demo",
                "description": "Default starter",
                "conditions": {"base": "K2CO3"},
            },
            "base_rules": [],
            "modifiers": [],
        }
    )
    assert db.metadata.get("name") == "Demo DB"
    assert db.check_applies({"Ar-Br": True}) is True
    matched = db.find_matching_rule({"Ar-Br": True})
    assert matched is not None
    assert matched.conditions == {"base": "K2CO3"}
    assert matched.rule_id == "DEF_demo"


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_feature_analyzer_includes_taxonomy_motif_ids():
    analyzer = FeatureAnalyzer()
    features = analyzer.analyze_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")
    assert features.get("Ar-Br") is True
    assert features.get("Ar-NH2") is True
