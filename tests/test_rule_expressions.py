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


@pytest.mark.skipif(not rdkit_available(), reason="RDKit not available")
def test_feature_analyzer_includes_taxonomy_motif_ids():
    analyzer = FeatureAnalyzer()
    features = analyzer.analyze_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")
    assert features.get("Ar-Br") is True
    assert features.get("Ar-NH2") is True

