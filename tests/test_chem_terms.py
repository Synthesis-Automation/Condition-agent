from chemtools.taxonomy import load_registry, reset_registry
from chemtools.taxonomy.terms import evaluate_terms
from chemtools.util.boolean_expr import evaluate


def test_registry_loads_chem_terms():
    reset_registry()
    registry = load_registry()
    assert hasattr(registry, "chem_terms")
    assert "electron_poor_arene" in registry.chem_terms


def test_evaluate_terms_over_feature_map():
    features = {
        "aromatic_present": True,
        "strong_ewg_present": True,
        "aldehyde_present": True,
        "carbonyl_present": True,
        "tert_butyl_present": False,
        "quaternary_carbon_present": True,
        "ortho_substitution_present": False,
    }
    result = evaluate_terms(
        features,
        term_ids=[
            "electron_poor_arene",
            "electron_poor_aldehyde",
            "sterically_bulky_carbonyl",
        ],
    )
    assert result["electron_poor_arene"] is True
    assert result["electron_poor_aldehyde"] is True
    assert result["sterically_bulky_carbonyl"] is True


def test_boolean_expr_supports_grouping_and_token_parentheses():
    values = {"ArB(OH)2_present": True, "n": 2}
    assert evaluate("ArB(OH)2_present AND n >= 2", values) is True
    assert evaluate("(ArB(OH)2_present OR missing_present) AND n >= 2", values) is True
    assert evaluate("NOT ArB(OH)2_present", values) is False

