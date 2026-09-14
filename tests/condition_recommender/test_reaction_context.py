"""Planning context preserves graph gates and uncertainty."""

from dataclasses import replace

import pytest

from condition_recommender.reaction_context import (
    classify_planned_reaction,
    context_core,
    load_reaction_context_policy,
    reaction_sides,
)

QUERY = "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1"


def test_source_and_departure_alternative_uses_shared_graph():
    result = classify_planned_reaction(QUERY, "Brc1ccccc1.C#N>>N#Cc1ccccc1")
    assert result.kind == "precursor_alternative"
    assert result.inputs_changed and result.advisory_only
    assert result.core_level == "retained_local"


def test_partner_order_and_map_numbers_do_not_change_input_identity():
    result = classify_planned_reaction(
        "[CH3:1][Br:2].[NH3:3]>>[CH3:1][NH2:3]",
        "[NH3:13].[Br:12][CH3:11]>>[NH2:13][CH3:11]",
    )
    assert result.kind == "supplied_reaction" and not result.inputs_changed


def test_same_product_different_bond_construction_is_separate_route():
    result = classify_planned_reaction("CCBr.N>>CCN", "CC#N>>CCN")
    assert result.kind == "route_alternative"
    assert "DIFFERENT_TRANSFORMATION_SAME_PRODUCT" in result.reasons


def test_different_product_never_supplies_original_conditions():
    result = classify_planned_reaction("CCBr.N>>CCN", "CCBr.N>>C=C")
    assert result.kind == "different_product"


@pytest.mark.parametrize("reason", ["AMBIGUOUS_CORRESPONDENCE", "CONFLICTING_EDITS"])
def test_unavailable_core_is_unresolved_even_for_identical_molecules(reason):
    core = replace(context_core(QUERY), levels=(), unavailable_reasons=(reason,))
    result = classify_planned_reaction(QUERY, QUERY, query_core=core)
    assert result.kind == "unresolved" and reason in result.reasons


@pytest.mark.parametrize("reaction", ["CCBr", ">>CCN", "CCBr>>", "bad>>CCN"])
def test_partial_or_invalid_input_is_rejected(reaction):
    with pytest.raises(ValueError):
        reaction_sides(reaction)


def test_definition_limits_are_versioned():
    policy = load_reaction_context_policy()
    assert policy["definition_id"] == "reaction_context.v1"
    assert policy["alternatives"] == 4
