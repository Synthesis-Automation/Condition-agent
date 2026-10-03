"""Exact-to-broad reaction evidence, source binding and bounded counts."""

from dataclasses import replace

import pytest

from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.reaction_precedent_support import (
    assess_reaction_precedent_support, qualify_reaction_precedent, reaction_evidence_core,
)
from condition_recommender.shared_core_index import build_shared_core_index, load_shared_core_index
from tests.condition_recommender.test_shared_core_retrieval import QUERY, record


def sources(tmp_path, reactions):
    index = build_generic_index([record(i, reaction) for i, reaction in enumerate(reactions, 1)])
    path = tmp_path / "shared.sqlite"
    build_shared_core_index(index, path)
    return index, load_shared_core_index(path, index)


@pytest.mark.parametrize("reaction,level", [
    (QUERY, "whole_reaction"),
    ("Ic1ccc(CC)cc1.N#C[Cu]>>N#Cc1ccc(CC)cc1", "observed_local"),
    ("Brc1ccccc1.N#C[Cu]>>N#Cc1ccccc1", "retained_local"),
    ("Brc1ncccc1.N#C[Cu]>>N#Cc1ncccc1", "retained_typed"),
])
def test_same_condition_graph_ladder_grades_actual_observations(tmp_path, reaction, level):
    index, shared = sources(tmp_path, [reaction])
    result = assess_reaction_precedent_support(QUERY, index, shared)
    assert result.status == "supported"
    assert result.strongest_level == level
    assert result.qualified_count == result.distinct_reference_count == 1
    match = result.matches[0]
    assert match.reaction_id == "reaction-1" and match.observation_id == "obs-1"
    assert match.yield_pct == 70 and match.resolved_recipe["temperature_c"] == 25
    assert bool(match.differences) == (level != "whole_reaction")


def test_product_identity_and_family_annotations_cannot_rescue_wrong_edits(tmp_path):
    index, shared = sources(tmp_path, ["NC(=O)c1ccccc1>>N#Cc1ccccc1"])
    result = assess_reaction_precedent_support(QUERY, index, shared)
    assert not result.matches and result.status == "no_qualified_precedent"
    assert result.exclusions == (("DIFFERENT_TRANSFORMATION_SAME_PRODUCT", 1),)


def test_map_and_partner_order_invariance(tmp_path):
    index, shared = sources(tmp_path, [QUERY])
    reordered = "[Cu]C#N.c1ccc(I)cc1>>c1ccc(C#N)cc1"
    result = assess_reaction_precedent_support(reordered, index, shared)
    assert result.strongest_level == "whole_reaction"
    assert result.matches == assess_reaction_precedent_support(QUERY, index, shared).matches


def test_exact_graph_grade_preserves_specified_stereochemistry():
    first = "CCBr.N[C@H](C)CC>>CCN[C@H](C)CC"
    second = "CCBr.N[C@@H](C)CC>>CCN[C@@H](C)CC"
    left, right = reaction_evidence_core(first), reaction_evidence_core(second)
    assert left and right
    evidence = qualify_reaction_precedent(
        left, right, reaction_id="opposite", reference_id="ref", reaction_smiles=second, source="fixture",
    )
    assert evidence is None or evidence.level != "whole_reaction"


def test_candidate_cap_and_display_cap_do_not_claim_exhaustive_support(tmp_path):
    index, shared = sources(tmp_path, [QUERY] * 3)
    result = assess_reaction_precedent_support(QUERY, index, shared, candidate_limit=2, match_limit=1)
    assert result.candidate_truncated and result.matches_truncated
    assert result.candidate_count == result.qualified_count == result.distinct_reference_count == 2
    assert len(result.matches) == 1


def test_projection_binding_errors_are_not_empty_results(tmp_path):
    index, shared = sources(tmp_path, [QUERY])
    with pytest.raises(ValueError, match="ARTIFACT_MISMATCH"):
        assess_reaction_precedent_support(QUERY, index, replace(shared, source_identity="wrong"))


@pytest.mark.parametrize("query", ["invalid", "CC>>CC", "[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]"])
def test_invalid_and_unresolved_queries_do_not_acquire_precedents(tmp_path, query):
    index, shared = sources(tmp_path, [QUERY])
    result = assess_reaction_precedent_support(query, index, shared)
    assert result.status == "unresolved" and not result.matches


@pytest.mark.parametrize("arguments", [{"candidate_limit": True}, {"candidate_limit": 513}, {"match_limit": 0}])
def test_limits_are_explicit_and_bounded(tmp_path, arguments):
    index, shared = sources(tmp_path, [QUERY])
    with pytest.raises(ValueError):
        assess_reaction_precedent_support(QUERY, index, shared, **arguments)
