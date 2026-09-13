"""General edit-graph coverage with state, topology and direction controls."""

from copy import deepcopy
from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_reaction_core import (
    LEVELS,
    build_shared_reaction_core,
    compare_reaction_cores,
)


def project(reaction):
    analysis = featurize_reaction(reaction)
    return build_shared_reaction_core(
        reaction,
        asdict(analysis.reaction_signature) if analysis.reaction_signature else {},
        asdict(analysis.reaction_core) if analysis.reaction_core else {},
    )


@pytest.mark.parametrize(
    "left,right",
    [
        ("CC=O>>CCO", "CCC=O>>CCCO"),
        ("CCO>>CC=O", "CCCO>>CCC=O"),
        ("CC(C)=O>>CC(C)O", "CCC(C)=O>>CCC(C)O"),
        (
            "[CH3:1][CH:2]=[CH:3][CH3:4]>>[CH3:1][CH2:2][CH2:3][CH3:4]",
            "[CH3:5][CH2:1][CH:2]=[CH:3][CH3:4]>>[CH3:5][CH2:1][CH2:2][CH2:3][CH3:4]",
        ),
        (
            "[CH3:1][S:2][S:3][CH3:4]>>[CH3:1][SH:2].[SH:3][CH3:4]",
            "[CH3:5][CH2:1][S:2][S:3][CH3:4]>>[CH3:5][CH2:1][SH:2].[SH:3][CH3:4]",
        ),
    ],
)
def test_edit_graph_generalizes_environment_without_join_assumptions(left, right):
    query, precedent = project(left), project(right)
    assert tuple(level.level for level in query.levels) == LEVELS
    assert compare_reaction_cores(query, precedent).eligible
    assert not query.levels[-1].generalized_edit_indices


@pytest.mark.parametrize(
    "left,right",
    [
        ("CC=O>>CCO", "CCO>>CC=O"),
        ("CC=O>>CCO", "CC(C)=O>>CC(C)O"),
        ("CC=O>>CCO", "CC#N>>CCN"),
        (
            "[NH2:1][CH2:2][CH2:3][CH2:4][CH2:5][Br:6]>>[NH:1]1[CH2:2][CH2:3][CH2:4][CH2:5]1",
            "[NH2:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:7][Br:6]>>[NH:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:7]1",
        ),
        ("[CH3:1][NH3+:2]>>[CH3:1][NH2:2]", "[CH3:1][NH2:2]>>[CH3:1][NH3+:2]"),
    ],
)
def test_direction_state_and_ring_closure_topology_remain_protected(left, right):
    query, precedent = project(left), project(right)
    assert query.levels and precedent.levels
    assert not compare_reaction_cores(query, precedent).eligible


def test_multiple_events_preserve_every_edit_and_separation():
    reaction = "[O:1]=[CH:2][CH2:3][CH2:4][CH:5]=[O:6]>>[OH:1][CH2:2][CH2:3][CH2:4][CH2:5][OH:6]"
    analysis = featurize_reaction(reaction)
    signature, core = (
        asdict(analysis.reaction_signature),
        asdict(analysis.reaction_core),
    )
    assert signature["event_count"] > 1
    projected = project(reaction)
    assert tuple(level.level for level in projected.levels) == LEVELS
    assert not projected.levels[-1].generalized_edit_indices
    assert not compare_reaction_cores(projected, project("CC=O>>CCO")).eligible
    # Corrupt actual observed order; no relaxed tier may hide the conflict.
    broken = deepcopy(signature)
    broken["edits"][0]["old_order"] = "TRIPLE"
    assert not build_shared_reaction_core(reaction, broken, core).levels


def test_product_only_seed_cannot_turn_dehydration_into_cyanation():
    dehydration = project("NC(=O)c1ccccc1>>N#Cc1ccccc1")
    cyanation = project("Brc1ccccc1.C#N>>N#Cc1ccccc1")
    assert len(dehydration.levels) == 3
    assert not compare_reaction_cores(dehydration, cyanation).eligible


def test_general_graph_map_and_component_order_invariance():
    first = "[CH3:1][S:2][S:3][CH3:4]>>[CH3:1][SH:2].[SH:3][CH3:4]"
    second = "[CH3:14][S:13][S:12][CH3:11]>>[CH3:14][SH:13].[SH:12][CH3:11]"
    assert [level.key for level in project(first).levels] == [
        level.key for level in project(second).levels
    ]


def test_reduction_product_stereochemistry_is_not_erased_at_broad_level():
    left = "[CH3:1][C:2](=[O:3])[CH2:4][CH3:5]>>[CH3:1][C@H:2]([OH:3])[CH2:4][CH3:5]"
    right = left.replace("[C@H:2]", "[C@@H:2]")
    query, precedent = project(left), project(right)
    assert query.levels and precedent.levels
    assert not compare_reaction_cores(query, precedent).eligible


def test_hydrogen_direction_conflict_cannot_pass_general_graph_gate():
    reaction = "CC=O>>CCO"
    analysis = featurize_reaction(reaction)
    signature, core = (
        asdict(analysis.reaction_signature),
        asdict(analysis.reaction_core),
    )
    edit = next(
        item for item in signature["edits"] if item["edit_type"] == "hydrogen_change"
    )
    edit["old_order"], edit["new_order"] = "SINGLE", None
    projected = build_shared_reaction_core(reaction, signature, core)
    assert not projected.levels
    assert projected.unavailable_reasons == ("HYDROGEN_EDIT_CONTRADICTS_GRAPH",)
