"""Reactant views preserve observed before graphs without product leakage."""

from dataclasses import asdict
import json

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_core_graph import canonical_graph
from reactive_taxonomy.shared_core_reactants import reactant_core_key
from reactive_taxonomy.shared_reaction_core import build_shared_reaction_core


def project(reaction: str):
    analysis = featurize_reaction(reaction)
    return build_shared_reaction_core(
        reaction, asdict(analysis.reaction_signature), asdict(analysis.reaction_core)
    )


def test_reactant_projection_is_invariant_to_mapping_and_partner_order():
    first = project("[CH3:1][Br:2].[NH3:3]>>[CH3:1][NH2:3]")
    second = project("[NH3:42].[Br:21][CH3:10]>>[NH2:42][CH3:10]")
    assert first.levels and second.levels
    assert [reactant_core_key(x.graph_payload) for x in first.levels] == [
        reactant_core_key(x.graph_payload) for x in second.levels
    ]


def test_qualified_leaving_group_and_source_generalization_is_preserved():
    iodide = project("Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1")
    mesylate = project("CS(=O)(=O)Oc1ccccc1.C#N>>N#Cc1ccccc1")
    assert len(iodide.levels) == len(mesylate.levels) == 3
    assert reactant_core_key(iodide.levels[0].graph_payload) != reactant_core_key(
        mesylate.levels[0].graph_payload
    )
    assert reactant_core_key(iodide.levels[1].graph_payload) == reactant_core_key(
        mesylate.levels[1].graph_payload
    )


def fixture_graph(*, after="O", before="C", distance=2, stereo="R"):
    labels = [
        ["center", [before, stereo], [after]],
        ["center", ["N"], [after]],
        ["before", ["C"]],
        ["after", [after]],
        ["edit", "formed", None, "SINGLE"],
    ]
    edges = [
        (0, 2, ["before", "SINGLE", "STEREONONE"]),
        (0, 1, ["before", "center_distance", distance]),
        (1, 3, ["after", "DOUBLE", "STEREONONE"]),
        (4, 0, "endpoint"),
        (4, 1, "endpoint"),
    ]
    return canonical_graph(labels, edges)


def test_product_states_context_and_events_do_not_enter_reactant_key():
    assert reactant_core_key(fixture_graph(after="O")) == reactant_core_key(
        fixture_graph(after="S")
    )


@pytest.mark.parametrize("change", [{"before": "N"}, {"distance": 3}, {"stereo": "S"}])
def test_before_state_stereo_and_joint_site_geometry_are_preserved(change):
    assert reactant_core_key(fixture_graph()) != reactant_core_key(
        fixture_graph(**change)
    )


@pytest.mark.parametrize("encoded", ["*", "[99*]", "C"])
def test_invalid_incidence_colors_are_rejected(encoded):
    with pytest.raises(ValueError, match="core graph color"):
        reactant_core_key(
            json.dumps([[json.dumps(["atom", ["before", ["C"]]])], encoded])
        )
