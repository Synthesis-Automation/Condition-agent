"""Deterministic limits on local graph projection, including dense symmetric centers."""

from dataclasses import asdict, replace

import pytest
from rdkit import Chem

from reactive_taxonomy import featurize_reaction, shared_core_graph as graphs, shared_reaction_core
from reactive_taxonomy.shared_reaction_core import SharedReactionCore, build_shared_reaction_core


def _project_centers(smiles, reverse=False):
    molecules = {index: Chem.MolFromSmiles(value) for index, value in enumerate(smiles)}
    points = [(component, atom.GetIdx()) for component, mol in molecules.items() for atom in mol.GetAtoms()]
    transitions = {}
    for component, atom in reversed(points) if reverse else points:
        state = {"component_index": component, "atom_index": atom}
        transitions[component, atom] = {"before_state": state, "after_state": state}
    return graphs.observed_edit_graph(
        molecules, molecules, transitions, (), radius=0, limit=128, ports={},
        removed=frozenset(), generalized=frozenset(),
    )


@pytest.mark.parametrize(("pairs", "vertices", "reason"), [
    (6, 13, None),
    (5, 13, "LOCAL_CORE_CENTER_DISTANCE_PAIR_LIMIT"),
    (6, 12, "LOCAL_CORE_INCIDENCE_VERTEX_LIMIT"),
])
def test_dense_graph_limit_boundaries_preserve_sparse_graphs_and_order(monkeypatch, pairs, vertices, reason):
    # Three carbon centers in one chain have six before/after distance edges.
    # Four molecular edges and three center nodes bring the incidence count to 13.
    policy = replace(graphs.shared_core_availability(), maximum_center_distance_pairs=pairs,
                     maximum_incidence_vertices=vertices)
    monkeypatch.setattr(graphs, "shared_core_availability", lambda: policy)
    assert _project_centers(["C", "C", "C"]) == _project_centers(["C", "C", "C"], reverse=True)
    if not reason:
        assert _project_centers(["CCC"]) == _project_centers(["CCC"], reverse=True)
    else:
        def unexpected(*args):
            pytest.fail("Over-budget graph reached RDKit canonicalization")
        monkeypatch.setattr(graphs, "canonical_graph", unexpected)
        for reverse in (True, False):
            with pytest.raises(graphs.SharedCoreProjectionLimit, match=reason):
                _project_centers(["CCC"], reverse=reverse)


def dense_reduction():
    """Return one graph with eighteen marked centers on a connected carbon chain."""
    size = 18
    left, right = [], []
    for i in range(1, size + 1):
        end = i in {1, size}
        left.append(f"[CH{'2' if end else ''}:{i}]")
        right.append(f"[CH{'3' if end else '2'}:{i}]")
        if i % 2:
            left.append("=")
    return "".join(left) + ">>" + "".join(right)


def test_dense_observation_retains_scientific_identity_and_explains_unavailable_local_keys():
    reaction = dense_reduction()
    analysis = featurize_reaction(reaction)
    projection = build_shared_reaction_core(
        reaction, asdict(analysis.reaction_signature), asdict(analysis.reaction_core),
    )
    assert projection.reaction_identity and projection.input_identity and projection.product_identity
    assert projection.observation_core_id == analysis.reaction_core.core_id
    assert projection.levels == ()
    assert projection.unavailable_reasons == ("LOCAL_CORE_CENTER_DISTANCE_PAIR_LIMIT",)
    assert SharedReactionCore.from_dict(projection.to_dict()) == projection


def test_qualified_level_survives_unavailable_relaxed_projection(monkeypatch):
    reaction = "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1"
    analysis = featurize_reaction(reaction)
    args = reaction, asdict(analysis.reaction_signature), asdict(analysis.reaction_core)
    expected = build_shared_reaction_core(*args)
    project = graphs.observed_edit_graph
    def limited(*args, **kwargs):
        if kwargs["ports"]:
            raise graphs.SharedCoreProjectionLimit("LOCAL_CORE_INCIDENCE_VERTEX_LIMIT")
        return project(*args, **kwargs)
    monkeypatch.setattr(shared_reaction_core, "observed_edit_graph", limited)
    actual = build_shared_reaction_core(*args)
    assert actual.levels == expected.levels[:1]
    assert actual.unavailable_reasons == ("LOCAL_CORE_INCIDENCE_VERTEX_LIMIT",)
