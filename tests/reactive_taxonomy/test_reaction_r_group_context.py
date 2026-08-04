"""Functional-group overlays for graph-provenanced R subgraphs."""

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction


PIPERAZINE_SUBSTITUTION = (
    "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O."
    "CC(C)c1nc(Cl)nc(Cl)c1Br>>"
    "CC(C)c1nc(Cl)nc(N2CCN(c3cc4c(cc3F)c(=O)c(C(=O)O)"
    "cn4C3CC3)CC2)c1Br"
)


def _shared_piperazine_context():
    analysis = featurize_reaction(PIPERAZINE_SUBSTITUTION)
    assert analysis.interpretation is not None
    context = next(
        value
        for value in analysis.interpretation.r_group_functional_contexts
        if value.component_index == 0
        and value.remote_class == "ring_aliphatic"
    )
    return analysis, context


def test_shared_r_scaffold_contains_unchanged_functional_groups() -> None:
    analysis, context = _shared_piperazine_context()

    assert analysis.interpretation is not None
    assert analysis.interpretation.schema_version == "7.1"
    assert context.context_id.startswith("RGFC1:")
    assert context.continuity == "retained"
    assert len(context.attachment_profile_ids) == 1
    assert [group.motif_id for group in context.functional_groups] == [
        "tertiary_amine",
        "aryl_halide",
        "carboxylic_acid",
    ]
    assert all(
        set(group.atom_indices).issubset(context.remote_atom_indices)
        for group in context.functional_groups
    )


def test_shared_r_scaffold_retains_each_port_distance() -> None:
    _, context = _shared_piperazine_context()
    groups = {group.motif_id: group for group in context.functional_groups}

    assert [
        distance.bond_distance
        for distance in groups["tertiary_amine"].port_distances
    ] == [2, 2]
    assert groups["tertiary_amine"].distance_to_reactive_site == 3
    assert [
        distance.bond_distance
        for distance in groups["carboxylic_acid"].port_distances
    ] == [9, 9]
    assert groups["carboxylic_acid"].distance_to_reactive_site == 10


def test_consumed_handle_is_not_an_r_group_functional_context() -> None:
    analysis = featurize_reaction(
        "Brc1ccc(C#N)cc1.CN>>CNc1ccc(C#N)cc1"
    )
    assert analysis.interpretation is not None
    motif_ids = {
        group.motif_id
        for context in analysis.interpretation.r_group_functional_contexts
        for group in context.functional_groups
    }

    assert "nitrile" in motif_ids
    assert "aryl_halide" not in motif_ids


def test_r_group_functional_context_serializes_only_with_interpretation() -> None:
    analysis, context = _shared_piperazine_context()
    payload = analysis.to_dict()

    assert "r_group_functional_contexts" not in asdict(analysis.observation)
    serialized = payload["interpretation"]["r_group_functional_contexts"]
    assert serialized
    assert serialized[0]["context_id"].startswith("RGFC1:")
    assert any(
        value["remote_subgraph_id"] == context.remote_subgraph_id
        for value in serialized
    )
