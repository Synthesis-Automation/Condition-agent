"""Template-free reaction-core projection regressions."""

from __future__ import annotations

from typing import Any

from reactive_taxonomy import featurize_reaction


RXNMAPPER_STYLE_ACETAL = (
    "[CH3:1][OH:2].O[CH3:5]."
    "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    ">>"
    "[CH3:1][O:2][CH:3]([O:4][CH3:5])"
    "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
)

CHEMICALLY_CURATED_ACETAL = (
    "[CH3:1][OH:2].[OH:13][CH3:5]."
    "[CH:3](=[O:4])[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
    ">>"
    "[CH3:1][O:2][CH:3]([O:13][CH3:5])"
    "[c:6]1[cH:7][cH:8][cH:9][cH:10][c:11]1[F:12]"
)


def _unexpected_template_load(*args: Any, **kwargs: Any) -> None:
    raise AssertionError("reaction-core construction consulted a template")


def test_acetal_core_is_built_without_template_or_grammar_lookup(
    monkeypatch: Any,
) -> None:
    monkeypatch.setattr(
        "reactive_taxonomy.reaction_templates.load_reaction_template_registry",
        _unexpected_template_load,
    )

    result = featurize_reaction(RXNMAPPER_STYLE_ACETAL)
    core = result.reaction_core

    assert result.named_family is None
    assert core is not None
    assert core.core_id.startswith("RCP1:")
    assert core.exact_core_key.startswith("RCX1:")
    assert core.typed_core_key.startswith("RCT1:")
    assert core.generic_core_key.startswith("RCG1:")
    assert core.center_transition_key.startswith("RCS1:")
    assert core.generic_label == (
        "C(H)(Ar)(=O) → C(H)(Ar)(O-R)2"
    )
    assert len(core.centers) == 1
    assert core.centers[0].atom_map_number == 3
    assert core.centers[0].before_state is not None
    assert core.centers[0].after_state is not None


def test_acetal_boundary_is_derived_from_removed_fluoroaryl_graph() -> None:
    core = featurize_reaction(RXNMAPPER_STYLE_ACETAL).reaction_core

    assert core is not None
    aryl_boundaries = [
        boundary
        for boundary in core.boundaries
        if boundary.boundary_class == "aryl"
    ]
    assert len(aryl_boundaries) == 2
    assert {boundary.side for boundary in aryl_boundaries} == {
        "reactant",
        "product",
    }
    assert all(
        boundary.continuity == "retained"
        for boundary in aryl_boundaries
    )
    assert all(
        boundary.fragment_smiles == "Fc1ccccc1"
        for boundary in aryl_boundaries
    )
    assert all(
        boundary.functional_group_ids == ("aryl_halide",)
        for boundary in aryl_boundaries
    )


def test_center_transition_is_robust_to_acetal_oxygen_origin_mapping() -> None:
    mapper_core = featurize_reaction(
        RXNMAPPER_STYLE_ACETAL
    ).reaction_core
    curated_core = featurize_reaction(
        CHEMICALLY_CURATED_ACETAL
    ).reaction_core

    assert mapper_core is not None
    assert curated_core is not None
    assert mapper_core.exact_core_key != curated_core.exact_core_key
    assert mapper_core.typed_core_key != curated_core.typed_core_key
    assert (
        mapper_core.center_transition_key
        == curated_core.center_transition_key
    )
    assert mapper_core.generic_label == curated_core.generic_label
    assert (
        "REACTION_CORE_BOUNDARY_CONTINUITY_UNRESOLVED"
        in mapper_core.warnings
    )
    assert curated_core.warnings == ()


def test_removed_heteroaryl_and_alkyl_graphs_receive_distinct_boundaries() -> None:
    heteroaryl = featurize_reaction(
        "[CH:1](=[O:2])[c:3]1[n:4][cH:5][cH:6][cH:7][cH:8]1"
        ">>"
        "[CH2:1]([OH:2])[c:3]1[n:4][cH:5][cH:6][cH:7][cH:8]1"
    ).reaction_core
    alkyl = featurize_reaction(
        "[CH3:3][C:1](=[O:2])[CH3:4]"
        ">>"
        "[CH3:3][CH:1]([OH:2])[CH3:4]"
    ).reaction_core

    assert heteroaryl is not None
    assert alkyl is not None
    assert {
        boundary.boundary_class for boundary in heteroaryl.boundaries
    } == {"heteroaryl"}
    assert heteroaryl.generic_label == (
        "C(H)(HetAr)(=O) → C(H)2(HetAr)(O-H)"
    )
    assert {
        boundary.boundary_class for boundary in alkyl.boundaries
    } == {"alkyl"}
    assert alkyl.generic_label == (
        "C(R)2(=O) → C(H)(R)2(O-H)"
    )


def test_reaction_core_identity_is_reactant_order_invariant() -> None:
    forward = featurize_reaction(
        "[CH3:1][Br:2].[NH2:3][CH3:4]"
        ">>"
        "[CH3:1][NH:3][CH3:4]"
    ).reaction_core
    reversed_order = featurize_reaction(
        "[NH2:3][CH3:4].[CH3:1][Br:2]"
        ">>"
        "[CH3:1][NH:3][CH3:4]"
    ).reaction_core

    assert forward is not None
    assert reversed_order is not None
    assert forward.core_id == reversed_order.core_id
    assert forward.exact_core_key == reversed_order.exact_core_key
    assert (
        forward.center_transition_key
        == reversed_order.center_transition_key
    )


def test_core_survives_incomplete_unknown_family_without_signature() -> None:
    result = featurize_reaction(
        "[OH:10][C:2](=[O:1])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1"
        ">>"
        "[O:1]=[C:2]([F:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1"
    )
    core = result.reaction_core

    assert result.reaction_signature is None
    assert result.named_family is None
    assert result.transformation_class is None
    assert result.reaction_completeness is not None
    assert result.reaction_completeness.status == "incomplete"
    assert core is not None
    assert core.generic_label == (
        "C(Ar)(O-H)(=O) → C(Ar)(F)(=O)"
    )
    assert core.edit_tokens == (
        "broken:C-O:SINGLE>NONE",
        "formed:C-F:NONE>SINGLE",
    )


def test_reaction_core_serializes_but_does_not_invent_unmapped_observation() -> None:
    mapped_payload = featurize_reaction(RXNMAPPER_STYLE_ACETAL).to_dict()
    unmapped = featurize_reaction(
        "CO.CO.O=Cc1ccccc1F>>COC(OC)c1ccccc1F"
    )

    assert mapped_payload["schema_version"] == "3.4"
    assert mapped_payload["reaction_core"]["schema_version"] == "1.0"
    assert mapped_payload["reaction_core"]["algorithm_version"] == (
        "reaction_core_projection.v1"
    )
    assert unmapped.reaction_core is None
