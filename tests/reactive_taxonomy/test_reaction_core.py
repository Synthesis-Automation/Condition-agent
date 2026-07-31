"""Template-free reaction-core projection regressions."""

from __future__ import annotations

from typing import Any

from reactive_taxonomy import build_reaction_core_projection, featurize_reaction
from reactive_taxonomy.reaction_core import (
    ReactionCoreProjection,
    build_reaction_core_projection as package_build_reaction_core_projection,
)


def test_reaction_core_package_preserves_the_public_constructor() -> None:
    """The focused package remains the sole implementation import surface."""
    assert (
        package_build_reaction_core_projection
        is build_reaction_core_projection
    )
    assert ReactionCoreProjection.__name__ == "ReactionCoreProjection"


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
MAPPED_DECARBOXYLATIVE_COUPLING = (
    "O=[C:15](O)[CH2:2][CH3:1]."
    "[cH:3]1[cH:4][c:5]([CH3:6])[n:7][c:8]2[cH:9][cH:10]"
    "[c:11]([F:12])[cH:13][c:14]12"
    ">>[CH3:1][CH2:2][c:3]1[cH:4][c:5]([CH3:6])[n:7]"
    "[c:8]2[cH:9][cH:10][c:11]([F:12])[cH:13][c:14]12"
)
MAPPED_SECONDARY_DECARBOXYLATIVE_COUPLING = (
    "[CH3:1][CH:2]([CH3:3])[C:4](=O)O."
    "[cH:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1"
    ">>[CH3:1][CH:2]([CH3:3])[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1"
)
MAPPED_ARYL_DECARBOXYLATIVE_COUPLING = (
    "O=[C:13](O)[c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1."
    "[cH:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1"
    ">>[c:1]1(-[c:7]2[cH:8][cH:9][cH:10][cH:11][cH:12]2)"
    "[cH:2][cH:3][cH:4][cH:5][cH:6]1"
)
MAPPED_INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION = (
    "O=[C:19](O)[C:2](=[O:1])[c:3]1[cH:4][c:5]"
    "([C:6]([F:7])([F:8])[F:9])[cH:10][cH:11][c:12]1-"
    "[c:13]1[cH:14][cH:15][cH:16][n:17][cH:18]1"
    ">>[O:1]=[C:2]1[c:3]2[cH:4][c:5]"
    "([C:6]([F:7])([F:8])[F:9])[cH:10][cH:11][c:12]2-"
    "[c:13]2[cH:14][cH:15][cH:16][n:17][c:18]21"
)
MAPPED_HETEROARYL_AMIDE_FORMATION = (
    "[CH3:1][C:2]1([CH3:3])[CH2:4][CH2:5][CH:6]([NH2:7])"
    "[CH2:8][CH2:9]1."
    "[O:10]=[C:11]([OH:12])[c:13]1[cH:14][cH:15][cH:16][nH:17]1"
    ">>"
    "[CH3:1][C:2]1([CH3:3])[CH2:4][CH2:5][CH:6]"
    "([NH:7][C:11](=[O:10])[c:13]2[cH:14][cH:15][cH:16][nH:17]2)"
    "[CH2:8][CH2:9]1"
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
    assert core.core_id.startswith("RCP2:")
    assert core.exact_core_key.startswith("RCX2:")
    assert core.typed_core_key.startswith("RCT2:")
    assert core.shape_core_key.startswith("RSH2:")
    assert core.center_transition_key.startswith("RCS2:")
    assert core.generic_label == (
        "C(H)(Ar)(=O) + O(H)(R) → C(H)(Ar)(O-R)2"
    )
    primary = [
        transition
        for transition in core.atom_transitions
        if transition.role == "primary_center"
    ]
    assert len(primary) == 1
    assert primary[0].atom_map_number == 3
    assert primary[0].before_state is not None
    assert primary[0].after_state is not None
    assert len(core.atom_transitions) > len(primary)


def test_acetal_remote_subgraph_is_derived_from_removed_fluoroaryl_graph() -> None:
    core = featurize_reaction(RXNMAPPER_STYLE_ACETAL).reaction_core

    assert core is not None
    aryl_subgraphs = [
        subgraph
        for subgraph in core.remote_subgraphs
        if subgraph.remote_class == "aryl"
    ]
    assert len(aryl_subgraphs) == 2
    assert {subgraph.side for subgraph in aryl_subgraphs} == {
        "reactant",
        "product",
    }
    assert all(
        subgraph.continuity == "retained"
        for subgraph in aryl_subgraphs
    )
    assert all(
        subgraph.fragment_smiles == "Fc1ccccc1"
        for subgraph in aryl_subgraphs
    )
    assert all(
        subgraph.functional_group_ids == ("aryl_halide",)
        for subgraph in aryl_subgraphs
    )
    assert all(len(subgraph.attachment_ports) == 1 for subgraph in aryl_subgraphs)


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
    assert mapper_core.shape_core_key == curated_core.shape_core_key
    assert (
        mapper_core.center_transition_key
        == curated_core.center_transition_key
    )
    assert mapper_core.generic_label == curated_core.generic_label
    assert mapper_core.warnings == ()
    assert curated_core.warnings == ()


def test_single_center_label_includes_an_external_formed_bond_partner() -> None:
    forward = featurize_reaction(
        MAPPED_HETEROARYL_AMIDE_FORMATION
    ).reaction_core
    reversed_reactants = featurize_reaction(
        MAPPED_HETEROARYL_AMIDE_FORMATION.replace(
            MAPPED_HETEROARYL_AMIDE_FORMATION.split(">>")[0],
            ".".join(
                reversed(
                    MAPPED_HETEROARYL_AMIDE_FORMATION.split(">>")[0].split(".")
                )
            ),
            1,
        )
    ).reaction_core

    assert forward is not None
    assert reversed_reactants is not None
    expected = (
        "C(HetAr)(O-H)(=O) + N(H)2(Cycloalkyl) "
        "→ C(HetAr)(N-R)(=O)"
    )
    assert forward.generic_label == expected
    assert reversed_reactants.generic_label == expected


def test_removed_heteroaryl_and_alkyl_graphs_receive_distinct_types() -> None:
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
        subgraph.remote_class for subgraph in heteroaryl.remote_subgraphs
    } == {"heteroaryl"}
    assert heteroaryl.generic_label == (
        "C(H)(HetAr)(=O) → C(H)2(HetAr)(O-H)"
    )
    assert {
        subgraph.remote_class for subgraph in alkyl.remote_subgraphs
    } == {"alkyl"}
    assert alkyl.generic_label == (
        "C(R)2(=O) → C(H)(R)2(O-H)"
    )


def test_core_label_exposes_a_departing_carboxyl_neighbor() -> None:
    core = featurize_reaction(MAPPED_DECARBOXYLATIVE_COUPLING).reaction_core

    assert core is not None
    assert core.generic_label == (
        "C(H)2(R)(C(=O)(O-H)) + C(H)(:HetAr)2 → C(H)2(R)(ArC)"
    )
    assert any(
        transition.before_state is not None
        and transition.before_state.concise_label == "C(R)(O)(=O)"
        and transition.after_state is None
        for transition in core.atom_transitions
    )


def test_decarboxylative_abstraction_separates_motif_and_limiters() -> None:
    primary = featurize_reaction(
        MAPPED_DECARBOXYLATIVE_COUPLING
    ).reaction_core
    secondary = featurize_reaction(
        MAPPED_SECONDARY_DECARBOXYLATIVE_COUPLING
    ).reaction_core
    aryl = featurize_reaction(
        MAPPED_ARYL_DECARBOXYLATIVE_COUPLING
    ).reaction_core

    assert primary is not None and primary.abstraction is not None
    assert secondary is not None and secondary.abstraction is not None
    assert aryl is not None and aryl.abstraction is not None
    abstractions = (
        primary.abstraction,
        secondary.abstraction,
        aryl.abstraction,
    )
    assert {value.motif_key for value in abstractions} == {
        primary.abstraction.motif_key
    }
    assert {value.general_label for value in abstractions} == {
        "R–C(=O)OH + Ar–H → R–Ar"
    }
    assert "transfer_center:primary_alkyl" in (
        primary.abstraction.limiter_tokens
    )
    assert "transfer_center:secondary_alkyl" in (
        secondary.abstraction.limiter_tokens
    )
    assert "transfer_center:aryl" in aryl.abstraction.limiter_tokens
    assert primary.abstraction.limiter_label == (
        "R = R′–CH₂– (primary alkyl); Ar = HetAr"
    )


def test_decarboxylative_abstraction_respects_intramolecular_topology() -> None:
    result = featurize_reaction(
        MAPPED_INTRAMOLECULAR_DECARBOXYLATIVE_CYCLIZATION
    )
    core = result.reaction_core

    assert result.reaction_topology is not None
    assert result.reaction_topology.reaction_scope == "intramolecular"
    assert result.reaction_topology.formed_ring_sizes == (5,)
    assert core is not None and core.abstraction is not None
    assert core.abstraction.general_label == (
        "R–C(=O)OH + Ar–H → R–Ar; intramolecular, "
        "5-membered ring"
    )
    assert core.abstraction.limiter_label == (
        "acyl center = R′–C(=O)– (acyl); partner center = HetAr"
    )
    assert "topology:intramolecular_cyclization" in (
        core.abstraction.limiter_tokens
    )
    assert "transfer_center:acyl" in core.abstraction.limiter_tokens


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
    assert forward.shape_core_key == reversed_order.shape_core_key


def test_shape_key_separates_same_center_transition_with_different_handles() -> None:
    vinyl_suzuki = featurize_reaction(
        "[Br:9][c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1."
        "O[B:10](O)[CH:2]=[CH2:1]"
        ">>"
        "[CH2:1]=[CH:2][c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1"
    ).reaction_core
    heck = featurize_reaction(
        "[Br:9][c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1."
        "[CH2:1]=[CH2:2]"
        ">>"
        "[CH2:1]=[CH:2][c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1"
    ).reaction_core

    assert vinyl_suzuki is not None
    assert heck is not None
    assert (
        vinyl_suzuki.center_transition_key
        == heck.center_transition_key
    )
    assert vinyl_suzuki.shape_core_key != heck.shape_core_key
    assert "reactant:site:transfer_group" in vinyl_suzuki.participant_tokens
    assert "reactant:site:transfer_group" not in heck.participant_tokens


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

    assert mapped_payload["schema_version"] == "3.5"
    assert mapped_payload["reaction_core"]["schema_version"] == "2.1"
    assert mapped_payload["reaction_core"]["algorithm_version"] == (
        "reaction_core_projection.v7"
    )
    assert unmapped.reaction_core is None
