import pytest

from reactive_taxonomy.reaction_connectivity import (
    bond_state,
    build_connectivity_edit_graph,
    endpoint_absent_state,
    make_bond_transition,
    reaction_edits_from_connectivity_graph,
)
from reactive_taxonomy.reaction_models import (
    BondState,
    BondTransition,
    ReactionAtomReference,
)


def _atom(index: int, *, element: str = "C") -> ReactionAtomReference:
    return ReactionAtomReference(
        side="reactant",
        component_index=0,
        atom_index=index,
        atom_map_number=None,
        element=element,
        formal_charge=0,
        aromatic=False,
        hybridization="SP3",
        local_environment_id=f"test:{element}:{index}",
    )


def _formed(left: ReactionAtomReference, right: ReactionAtomReference):
    return make_bond_transition(
        atom_1=left,
        atom_2=right,
        before_state=bond_state(None),
        after_state=bond_state("SINGLE"),
        observation_scope="observed_product",
        evidence="test",
        confidence=1.0,
    )


def test_bond_state_rejects_contradictory_order() -> None:
    with pytest.raises(ValueError, match="cannot carry"):
        BondState("no_bond", "SINGLE")

    with pytest.raises(ValueError, match="requires an order"):
        BondState("bond", None)


def test_bond_transition_validates_delta_and_projection_scope() -> None:
    left = _atom(0)
    right = _atom(1)

    with pytest.raises(ValueError, match="must match"):
        BondTransition(
            atom_1=left,
            atom_2=right,
            before_state=bond_state("DOUBLE"),
            after_state=bond_state("SINGLE"),
            delta_units=1,
            observation_scope="observed_product",
            evidence="test",
            confidence=1.0,
        )

    with pytest.raises(ValueError, match="main_product_projection"):
        BondTransition(
            atom_1=left,
            atom_2=right,
            before_state=bond_state("SINGLE"),
            after_state=endpoint_absent_state(),
            delta_units=None,
            observation_scope="observed_product",
            evidence="test",
            confidence=1.0,
        )


def test_projected_attachment_has_no_delta_but_keeps_legacy_broken_view() -> None:
    transition = make_bond_transition(
        atom_1=_atom(0),
        atom_2=_atom(1, element="Br"),
        before_state=bond_state("SINGLE"),
        after_state=endpoint_absent_state(),
        observation_scope="main_product_projection",
        evidence="supplied_atom_mapping",
        confidence=1.0,
    )
    graph = build_connectivity_edit_graph(
        bond_transitions=(transition,),
        evidence="validated_atom_mapping",
        confidence=1.0,
    )

    assert transition.delta_units is None
    edits = reaction_edits_from_connectivity_graph(graph)
    assert len(edits) == 1
    assert edits[0].edit_type == "broken"
    assert edits[0].old_order == "SINGLE"
    assert edits[0].new_order is None


def test_canonical_graph_distinguishes_edit_incidence_not_only_counts() -> None:
    atoms = tuple(_atom(index) for index in range(4))
    path = build_connectivity_edit_graph(
        bond_transitions=(
            _formed(atoms[0], atoms[1]),
            _formed(atoms[1], atoms[2]),
            _formed(atoms[2], atoms[3]),
        ),
        evidence="test",
        confidence=1.0,
    )
    star = build_connectivity_edit_graph(
        bond_transitions=(
            _formed(atoms[0], atoms[1]),
            _formed(atoms[0], atoms[2]),
            _formed(atoms[0], atoms[3]),
        ),
        evidence="test",
        confidence=1.0,
    )

    assert path.shadow_key != star.shadow_key


def test_aromatic_transition_does_not_receive_integer_arithmetic() -> None:
    transition = make_bond_transition(
        atom_1=_atom(0),
        atom_2=_atom(1),
        before_state=bond_state("AROMATIC"),
        after_state=bond_state("SINGLE"),
        observation_scope="observed_product",
        evidence="test",
        confidence=1.0,
    )
    graph = build_connectivity_edit_graph(
        bond_transitions=(transition,),
        evidence="test",
        confidence=1.0,
    )

    assert transition.delta_units is None
    assert "UNSUPPORTED_BOND_DOMAIN" in graph.warnings


def test_large_symmetric_component_has_explicit_deterministic_overflow() -> None:
    center = _atom(0)
    leaves = tuple(_atom(index) for index in range(1, 10))
    transitions = tuple(_formed(center, leaf) for leaf in leaves)

    forward = build_connectivity_edit_graph(
        bond_transitions=transitions,
        evidence="test",
        confidence=1.0,
    )
    reversed_order = build_connectivity_edit_graph(
        bond_transitions=tuple(reversed(transitions)),
        evidence="test",
        confidence=1.0,
    )

    assert "CONNECTIVITY_CANONICALIZATION_OVERFLOW" in forward.warnings
    assert forward.shadow_key == reversed_order.shadow_key
