"""Exact route links must preserve identity, ambiguity and all components."""

import pytest

from reactive_taxonomy.reaction_sequence import connect_reaction_sequence


def test_links_ignore_serialization_maps_and_partner_order() -> None:
    steps = connect_reaction_sequence(
        (
            "CCO>>[CH3:1][CH:2]=[O:3]",
            "N.O=CC>O>CCN",
        )
    )
    assert steps[0].products == ("CC=O",)
    assert steps[1].carried_reactant_index == 1
    assert steps[1].reactants == ("N", "CC=O")
    assert steps[1].agents == ("O",)
    assert steps[0].carried_reactant_index is None


def test_unique_product_link_retains_unselected_products() -> None:
    steps = connect_reaction_sequence(("CCO>>O.CC=O", "CC=O.N>>CCN"))
    assert steps[0].product_index == 1
    assert steps[0].products == ("O", "CC=O")
    selected = connect_reaction_sequence(
        ("N.CCBr>>CCN.Br",), main_reactant_index=1, product_indices=(0,)
    )
    assert selected[0].carried_reactant_index == 1
    assert selected[0].products == ("CCN", "Br")


@pytest.mark.parametrize(
    "reactions, message",
    [
        ((), "nonempty"),
        ("CCO>>CC=O", "nonempty"),
        (("CCO>CC=O",), "expected"),
        (("CCO>>",), "expected"),
        (("C1>>C",), "Invalid SMILES"),
        (("C..O>>CO",), "Invalid SMILES"),
        (("CCO name>>CC=O",), "Invalid SMILES"),
        (("C>>N", "O>>CO"), "disconnected"),
        (("C>>N", "N.N>>NN"), "ambiguous duplicate"),
        (("C>>N.O", "N.O>>NO"), "ambiguous product"),
        (("C>>N.O",), "ambiguous product"),
        (("C>>C[NH3+]", "CN>>CNC"), "disconnected"),
        (("C>>[13CH4]", "C>>CC"), "disconnected"),
        (("C>>C[C@H](O)Cl", "C[C@@H](O)Cl>>CC(=O)Cl"), "disconnected"),
        (("C>>C[C@H](O)Cl", "CC(O)Cl>>CC(=O)Cl"), "disconnected"),
    ],
)
def test_invalid_disconnected_or_ambiguous_sequences_fail(reactions, message) -> None:
    with pytest.raises(ValueError, match=message):
        connect_reaction_sequence(reactions)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"main_reactant_index": True},
        {"main_reactant_index": 2},
        {"product_indices": ()},
        {"product_indices": (-1,)},
        {"product_indices": (False,)},
    ],
)
def test_invalid_selections_fail(kwargs) -> None:
    with pytest.raises(ValueError):
        connect_reaction_sequence(("CCO>>CC=O",), **kwargs)


def test_explicit_selection_cannot_override_disconnected_structures() -> None:
    with pytest.raises(ValueError, match="disconnected"):
        connect_reaction_sequence(("C>>N.O", "N>>CN"), product_indices=(1, None))
