"""Graph site identity without invented correspondence or fragment chemistry."""

from rdkit import Chem
import pytest

from reactive_taxonomy.reaction_identity_graphs import (
    provisional_product_site,
    retained_precursor_graph,
)


@pytest.mark.parametrize(
    "reaction",
    [
        "[CH3:1][Br:2].[NH3:3]>>[CH3:1][NH2:3]",
        "[CH3:1][Br:2].[OH2:3]>>[CH3:1][OH:3]",
        "[CH3:1][Br:2].[SH2:3]>>[CH3:1][SH:3]",
        "[CH3:1][CH:2]=[O:3]>>[CH3:1][CH2:2][OH:3]",
        "[CH3:1][O:2][CH3:3]>>[CH3:1][OH:2]",
    ],
)
def test_sites_cover_formation_order_changes_and_departure(reaction):
    key = provisional_product_site(reaction)
    assert key.startswith("SITE2:")
    reactants, product = reaction.split(">>")
    shuffled = ".".join(reversed(reactants.split("."))) + ">>" + product
    assert provisional_product_site(shuffled) == key


def test_site_identity_is_map_and_serialization_invariant():
    first = "[CH3:1][O:2][CH3:3]>>[CH3:1][OH:2]"
    second = "[CH3:90][O:80][CH3:70]>>[OH:80][CH3:90]"
    assert provisional_product_site(first) == provisional_product_site(second)


def test_site_distinguishes_non_equivalent_target_positions():
    first = "[CH3:1][CH:2]([Br:8])[CH2:3][OH:4]>>[CH3:1][CH2:2][CH2:3][OH:4]"
    second = "[CH2:1]([Br:8])[CH2:2][CH2:3][OH:4]>>[CH3:1][CH2:2][CH2:3][OH:4]"
    assert provisional_product_site(first) != provisional_product_site(second)


@pytest.mark.parametrize(
    "reaction",
    [
        "CCBr.N>>CCN",  # absent correspondence
        "[CH3:1][O:1]>>[CH3:1][OH:2]",  # duplicate maps
        "[CH3:1][Br:2]>>[CH3:1][OH:2]",  # contradictory element correspondence
        "[CH3:1][OH:2]>>[CH3:1][OH:2]",  # no edits
    ],
)
def test_unusable_or_noop_mapping_does_not_create_site(reaction):
    assert provisional_product_site(reaction) == ""


def test_retained_aromatic_fragment_is_graph_identity_not_capped_molecule():
    precursor = Chem.MolFromSmiles("[CH3:9][n:1]1[cH:2][cH:3][cH:4][cH:5]1")
    product = Chem.MolFromSmiles("[nH:1]1[cH:2][cH:3][cH:4][cH:5]1")
    before = Chem.MolToSmiles(precursor)
    key = retained_precursor_graph(precursor, product, "operator")
    assert key.startswith("SYN2:")
    assert Chem.MolToSmiles(precursor) == before
    # The observed N-H state is distinct; deletion does not invent that H.
    assert retained_precursor_graph(product, product, "operator") != key
