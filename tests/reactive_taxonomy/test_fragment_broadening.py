"""Explicit query edits preserve chemistry constraints and target membership."""

from rdkit import Chem
import pytest

from reactive_taxonomy.fragment_broadening import broadening_policy, propose_fragment_queries
from reactive_taxonomy.fragment_search import compile_fragment_query, fragment_embeddings


def matches(variant, smiles):
    query = compile_fragment_query(variant["query"], variant["query_format"], variant["topology"])
    return bool(fragment_embeddings(query, Chem.MolFromSmiles(smiles))[0])


def test_explicit_ring_relaxation_preserves_aromaticity_charge_and_ambiguity():
    result = propose_fragment_queries("c1ccncc1", "c1ccncc1")
    variant = result["variants"][0]
    assert variant["relaxations"] == ["ring_boundary"]
    assert matches(variant, "c1ccc2ncccc2c1")
    assert not matches(variant, "C1CCNCC1")
    assert not matches(variant, "c1cc[nH+]cc1")
    assert variant["alignment_ambiguous"]
    assert result == propose_fragment_queries("c1ccncc1", "c1ccncc1")


def test_peripheral_edit_keeps_ring_ketone_and_does_not_require_removed_group_absent():
    target = "CC1CCC(=O)CC1"
    variants = propose_fragment_queries(target, target)["variants"]
    smaller = [v for v in variants if v["relaxations"] == ["peripheral_context"]]
    assert smaller
    for variant in smaller:
        assert matches(variant, target)
        assert matches(variant, "O=C1CCCCC1")
        assert not matches(variant, "C1CCCCC1")


def test_known_core_carbon_nitrogen_choices_recover_carbazolone_analogue():
    target = "CC1(C)CC(C)(C)c2[nH]c3ccncc3c2C1=O"
    query = "O=C1CCCc2[nH]c3ccncc3c21"
    result = propose_fragment_queries(target, query, aromatic_atom_ids=[8, 9, 10, 11])
    variant = result["variants"][-1]
    assert variant["relaxations"] == ["ring_boundary", "aromatic_carbon_nitrogen"]
    assert matches(variant, "O=C1CCCc2[nH]c3ccccc3c21")
    assert not matches(variant, "O=C1CCCc2oc3ccccc3c21")
    assert variant["target_validation"]["matches_target"]


@pytest.mark.parametrize("kwargs", [
    {"target_smiles": "CC", "query": "CO"},
    {"target_smiles": "c1cc[nH]c1", "query": "c1cc[nH]c1", "aromatic_atom_ids": [3]},
    {"target_smiles": "c1ccncc1", "query": "c1ccncc1", "aromatic_atom_ids": [True]},
    {"target_smiles": "c1ccncc1", "query": "c1ccncc1", "aromatic_atom_ids": [0, 0]},
    {"target_smiles": "c1ccncc1", "query": "c1ccncc1", "aromatic_atom_ids": [-1]},
])
def test_invalid_or_chemically_constrained_edits_rejected(kwargs):
    with pytest.raises(ValueError):
        propose_fragment_queries(**kwargs)


def test_custom_smarts_not_silently_rewritten_and_stereo_preserved():
    assert not propose_fragment_queries("CCO", "C[O,N]", "smarts", "subgraph")["variants"]
    result = propose_fragment_queries("N[C@@H](C)C(=O)O", "N[C@@H](C)C(=O)O")
    variant = result["variants"][0]
    assert matches(variant, "N[C@@H](C)C(=O)O")
    assert not matches(variant, "N[C@H](C)C(=O)O")
    assert broadening_policy()["definition_version"] == "fragment_broadening.v1@1.0"
