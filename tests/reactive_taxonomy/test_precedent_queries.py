"""Scaffold queries broaden substitution explicitly without inventing molecules."""

import pytest
from rdkit import Chem

from reactive_taxonomy.fragment_search import compile_fragment_query, validate_fragment_target
from reactive_taxonomy.precedent_queries import plan_precedent_queries
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts

TARGET = "CC1(CC(C)(C)C(c(c1n2C(C)=O)c3c2ccnc3)=O)C"
NH_CORE = "[H][n]1c2c(cncc2)c2C(CCCc12)=O"


def test_core_accepts_nh_and_acetyl_and_refinements_remain_nested():
    plan = plan_precedent_queries(TARGET)
    ladder = plan["ladders"][0]
    assert 2 <= len(ladder) <= 4
    first = compile_fragment_query(ladder[0]["query"], "smarts", "subgraph")
    assert validate_fragment_target(first, NH_CORE).matches_target
    for before, after in zip(ladder, ladder[1:]):
        assert set(before["target_atom_ids"]) < set(after["target_atom_ids"])
    for step in ladder:
        query = compile_fragment_query(step["query"], "smarts", "subgraph")
        assert validate_fragment_target(query, TARGET).matches_target
    assert not validate_fragment_target(compile_fragment_query(ladder[-1]["query"], "smarts", "subgraph"), NH_CORE).matches_target
    assert "unconstrained_hydrogen_count" in plan["relaxations"]


def test_serialization_invariance_and_retained_stereo_charge():
    target = "C[C@H](O)c1cc[nH+]cc1"
    plan = plan_precedent_queries(target)
    assert plan == plan_precedent_queries(Chem.MolToSmiles(Chem.MolFromSmiles(target), rootedAtAtom=4))
    step = plan["ladders"][0][-1]
    query = compile_fragment_query(step["query"], "smarts", "subgraph")
    assert validate_fragment_target(query, target).matches_target
    assert not validate_fragment_target(query, "C[C@@H](O)c1cc[nH+]cc1").matches_target
    assert not validate_fragment_target(query, "C[C@H](O)c1ccncc1").matches_target


@pytest.mark.parametrize("target", ["bad", "C.C", "[CH2]C", "*CC"])
def test_invalid_ambiguous_targets_fail_explicitly(target):
    with pytest.raises(ValueError):
        plan_precedent_queries(target)


def test_smarts_bond_order_is_not_an_alternative_bond():
    assert compile_fragment_query("[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-1", "smarts")
    with pytest.raises(ValueError, match="Alternative bond"):
        compile_fragment_query("C-,=C", "smarts")


@pytest.mark.parametrize("target", [
    "CC(C)C[C@@H]1CN2[C@@H](c3c(C2)cccc3)CC1=O",
    "CC(C)C[C@H]1CN2[C@H](c3c(C2)cccc3)CC1=O",
    "C[C@]1(O)CCC[C@@H]1F",
    "[2H][C@](F)(Cl)Br",
    "[13CH3][C@@H](O)c1ccccc1",
])
def test_stereo_ring_queries_survive_serialization_and_reject_inversion(target):
    molecule = Chem.MolFromSmiles(target)
    plan = plan_precedent_queries(target)
    assert plan["definition_version"] == "precedent_discovery.v1@1.2"
    assert plan == plan_precedent_queries(Chem.MolToSmiles(molecule, rootedAtAtom=4))
    canonical = Chem.MolFromSmiles(plan["target_smiles"])
    for ladder in plan["ladders"]:
        for step in ladder:
            query = compile_fragment_query(step["query"], "smarts", "subgraph")
            assert validate_fragment_target(query, target).matches_target
            assert step["target_atom_ids"] in canonical.GetSubstructMatches(
                query.molecule, useChirality=True, uniquify=False,
            )
        full = compile_fragment_query(ladder[-1]["query"], "smarts", "subgraph")
        for atom in molecule.GetAtoms():
            if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            inverted = Chem.Mol(molecule)
            inverted.GetAtomWithIdx(atom.GetIdx()).InvertChirality()
            assert not validate_fragment_target(full, Chem.MolToSmiles(inverted)).matches_target
    # Query creation must not inject stereo into shared cached atom templates.
    for expression in ("[#6;A;+0]", "[#6;A;+0;H1]", "[C;H1;+0]"):
        assert compile_smarts(expression).GetAtomWithIdx(0).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED


def test_focused_context_retains_ester_position_without_other_substituents():
    target = "COc1cc(C(=O)OCC)ccn1"
    plan = plan_precedent_queries(target)
    focused = plan["focused_queries"]
    assert 1 <= len(focused) <= 3
    ester_queries = [compile_fragment_query(step["query"], "smarts", "subgraph")
                     for step in focused if validate_fragment_target(
                         compile_fragment_query(step["query"], "smarts", "subgraph"),
                         "CCOC(=O)c1ccncc1").matches_target]
    assert ester_queries
    for query in ester_queries:
        assert validate_fragment_target(query, target).matches_target
        assert not validate_fragment_target(query, "CCOC(=O)c1ccccn1").matches_target
        assert not validate_fragment_target(query, "CCOC(=O)c1ccccc1").matches_target
    assert plan == plan_precedent_queries(Chem.MolToSmiles(Chem.MolFromSmiles(target), rootedAtAtom=3))


def test_focused_context_keeps_stereo_and_explicit_atom_correspondence():
    target = "COc1cc([C@H](O)C)ccn1"
    plan = plan_precedent_queries(target)
    molecule = Chem.MolFromSmiles(plan["target_smiles"])
    stereo_queries = []
    for step in plan["focused_queries"]:
        query = compile_fragment_query(step["query"], "smarts", "subgraph")
        assert step["target_atom_ids"] in molecule.GetSubstructMatches(
            query.molecule, useChirality=True, uniquify=False)
        if any(molecule.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
               for i in step["target_atom_ids"]):
            stereo_queries.append(query)
    assert stereo_queries
    assert all(not validate_fragment_target(query, "COc1cc([C@@H](O)C)ccn1").matches_target
               for query in stereo_queries)


@pytest.mark.parametrize("key,value", [
    ("max_focused_queries_per_core", 0), ("max_focused_queries_per_core", 4),
    ("focused_context_radius", True), ("focused_context_radius", 4),
    ("focused_context_priority", ["atom_count"]),
])
def test_focused_policy_rejects_invalid_definitions(monkeypatch, key, value):
    import json
    from reactive_taxonomy.precedent_queries import discovery_policy

    policy = discovery_policy()
    policy[key] = value
    monkeypatch.setattr(json, "loads", lambda _: policy)
    with pytest.raises(ValueError, match="policy"):
        discovery_policy()
