"""Fragment topology and local graph-change evidence, including negative controls."""

from dataclasses import asdict

import pytest
from rdkit import Chem

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.fragment_search import (
    classify_fragment_embedding, compile_fragment_query, fragment_embeddings,
    fragment_search_policy, indexed_product, project_fragment_evidence,
    validate_fragment_target,
)


ETHER = "[CH3:1][Br:2].[OH:3][CH3:4]>>[CH3:1][O:3][CH3:4]"
BRIDGED_TARGET = "O=C1c2ccccc2C[C@H]3O[C@@H](C4)CC[C@@H]4N13"


@pytest.mark.parametrize("query,target,format,topology,expected", [
    (BRIDGED_TARGET, BRIDGED_TARGET, "smiles", "preserve_rings", True),
    ("O=C1CCC2OC3CCC(C3)N12", BRIDGED_TARGET, "smiles", "preserve_rings", False),
    ("O=C1CCCC2OC3CCC(C3)N12", BRIDGED_TARGET, "smiles", "preserve_rings", False),
    ("c1ccccc1", "Cc1ccccc1", "smiles", "preserve_rings", True),
    ("c1ccccc1", "c1ccc2ccccc2c1", "smiles", "preserve_rings", False),
    ("c1ccccc1", "c1ccc2ccccc2c1", "smarts", "subgraph", True),
    ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "smiles", "subgraph", False),
    ("N[C@@H](C)C(=O)O", "NC(C)C(=O)O", "smiles", "subgraph", False),
    ("CO", "COC", "smiles", "preserve_rings", True),
])
def test_target_validation_uses_corpus_matching_semantics(query, target, format, topology, expected):
    compiled = compile_fragment_query(query, format, topology)
    result = validate_fragment_target(compiled, target)
    assert result.matches_target is expected
    assert result.query_id == compiled.query_id
    assert result.schema_version == "fragment_target_validation.v1"
    assert result.compiler_version == compiled.compiler_version
    assert result.definition_version == compiled.definition_version
    assert result == validate_fragment_target(compiled, Chem.MolToSmiles(Chem.MolFromSmiles(target)))


@pytest.mark.parametrize("target", ["", "broken", "C.O", None])
def test_target_validation_rejects_invalid_targets(target):
    with pytest.raises(ValueError):
        validate_fragment_target(compile_fragment_query("C"), target)


def classify(reaction: str, query: str, *, topology: str = "preserve_rings") -> list[dict]:
    observation = asdict(featurize_reaction(reaction).observation)
    evidence = project_fragment_evidence(reaction, observation)
    compiled = compile_fragment_query(query, topology=topology)
    canonical, order = indexed_product(reaction.split(">")[-1])
    matches, _ = fragment_embeddings(compiled, Chem.MolFromSmiles(canonical))
    return [classify_fragment_embedding(compiled, match, order, evidence[0]) for match in matches]


@pytest.mark.parametrize("query,product,topology,expected", [
    ("c1ccccc1", "Cc1ccc(Cl)cc1", "preserve_rings", True),
    ("c1ccccc1", "c1ccc2ccccc2c1", "preserve_rings", False),
    ("c1ccccc1", "c1ccc2ccccc2c1", "subgraph", True),
    ("C1CCCCC1", "C1CCC2(CC1)CCCC2", "preserve_rings", False),
    ("CO", "C1COCC1", "preserve_rings", False),
    ("CO", "C1COCC1", "subgraph", True),
    ("c1ccc2c(c1)COc1ccccc1-2", "CC(=O)c1ccc2c(c1)COc1ccccc1-2", "preserve_rings", True),
    ("CO", "C[OH2+]", "subgraph", False),
    ("[13CH3]O", "CO", "subgraph", False),
    ("[2H]CO", "[2H]CO", "preserve_rings", True),
    ("[2H]CO", "[3H]CO", "preserve_rings", False),
    ("[2H]CO", "CCO", "preserve_rings", False),
    ("[nH]1cccc1", "Cn1cccc1", "preserve_rings", False),
    ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "subgraph", False),
    ("N[C@@H](C)C(=O)O", "N[C@@H](C)C(=O)O", "subgraph", True),
    ("F/C=C/F", "F/C=C\\F", "subgraph", False),
    ("[nH0]1ccccc1", "n1ccccc1", "subgraph", True),
])
def test_query_constraints(query, product, topology, expected):
    compiled = compile_fragment_query(query, topology=topology)
    matches, _ = fragment_embeddings(compiled, Chem.MolFromSmiles(product))
    assert bool(matches) is expected


@pytest.mark.parametrize("query,format", [("CC.O", "smiles"), ("[C:1]", "smiles"),
                                          ("[$(C)]", "smarts"), ("c:c", "smarts"),
                                          ("[C;R]-[O;R]", "smarts"), ("broken", "smiles"),
                                          ("[CH3]", "smiles")])
def test_invalid_or_incomplete_queries_fail_actionably(query, format):
    with pytest.raises(ValueError):
        compile_fragment_query(query, format)


def test_explicit_smarts_and_identity():
    q = compile_fragment_query("[#6]-[O,N]", "smarts", "subgraph")
    assert fragment_embeddings(q, Chem.MolFromSmiles("CN"))[0]
    assert q.query_id == compile_fragment_query("[#6]-[O,N]", "smarts", "subgraph").query_id
    assert fragment_search_policy()["schema_version"] == "fragment_search_policy.v1"


def test_internal_formation_is_not_the_same_as_boundary_attachment():
    assert any("constructed" in row["relationships"] for row in classify(ETHER, "COC"))
    rows = classify(ETHER, "CO")
    assert any("constructed" in row["relationships"] for row in rows)
    assert any("boundary_changed" in row["relationships"] and "constructed" not in row["relationships"] for row in rows)


def test_oxidation_is_modification_and_unmapped_stays_unresolved():
    rows = classify("[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]", "CC=O")
    assert rows and all("modified" in r["relationships"] and "constructed" not in r["relationships"] for r in rows)
    rows = classify("CO.CBr>>COC", "COC")
    assert rows and all(r["relationships"] == ["unresolved"] for r in rows)


def test_retention_needs_complete_local_correspondence():
    rows = classify("[CH3:1][CH2:2][CH2:3][OH:4]>>[CH3:1][CH2:2][CH:3]=[O:4]", "CC")
    assert any(r["relationships"] == ["carried_through"] for r in rows)
    rows = classify("[CH3:1]C[OH:3]>>[CH3:1]C[OH:3]", "CCO")
    assert all("unresolved" in r["relationships"] for r in rows)


def test_conflicting_and_invalid_maps_never_claim_construction():
    observation = asdict(featurize_reaction(ETHER).observation)
    observation["warnings"] = ["MAPPING_EVIDENCE_CONFLICT"]
    assert project_fragment_evidence(ETHER, observation)[0]["evidence_status"] == "unresolved"
    observation["warnings"] = []
    observation["edit_hypotheses"] = [{"hypothesis_id": "alternative"}]
    assert project_fragment_evidence(ETHER, observation)[0]["evidence_status"] == "unresolved"
    duplicate = "[CH3:1][Br:2].[OH:1][CH3:4]>>[CH3:1][O:1][CH3:4]"
    observation["edit_hypotheses"] = []
    observation["input_reaction_smiles"] = duplicate
    assert project_fragment_evidence(duplicate, observation)[0]["evidence_status"] == "unresolved"


def test_canonical_atom_order_preserves_map_identity():
    original = "[OH:8][CH2:9][CH3:10]"
    canonical, order = indexed_product(original)
    assert order != tuple(range(3))
    before, after = Chem.MolFromSmiles(original), Chem.MolFromSmiles(canonical)
    for i, original_index in enumerate(order):
        assert before.GetAtomWithIdx(original_index).GetSymbol() == after.GetAtomWithIdx(i).GetSymbol()


def test_nonquery_chord_is_not_a_construction_witness():
    q = compile_fragment_query("CCC", topology="subgraph")
    evidence = {"evidence_status": "validated_supplied_mapping", "uncertain_atoms": [],
                "changes": [{"kind": "formed", "product_atoms": [0, 2]}]}
    result = classify_fragment_embedding(q, (0, 1, 2), (0, 1, 2), evidence)
    assert "constructed" not in result["relationships"]
    assert "modified" in result["relationships"]


def test_embedding_cap_is_explicit():
    q = compile_fragment_query("C", topology="subgraph")
    matches, truncated = fragment_embeddings(q, Chem.MolFromSmiles("CCCC"), maximum=2)
    assert len(matches) == 2 and truncated


def test_mapped_ring_closure_supplies_an_internal_bond_witness():
    rows = classify("[Br:1][CH2:2][CH2:3][OH:4]>>[CH2:2]1[CH2:3][O:4]1", "C1CO1")
    assert rows and any("constructed" in r["relationships"] for r in rows)
    witnesses = [w for r in rows for w in r["witnesses"] if w["relationship"] == "constructed"]
    assert all(w["before"]["state_kind"] == "no_bond" for w in witnesses)


def test_empty_product_is_excluded():
    with pytest.raises(ValueError, match="Invalid product"):
        indexed_product("")
