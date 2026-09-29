"""Target-derived query validity, provenance, chemistry context and bounds."""

import pytest
from rdkit import Chem

from reactive_taxonomy.fragment_search import compile_fragment_query, fragment_embeddings
from reactive_taxonomy.search_fragments import search_fragment_policy, suggest_search_fragments


TARGET = "CC(=O)c1ccc2c(c1)COc1ccccc1-2"
CORE = "c1ccc2c(c1)COc1ccccc1-2"


@pytest.mark.parametrize("smiles", [
    TARGET, "c1ccc2ccccc2c1", "CC1CCC2(CC1)CCCC2", "CC1CC2CCC1C2",
    "CN1CCCC1=O", "CCOc1ccc(-c2ncccc2)cc1", "CC(=O)NCCO",
    "CCC[C@H](O)CCCN", "CC/C=C/CCO", "C[N+](C)(C)CCO",
    "[13CH3]COCCO", "[2H]C(O)CCN", "CC[n+]1ccccc1", "C", "c1ccccc1-c1ccccc1",
])
def test_all_candidates_are_exact_target_regions(smiles):
    result = suggest_search_fragments(smiles)
    target = Chem.MolFromSmiles(result.target_smiles)
    original = Chem.MolFromSmiles(smiles)
    assert 1 <= len(result.candidates) <= 5
    for reference in result.target_atoms:
        atom = original.GetAtomWithIdx(reference.input_atom_index)
        assert (reference.element, reference.isotope, reference.charge) == (
            atom.GetSymbol(), atom.GetIsotope(), atom.GetFormalCharge())
    for candidate in result.candidates:
        query = compile_fragment_query(candidate.query, candidate.query_format, candidate.topology)
        matches, _ = fragment_embeddings(query, target, maximum=10000)
        assert candidate.query_atom_target_ids in matches
        selected = set(candidate.target_atom_ids)
        assert set(candidate.query_atom_target_ids) == selected
        assert selected.isdisjoint(candidate.omitted_atom_ids)
        assert selected | set(candidate.omitted_atom_ids) == set(range(target.GetNumAtoms()))
        for boundary in candidate.boundaries:
            assert boundary.retained_atom_id in selected
            assert boundary.omitted_atom_id not in selected
            assert str(target.GetBondBetweenAtoms(
                boundary.retained_atom_id, boundary.omitted_atom_id).GetBondType()) == boundary.bond_type
        # Every retained ring bond neighbourhood is complete, including spiro/bridged rings.
        for atom_id in selected:
            for bond in target.GetAtomWithIdx(atom_id).GetBonds():
                if bond.IsInRing():
                    assert bond.GetOtherAtomIdx(atom_id) in selected
    assert result.to_dict()["candidates"][0]["target_atom_ids"] == list(result.candidates[0].target_atom_ids)


def test_cyclic_ether_core_is_first_and_context_is_not_an_invented_aldehyde():
    result = suggest_search_fragments(TARGET)
    assert result.candidates[0].query == CORE
    assert result.candidates[0].kind == "ring_system"
    assert len(result.candidates[0].boundaries) == 1
    assert any(c.query == result.target_smiles for c in result.candidates)
    for candidate in result.candidates:
        if "=O" in candidate.query:
            assert "CC(=O)" in candidate.query


def test_carbonyl_stereo_charge_and_isotope_context_survive():
    lactam = suggest_search_fragments("CCC1CCNC1=O").candidates[0]
    assert any(b.GetBondTypeAsDouble() == 2 and
               8 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())
               for b in Chem.MolFromSmiles(lactam.query).GetBonds())
    stereo = suggest_search_fragments("CC[C@H](O)CCCCCCN")
    assert any("@" in c.query for c in stereo.candidates)
    for candidate in stereo.candidates:
        if "@" in candidate.query:
            opposite = Chem.MolFromSmiles(candidate.query)
            for atom in opposite.GetAtoms():
                atom.InvertChirality()
            assert not fragment_embeddings(compile_fragment_query(candidate.query), opposite)[0]
    charged = suggest_search_fragments("C[N+](C)(C)CCO")
    assert any("[N+]" in c.query for c in charged.candidates)
    assert any("[13CH3]" in c.query for c in suggest_search_fragments("[13CH3]COCCN").candidates)


def test_common_ring_caution_is_not_an_invented_frequency_estimate():
    result = suggest_search_fragments("c1ccccc1-c1ccccc1")
    assert all("LOW_STRUCTURAL_SPECIFICITY_MAY_BE_BROAD" in c.cautions for c in result.candidates)
    assert all(c.structural_priority == 0 for c in result.candidates)
    assert "not rarity" in result.limitations[1]


def test_canonical_identity_and_original_input_projection():
    mol = Chem.MolFromSmiles(TARGET)
    reversed_mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    alternate = Chem.MolToSmiles(reversed_mol, canonical=False)
    for i, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetAtomMapNum(i)
    mapped = Chem.MolToSmiles(mol, canonical=False)
    baseline = suggest_search_fragments(TARGET)
    for value in (alternate, mapped):
        result = suggest_search_fragments(value)
        assert result.target_id == baseline.target_id
        assert result.candidates == baseline.candidates
        original = Chem.MolFromSmiles(value)
        target = Chem.MolFromSmiles(result.target_smiles)
        order = {r.atom_id: r.input_atom_index for r in result.target_atoms}
        for bond in target.GetBonds():
            original_bond = original.GetBondBetweenAtoms(order[bond.GetBeginAtomIdx()], order[bond.GetEndAtomIdx()])
            assert original_bond.GetBondType() == bond.GetBondType()


def test_explicit_selection_extracts_exactly_and_rejects_conflicting_context():
    suggested = suggest_search_fragments(TARGET)
    candidate = suggested.candidates[0]
    result = suggest_search_fragments(TARGET, selected_atom_ids=list(candidate.target_atom_ids))
    assert result.selection_mode == "selected_atoms"
    assert len(result.candidates) == 1
    assert result.candidates[0].query == candidate.query
    assert result.candidates[0].target_atom_ids == candidate.target_atom_ids
    with pytest.raises(ValueError, match="complete ring"):
        suggest_search_fragments(TARGET, selected_atom_ids=list(candidate.target_atom_ids[:-1]))
    with pytest.raises(ValueError, match="connected"):
        suggest_search_fragments("CCCO", selected_atom_ids=[0, 3])
    with pytest.raises(ValueError, match="valence/stereo"):
        suggest_search_fragments("CC(=O)C", selected_atom_ids=[1, 2])


def test_equivalent_explicit_occurrences_have_distinct_provenance_ids():
    mol = Chem.MolFromSmiles(suggest_search_fragments("c1ccccc1-c1ccccc1").target_smiles)
    first, second = mol.GetRingInfo().AtomRings()
    a = suggest_search_fragments(Chem.MolToSmiles(mol), selected_atom_ids=list(first)).candidates[0]
    b = suggest_search_fragments(Chem.MolToSmiles(mol), selected_atom_ids=list(second)).candidates[0]
    assert a.query == b.query
    assert a.candidate_id != b.candidate_id


def test_extraction_must_not_silently_erase_stereo_when_caps_make_arms_equivalent():
    smiles = "CC[C@H](O)CCCCCCN"
    result = suggest_search_fragments(smiles)
    target = Chem.MolFromSmiles(result.target_smiles)
    center = next(a for a in target.GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
    selected = [center.GetIdx()] + [a.GetIdx() for a in center.GetNeighbors()]
    with pytest.raises(ValueError, match="loses specified stereochemistry"):
        suggest_search_fragments(smiles, selected_atom_ids=selected)
    assert result.rejected_count > 0
    assert any(c.query == result.target_smiles for c in result.candidates)


@pytest.mark.parametrize("smiles,kwargs", [
    ("", {}), ("bad", {}), ("C.O", {}), ("C*", {}), ("[CH3]", {}), ("C" * 201, {}),
    (TARGET, {"limit": True}), (TARGET, {"limit": 6}), (TARGET, {"limit": 0}),
    (TARGET, {"selected_atom_ids": []}), (TARGET, {"selected_atom_ids": [1, 1]}),
    (TARGET, {"selected_atom_ids": [True]}), (TARGET, {"selected_atom_ids": [1000]}),
])
def test_invalid_or_ambiguous_input_is_rejected(smiles, kwargs):
    with pytest.raises(ValueError):
        suggest_search_fragments(smiles, **kwargs)


def test_bounds_are_explicit_and_small_molecules_do_not_force_five_candidates():
    assert len(suggest_search_fragments("CO").candidates) == 1
    limited = suggest_search_fragments(TARGET, limit=1)
    assert len(limited.candidates) == 1 and limited.output_truncated
    oversized = suggest_search_fragments("C" * 101)
    assert not oversized.candidates and oversized.rejected_count == 1
    truncated = suggest_search_fragments("CO" * 90)
    assert truncated.generation_truncated
    assert len(truncated.candidates) <= 5


def test_policy_schema_and_weights_are_validated(monkeypatch):
    import reactive_taxonomy.search_fragments as module
    policy = search_fragment_policy()
    assert policy["definition_version"] == "search_fragments.v1@1.0"
    policy["feature_weights"]["heteroatoms"] = -1
    monkeypatch.setattr(module.json, "loads", lambda _: policy)
    with pytest.raises(ValueError, match="nonnegative"):
        search_fragment_policy()
