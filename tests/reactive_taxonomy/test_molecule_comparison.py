"""Strict cores, bounded ambiguity and stereo are not reaction correspondence."""

from types import SimpleNamespace

import pytest
from rdkit import Chem

from reactive_taxonomy import compare_molecules
import reactive_taxonomy.molecule_comparison as module


def test_equivalent_smiles_preserve_identity_and_skip_mcs(monkeypatch):
    monkeypatch.setattr(module.rdFMCS, 'FindMCS', lambda *args: pytest.fail('Identical graph needs no MCS'))
    result = compare_molecules('[OH:9][CH2:5][CH3:3]', 'CCO')
    assert result.same_constitution
    assert result.stereo_relationship == 'same_specification'
    assert result.left_coverage == result.right_coverage == 1
    assert result.core_method == 'exact_constitution'
    assert result.alignment_count_observed == 1
    assert 'ATOM_MAP_LABELS_IGNORED_NOT_REACTION_CORRESPONDENCE' in result.warnings


def test_symmetric_aromatic_core_retains_multiple_attachments():
    result = compare_molecules('Cc1ccccc1', 'Clc1ccccc1')
    assert result.core_atom_count == 6
    assert result.alignment_count_observed == 12
    assert result.alignment_ambiguous and result.alignments_truncated
    assert not result.embeddings_truncated
    assert result.stereo_relationship == 'not_compared_different_graphs'
    assert len(result.alignments) == 3
    for alignment in result.alignments:
        assert len(alignment.left_only_atom_ids) == len(alignment.right_only_atom_ids) == 1
        assert len(alignment.left_boundaries) == len(alignment.right_boundaries) == 1
        assert len(alignment.atom_pairs) == 6
    assert compare_molecules('Cc1ccccc1', 'Clc1ccccc1') == result


@pytest.mark.parametrize('left,right', [
    ('C', 'N'), ('[13CH4]', 'C'), ('N', '[NH4+]'),
    ('c1ccccc1', 'C1CCCCC1'), ('C1CCCCC1', 'C1CCCC1'),
    ('c1ccccc1', 'c1ccc2ccccc2c1'),
])
def test_elements_isotopes_charge_aromaticity_and_complete_rings_are_strict(left, right):
    result = compare_molecules(left, right)
    assert result.status == 'no_verified_core'
    assert result.core_atom_count == 0 and not result.alignments
    assert not result.same_constitution


def test_positional_isomers_are_not_a_full_match():
    result = compare_molecules('Cc1ccccc1O', 'Cc1ccc(O)cc1')
    assert not result.same_constitution
    assert 6 <= result.core_atom_count < len(result.left.atoms)
    assert result.left_coverage < 1 and result.right_coverage < 1


@pytest.mark.parametrize('right,relationship', [
    ('C[C@H](O)Cl', 'same_specification'),
    ('C[C@@H](O)Cl', 'different_specified_stereochemistry'),
    ('CC(O)Cl', 'different_or_incomplete_specification'),
])
def test_stereo_relationship_distinguishes_missing_assignment(right, relationship):
    result = compare_molecules('C[C@H](O)Cl', right)
    assert result.same_constitution
    assert result.stereo_relationship == relationship
    assert result.core_atom_count == 4


def test_ez_isomers_and_unspecified_alkene_are_distinguished():
    assert compare_molecules('C/C=C/C', 'C/C=C\\C').stereo_relationship == 'different_specified_stereochemistry'
    result = compare_molecules('C/C=C/C', 'CC=CC')
    assert result.stereo_relationship == 'different_or_incomplete_specification'
    assert 'UNSPECIFIED_STEREOCHEMISTRY' in result.warnings


def test_selected_core_skips_search_and_does_not_auto_broaden(monkeypatch):
    monkeypatch.setattr(module.rdFMCS, 'FindMCS', lambda *args: pytest.fail('Explicit core needs no MCS'))
    result = compare_molecules('Cc1ccccc1', 'Clc1ccccc1', 'c1ccccc1')
    assert result.core_atom_count == 6 and result.core_method == 'supplied_core'
    missing = compare_molecules('Cc1ccccc1', 'Clc1ccccc1', 'c1ccncc1')
    assert missing.status == 'no_verified_core'
    assert not missing.alignments


def test_supplied_core_cannot_override_contradictory_charge_or_stereo():
    result = compare_molecules('CN', 'C[NH3+]', 'CN')
    assert not result.alignments
    result = compare_molecules('C[C@H](O)Cl', 'C[C@@H](O)Cl', 'C[C@H](O)Cl')
    assert not result.alignments


def test_mcs_timeout_retains_partial_evidence_and_warning(monkeypatch):
    monkeypatch.setattr(module.rdFMCS, 'FindMCS', lambda *args: SimpleNamespace(
        smartsString='[#6]-[#6]', canceled=True,
    ))
    result = compare_molecules('CCO', 'CCN')
    assert result.status == 'partial_timeout' and result.search_timed_out
    assert result.core_atom_count == 2
    assert 'MCS_TIMEOUT_CORE_NOT_PROVEN_MAXIMAL' in result.warnings


def test_embedding_limit_does_not_claim_exhaustive_counts(monkeypatch):
    policy = module.molecule_inspection_policy()
    policy['max_embeddings'] = 1
    monkeypatch.setattr(module, 'molecule_inspection_policy', lambda: policy)
    result = compare_molecules('Cc1ccccc1', 'Clc1ccccc1')
    assert result.embeddings_truncated and result.alignment_ambiguous
    assert 'EMBEDDING_LIMIT_ALIGNMENT_COUNTS_ARE_LOWER_BOUNDS' in result.warnings


def test_reembedding_rejects_induced_bond_conflicts():
    left, right = Chem.MolFromSmiles('C1CCC1'), Chem.MolFromSmiles('C1CCCC1')
    assert not module._compatible_alignment(left, right, (0, 1, 2, 3), (0, 1, 2, 3))


@pytest.mark.parametrize('seconds', [0, 6, True, 1.5])
def test_timeout_is_bounded(seconds):
    with pytest.raises(ValueError, match='timeout_seconds'):
        compare_molecules('CCO', 'CCN', timeout_seconds=seconds)
