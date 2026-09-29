"""Graph provenance, focused inspection and explicitly scoped stereo evidence."""

from dataclasses import asdict
import json

import pytest
from rdkit import Chem

from reactive_taxonomy import inspect_reactive_sites
from reactive_taxonomy.api import detect_reactive_site_hypotheses
from reactive_taxonomy.environments import build_site_environment
from reactive_taxonomy.molecular_motifs import detect_molecular_motifs
from reactive_taxonomy.molecule_inspection import molecule_inspection_policy
from reactive_taxonomy.search_fragments import suggest_search_fragments


def test_atom_ids_agree_with_fragments_and_parsed_input():
    smiles = '[OH:5][C@@H:4]([CH3:8])[CH2:3][Br:9]'
    result = inspect_reactive_sites(smiles)
    fragments = suggest_search_fragments(smiles)
    assert result.molecule.canonical_smiles == fragments.target_smiles
    assert [(a.atom_id, a.input_atom_index) for a in result.molecule.atoms] == [
        (a.atom_id, a.input_atom_index) for a in fragments.target_atoms
    ]
    original = Chem.MolFromSmiles(smiles)
    for atom in result.molecule.atoms:
        assert atom.element == original.GetAtomWithIdx(atom.input_atom_index).GetSymbol()
    assert 'ATOM_MAP_LABELS_IGNORED_NOT_REACTION_CORRESPONDENCE' in result.warnings
    json.dumps(result.to_dict())


def test_focused_descriptors_preserve_domain_result_and_other_sites():
    all_sites = inspect_reactive_sites('BrCCCCCCO')
    bromine = next(a.atom_id for a in all_sites.molecule.atoms if a.element == 'Br')
    result = inspect_reactive_sites('BrCCCCCCO', [bromine], radius=0)
    assert result.selected_atom_ids == result.neighborhood_atom_ids == (bromine,)
    assert result.sites and result.other_sites
    assert all(bromine in site.atom_indices for site in result.sites)
    assert all(bromine not in site.atom_indices for site in result.other_sites)
    molecule = Chem.MolFromSmiles(result.molecule.canonical_smiles)
    motifs = detect_molecular_motifs(molecule, 0)
    hypotheses = detect_reactive_site_hypotheses(molecule)
    assert set(s.hypothesis_id for s in result.sites + result.other_sites) == set(s.hypothesis_id for s in hypotheses)
    for site, environment in zip(result.sites, result.environments):
        assert asdict(environment) == asdict(build_site_environment(molecule, site, motifs))
        assert environment.reactivity_profile.definition_versions
    assert any('not rates' in limitation for limitation in result.limitations)


@pytest.mark.parametrize('smiles,kind,specified', [
    ('CC(O)Cl', 'Atom_Tetrahedral', 'Unspecified'),
    ('C[C@H](O)Cl', 'Atom_Tetrahedral', 'Specified'),
    ('CC=CC', 'Bond_Double', 'Unspecified'),
    ('C/C=C/C', 'Bond_Double', 'Specified'),
])
def test_potential_stereo_includes_unspecified_double_bonds(smiles, kind, specified):
    result = inspect_reactive_sites(smiles)
    assert any(s.stereo_type == kind and s.specification == specified for s in result.molecule.stereochemistry)
    assert ('UNSPECIFIED_STEREOCHEMISTRY' in result.warnings) == (specified == 'Unspecified')


def test_no_site_is_not_an_inertness_claim():
    result = inspect_reactive_sites('C', [0], radius=0)
    assert not result.sites
    assert 'NO_DETECTED_SITE_IN_REGION_NOT_PROOF_OF_INERTNESS' in result.warnings


@pytest.mark.parametrize('selected', [[], [True], [-1], [4], [0, 0], '0'])
def test_invalid_atom_selections_are_not_silently_repaired(selected):
    with pytest.raises(ValueError, match='selected_atom_ids'):
        inspect_reactive_sites('CCO', selected)


@pytest.mark.parametrize('smiles', ['', 'bad_smiles', 'CCO.[Na+]', '*CC', '[CH3]', 'C' * 201])
def test_invalid_or_unsupported_graphs_fail_explicitly(smiles):
    with pytest.raises(ValueError):
        inspect_reactive_sites(smiles)


def test_policy_validation_rejects_silent_relaxation(monkeypatch):
    import reactive_taxonomy.molecule_inspection as module
    policy = molecule_inspection_policy()
    assert policy['default_timeout_seconds'] == 2
    policy['complete_rings_only'] = False
    monkeypatch.setattr(module.json, 'loads', lambda _: policy)
    with pytest.raises(ValueError, match='semantics'):
        molecule_inspection_policy()
