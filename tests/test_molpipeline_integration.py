import numpy as np
import pytest

from chemtools.integrations import molpipeline as mp

molpipeline = pytest.importorskip("molpipeline")


def _build_role_aggregator() -> mp.MolPipelineRoleAggregator:
    return mp.MolPipelineRoleAggregator(
        {
            "LIGAND": mp.build_morgan_pipeline(n_bits=8, radius=1, n_jobs=1),
            "BASE": mp.build_physchem_pipeline(["HeavyAtomMolWt"], n_jobs=1),
            "SOLVENT": mp.build_physchem_pipeline(["TPSA"], n_jobs=1),
        }
    )


def test_collect_role_smiles_basic():
    reaction = {
        "reagents": [
            {"role": "ligand", "smiles": "c1ccccc1"},
            {"role": "base", "smiles": "CCO"},
        ]
    }
    role_smiles = mp.collect_role_smiles(reaction)
    assert role_smiles["LIGAND"] == ["c1ccccc1"]
    assert role_smiles["BASE"] == ["CCO"]


def test_role_aggregator_features_and_concatenation():
    aggregator = _build_role_aggregator()
    reaction = {
        "reagents": [
            {"role": "ligand", "smiles": "c1ccccc1"},
            {"role": "base", "smiles": "CCO"},
        ]
    }
    features = aggregator.featurize_roles(reaction=reaction)
    assert set(features) == {"LIGAND", "BASE", "SOLVENT"}
    assert features["LIGAND"].shape == (8,)
    assert features["BASE"].shape == (1,)
    assert np.allclose(features["SOLVENT"], 0.0)

    concatenated = aggregator.concatenate(reaction=reaction)
    assert concatenated.shape == (10,)
    # Last element corresponds to SOLVENT TPSA descriptor (missing -> zero fill)
    assert np.isclose(concatenated[-1], 0.0)


def test_role_aggregator_role_smiles_override():
    aggregator = _build_role_aggregator()
    role_smiles = {
        "LIGAND": ["c1ccccc1"],
        "BASE": ["CCO"],
    }
    features = aggregator.featurize_roles(role_smiles=role_smiles)
    assert features["SOLVENT"].shape == (1,)  # zeros fallback
    concat = aggregator.concatenate(role_smiles=role_smiles)
    assert concat.shape == (10,)
