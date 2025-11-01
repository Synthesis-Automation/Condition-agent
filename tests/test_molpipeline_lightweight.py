import pytest

pytest.importorskip("molpipeline")

from chemtools.integrations import molpipeline as mp
from chemtools.featurizers import molpipeline as mp_feats
from chemtools.smiles import normalize


def _require_auto_to_mol(value: str):
    mol = mp.try_auto_to_mol(value)
    if mol is None:
        pytest.skip("MolPipeline AutoToMol cannot parse this input in the current environment.")
    return mol


def test_try_auto_to_mol_smiles_roundtrip():
    mol = mp.try_auto_to_mol("c1ccccc1O")
    assert mol is not None
    assert mol.GetNumAtoms() >= 6


def test_try_auto_to_mol_inchi_supports_normalize():
    inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
    _require_auto_to_mol(inchi)
    result = normalize(inchi)
    assert result["smiles_norm"]
    assert result["smiles_norm"].startswith("CC")


def test_morgan_fingerprint_dense_shape():
    vectors = mp_feats.morgan_fingerprint("c1ccccc1", return_sparse=False)
    assert vectors.shape == (1, 2048)
    assert vectors.dtype == float


def test_morgan_fingerprint_sparse_shape():
    matrix = mp_feats.morgan_fingerprint(["c1ccccc1"])
    assert matrix.shape[0] == 1
    assert matrix.shape[1] == 2048


def test_physchem_features_subset():
    descriptors = mp_feats.physchem_features(
        ["c1ccccc1"],
        descriptor_list=["TPSA", "MolLogP"],
    )
    assert descriptors.shape == (1, 2)
    assert descriptors.dtype == float
