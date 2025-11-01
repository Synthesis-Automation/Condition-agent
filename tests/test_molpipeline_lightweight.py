import pytest

pytest.importorskip("molpipeline")

from chemtools.integrations import molpipeline as mp
from chemtools.featurizers import molpipeline as mp_feats
from chemtools.featurizers import molecular as mol_feat
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


def test_molecular_featurize_includes_molpipeline_vectors():
    result = mol_feat.featurize(
        "Brc1ccccc1",
        "Nc1ccccc1",
    )
    mp_payload = result.get("molpipeline")
    assert isinstance(mp_payload, dict)
    settings = mp_payload.get("settings")
    assert isinstance(settings, dict)
    bits = settings.get("morgan_bits")
    assert isinstance(bits, int) and bits > 0
    descriptors = settings.get("physchem_descriptors")
    assert isinstance(descriptors, list) and len(descriptors) > 0

    elec_vec = mp_payload.get("electrophile")
    nuc_vec = mp_payload.get("nucleophile")
    assert isinstance(elec_vec, dict)
    assert isinstance(nuc_vec, dict)

    for vec in (elec_vec, nuc_vec):
        morgan = vec.get("morgan_fp")
        physchem = vec.get("physchem")
        physchem_map = vec.get("physchem_map")
        assert isinstance(morgan, list)
        assert len(morgan) == bits
        assert all(isinstance(x, float) for x in morgan)
        assert isinstance(physchem, list)
        assert len(physchem) == len(descriptors)
        assert all(isinstance(x, float) for x in physchem)
        assert isinstance(physchem_map, dict)
        assert set(physchem_map.keys()) == set(descriptors)
        assert all(isinstance(val, float) for val in physchem_map.values())
