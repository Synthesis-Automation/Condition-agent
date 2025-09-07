from chemtools.smiles import normalize
def test_normalize_basic():
    res = normalize("Brc1ccc(F)cc1.O")
    assert res["largest_smiles"]
    assert res["smiles_norm"]
