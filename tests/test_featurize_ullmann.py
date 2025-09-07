from chemtools.featurizers.ullmann import featurize
def test_featurize_ullmann():
    f = featurize("Brc1ccc(F)cc1","Nc1ccccc1")
    assert f["LG"] in ("Br","Cl","I","OTf","UNK")
    assert "bin" in f
