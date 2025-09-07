from chemtools.router import detect_family
def test_detect_family_ullmann():
    out = detect_family(["Brc1ccc(F)cc1", "Nc1ccccc1"])
    assert "family" in out and "confidence" in out
