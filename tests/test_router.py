from chemtools.router import detect_family, detect_family_from_reaction


def test_detect_family_ullmann():
    out = detect_family(["Brc1ccc(F)cc1", "Nc1ccccc1"])
    assert "family" in out and "confidence" in out


def test_detect_family_from_reaction_buchwald_pd():
    rsmi = "Brc1ccc(F)cc1.Nc1ccccc1>[Pd]>c1ccc(Nc2ccccc2)cc1"
    out = detect_family_from_reaction(rsmi, use_rxn_insight=False)
    assert out["family"] == "Buchwald_CN"
    assert out["rule"]["hits"].get("catalyst_pd") is True
    assert "Pd" in out.get("catalysts", {}).get("metals", [])


def test_detect_family_from_reaction_ullmann_cu():
    rsmi = "Brc1ccc(F)cc1.Nc1ccccc1>[Cu]I>c1ccc(Nc2ccccc2)cc1"
    out = detect_family_from_reaction(rsmi, use_rxn_insight=False)
    assert out["family"] == "Ullmann_CN"
    assert out["rule"]["hits"].get("catalyst_cu") is True
    assert "Cu" in out.get("catalysts", {}).get("metals", [])
