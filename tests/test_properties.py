from chemtools.properties import lookup


def test_lookup_base_k3po4():
    res = lookup("7778-53-2")
    assert res["found"] is True
    rec = res["record"]
    assert rec["uid"] == "7778-53-2"
    assert rec.get("role") == "BASE"
    assert abs(rec.get("pKa_DMSO", 0.0) - 30.0) < 1e-6


def test_lookup_solvent_water_kt():
    res = lookup("Water")
    assert res["found"] is True
    kt = res["record"].get("KT")
    assert isinstance(kt, dict)
    assert kt.get("alpha") and kt.get("beta") and kt.get("pi*")


def test_lookup_base_koh_by_token():
    res = lookup("KOH")
    assert res["found"] is True
    rec = res["record"]
    assert rec["uid"] == "1310-58-3"
    assert rec.get("pKa_water") == 15.7

