import pytest


pytest.importorskip("drfp")


def _fake_rows():
    # Two Ullmann C–N candidates; R1 matches categorical bin, R2 does not.
    # R2's reaction_smiles equals the query used in the test to yield DRFP=1.0.
    return [
        {
            "reaction_id": "R1",
            "rxn_type": "Ullmann C–N",
            "yield_value": 50.0,
            "T_C": None,
            "time_h": None,
            "condition_core": "Cu/L1",
            "base_uid": "B1",
            "solvent_uid": "S1",
            "features": {"LG": "Br", "nuc_class": "aniline", "bin": "LG:Br|NUC:aniline"},
            "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>N(c1ccccc1)c1ccccc1",
        },
        {
            "reaction_id": "R2",
            "rxn_type": "Ullmann C–N",
            "yield_value": 50.0,
            "T_C": None,
            "time_h": None,
            "condition_core": "Cu/L2",
            "base_uid": "B2",
            "solvent_uid": "S2",
            "features": {"LG": "Cl", "nuc_class": "phenol", "bin": "LG:Cl|NUC:phenol"},
            "reaction_smiles": "Clc1ccccc1.Oc1ccccc1>>Oc1ccccc1c2ccccc2",
        },
    ]


def test_knn_drfp_integration_reranks(monkeypatch):
    import chemtools.precedent as prec

    # Use synthetic dataset and clear caches
    monkeypatch.setattr(prec, "_load", lambda: _fake_rows())
    prec._knn_cached.cache_clear()

    # Query equals R2 reaction_smiles to maximize DRFP similarity for R2.
    query_rsmi = _fake_rows()[1]["reaction_smiles"]
    features = {"LG": "Br", "nuc_class": "aniline", "bin": "LG:Br|NUC:aniline"}
    relax = {"use_drfp": True, "reaction_smiles": query_rsmi, "drfp_weight": 0.9, "strict_bin": False, "min_candidates": 2}
    out = prec.knn("Ullmann_CN", features, k=2, relax=relax)
    ids = [p.get("reaction_id") for p in (out.get("precedents") or [])]
    assert ids and ids[0] == "R2"

