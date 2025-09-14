import numpy as np


def _fake_rows():
    # Two candidates in the same family. R1 matches categorical bin; R2 does not.
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
            "reaction_smiles": "LOW_RS",  # encodes to all zeros
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
            "reaction_smiles": "HIGH_RS",  # encodes to all ones
        },
    ]


def _patch_drfp(monkeypatch):
    import chemtools.reaction_similarity as rs

    def fake_available():
        return True

    def fake_encode(rsmi: str, n_bits: int = 4096, radius: int = 3):
        # Very small vectors are fine; Tanimoto stays consistent.
        if "HIGH" in rsmi:
            return np.ones(8, dtype="uint8")
        return np.zeros(8, dtype="uint8")

    monkeypatch.setattr(rs, "drfp_available", fake_available)
    monkeypatch.setattr(rs, "encode_drfp", fake_encode)
    monkeypatch.setattr(rs, "encode_drfp_cached", fake_encode)


def test_knn_drfp_reranks_neighbors(monkeypatch):
    import chemtools.precedent as prec

    # Use synthetic dataset
    monkeypatch.setattr(prec, "_load", lambda: _fake_rows())
    # Ensure cache is fresh
    prec._knn_cached.cache_clear()

    # Force DRFP path
    _patch_drfp(monkeypatch)

    features = {"LG": "Br", "nuc_class": "aniline", "bin": "LG:Br|NUC:aniline"}
    relax = {
        "use_drfp": True,
        "reaction_smiles": "HIGH_QUERY",  # encodes to ones
        "drfp_weight": 0.9,  # DRFP dominates
        "strict_bin": False,
        "min_candidates": 2,
    }
    out = prec.knn("Ullmann_CN", features, k=2, relax=relax)
    ids = [p.get("reaction_id") for p in (out.get("precedents") or [])]
    # DRFP should push R2 (HIGH sim) above R1 (LOW sim), despite worse categorical features
    assert ids and ids[0] == "R2"


def test_knn_no_drfp_prefers_categorical(monkeypatch):
    import chemtools.precedent as prec

    monkeypatch.setattr(prec, "_load", lambda: _fake_rows())
    prec._knn_cached.cache_clear()

    features = {"LG": "Br", "nuc_class": "aniline", "bin": "LG:Br|NUC:aniline"}
    relax = {"use_drfp": False, "strict_bin": False, "min_candidates": 2}
    out = prec.knn("Ullmann_CN", features, k=2, relax=relax)
    ids = [p.get("reaction_id") for p in (out.get("precedents") or [])]
    # Without DRFP, exact bin match (R1) should come first
    assert ids and ids[0] == "R1"
