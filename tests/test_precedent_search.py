from chemtools.precedent import search


def test_knn_precedent_limit_respects_relax_override(monkeypatch) -> None:
    rows = []
    for i in range(40):
        rows.append(
            {
                "reaction_id": f"r:{i}",
                "dataset_reaction_id": "C_N_Coupling",
                "rxn_type": "C_N_Coupling",
                "reaction_smiles": f"A{i}.B{i}>>P{i}",
                "features": {"bin": "BIN_TEST", "LG": "Br", "nuc_class": "amine"},
                "yield_value": 80.0 - (i * 0.1),
                "condition_core": "Pd/SPhos",
                "base_uid": "K3PO4",
                "solvent_uid": "THF",
                "reagents": [],
                "solvents": [],
                "source_file": "C_N_Coupling",
                "source_group": "literature",
                "reference": "",
                "conditions": {},
                "catalyst": {},
                "full_system": None,
                "catalytic_system": [],
                "T_C": None,
                "time_h": None,
            }
        )

    monkeypatch.setattr(search, "_load", lambda: rows)
    monkeypatch.setattr(search, "_load_selective", lambda families=None: rows)
    search._knn_cached.cache_clear()

    out = search.knn(
        family="C_N_Coupling",
        features={"bin": "BIN_TEST", "LG": "Br", "nuc_class": "amine"},
        k=15,
        relax={
            "filter_by_reagent_database": False,
            "selective_loading": False,
            "precedent_limit": 25,
        },
    )

    precedents = list(out.get("precedents", []) or [])
    assert len(precedents) == 25


def test_knn_limits_on_demand_drfp_reranking(monkeypatch) -> None:
    rows = []
    for i in range(30):
        rows.append(
            {
                "reaction_id": f"r:{i}",
                "dataset_reaction_id": "Suzuki_miyaura",
                "rxn_type": "Suzuki_miyaura",
                "reaction_smiles": f"A{i}.B{i}>>P{i}",
                "features": {"bin": "BIN_TEST", "LG": "Br", "nuc_class": "amine"},
                "yield_value": 90.0 - i,
                "condition_core": "Pd/SPhos",
                "base_uid": "K3PO4",
                "solvent_uid": "THF",
                "reagents": [],
                "solvents": [],
                "source_file": "Suzuki_miyaura",
                "source_group": "literature",
                "reference": "",
                "conditions": {},
                "catalyst": {},
                "full_system": None,
                "catalytic_system": [],
                "T_C": None,
                "time_h": None,
                "precomputed": {},
            }
        )

    calls = {"count": 0}

    class FakeRS:
        @staticmethod
        def drfp_available():
            return True

        @staticmethod
        def encode_drfp_cached(rsmi, n_bits=4096, radius=3):
            calls["count"] += 1
            return [1, 0, 1]

        @staticmethod
        def tanimoto(a, b):
            return 0.8

    monkeypatch.setattr(search, "_load", lambda: rows)
    monkeypatch.setattr(search, "_load_selective", lambda families=None: rows)
    monkeypatch.setattr(search, "_family_text", lambda family: family)
    monkeypatch.setattr(search, "rs", FakeRS)
    monkeypatch.setattr(search, "_drfp_storage_available", False)
    search._knn_cached.cache_clear()

    out = search.knn(
        family="Suzuki_miyaura",
        features={"bin": "BIN_TEST", "LG": "Br", "nuc_class": "amine"},
        k=5,
        relax={
            "filter_by_reagent_database": False,
            "selective_loading": True,
            "use_drfp": True,
            "reaction_smiles": "A.B>>P",
            "drfp_rerank_limit": 7,
        },
    )

    # One query fingerprint + seven on-demand candidate fingerprints.
    assert calls["count"] == 8
    assert out["drfp_load_strategy"]["computed"] == 7
    assert out["drfp_rerank_limit"] == 7
    assert out["drfp_rerank_candidates"] == 30
