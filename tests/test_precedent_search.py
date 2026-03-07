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
