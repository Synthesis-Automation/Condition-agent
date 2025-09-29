def test_recommend_majority_core_and_constraints(monkeypatch):
    # Mock precedents to control votes
    precs = [
        {"reaction_id": "r1", "core": "Cu/L1", "base_uid": "B1", "solvent_uid": "S1", "T_C": 100, "time_h": 10},
        {"reaction_id": "r2", "core": "Cu/L1", "base_uid": "B2", "solvent_uid": "S1", "T_C": 110, "time_h": 12},
        {"reaction_id": "r3", "core": "Cu/L2", "base_uid": "B3", "solvent_uid": "S2", "T_C": 120, "time_h": 8},
    ]
    pack = {"prototype_id": "proto_x", "support": len(precs), "precedents": precs}

    import chemtools.precedent as prec
    import chemtools.recommend as rec

    monkeypatch.setattr(prec, "knn", lambda family, features, k=25, relax=None: pack)

    # Disallow S1 via constraints to force S2 pick even though S1 more frequent
    constraints = {"blacklist": ["S1"]}
    out = rec.recommend_from_reaction("Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1", k=5, relax={}, constraint_rules=constraints)

    assert out["recommendation"]["core"] == "Cu/L1"
    assert out["recommendation"]["solvent_uid"] == "S2"
    assert out["recommendation"]["T_C"] in (105.0, 110, 100)  # median or fallback

    role_pack = out.get("role_featurization")
    assert role_pack is None or isinstance(role_pack, dict)
    if isinstance(role_pack, dict):
        assignments = role_pack.get("assignments") or {}
        assert "electrophile" in assignments
        assert "nucleophile" in assignments
        elec_idx = assignments["electrophile"].get("index")
        reactants = role_pack.get("reactants") or []
        if isinstance(elec_idx, int) and 0 <= elec_idx < len(reactants):
            assert "smiles" in reactants[elec_idx]

    rule_features = out.get("rule_features")
    assert isinstance(rule_features, dict)
    assert rule_features.get("electrophile", {}).get("class")
    assert rule_features.get("nucleophile", {}).get("class")


def test_structured_recommendation_exposes_starting_materials(monkeypatch):
    pack = {"prototype_id": "proto_x", "support": 0, "precedents": []}

    import chemtools.precedent as prec
    import chemtools.recommend as rec

    monkeypatch.setattr(prec, "knn", lambda family, features, k=25, relax=None: pack)

    result = rec.recommend_conditions_structured(
        reaction="Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1",
        reaction_type=None,
        k=5,
        limit=3,
        relax={},
        constraints={},
    )

    starting = result.get("starting_materials")
    assert isinstance(starting, dict)
    assert "rule_features" in starting
    assert "role_featurization" in starting

