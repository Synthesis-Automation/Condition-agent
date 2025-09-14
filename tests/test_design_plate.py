def test_design_plate_24_across_cores(monkeypatch):
    # Build a fake precedent pack with 3 cores and some bases/solvents
    precs = []
    for i in range(10):
        precs.append({"reaction_id": f"rA{i}", "core": "Cu/L1", "base_uid": "B1", "solvent_uid": "S1", "T_C": 100, "time_h": 10})
    for i in range(5):
        precs.append({"reaction_id": f"rB{i}", "core": "Cu/L2", "base_uid": "B2", "solvent_uid": "S2", "T_C": 110, "time_h": 12})
    for i in range(3):
        precs.append({"reaction_id": f"rC{i}", "core": "Cu/L3", "base_uid": "B3", "solvent_uid": "S3", "T_C": 90, "time_h": 8})
    pack = {"prototype_id": "proto_x", "support": len(precs), "precedents": precs}

    import chemtools.precedent as prec
    import chemtools.recommend as rec

    monkeypatch.setattr(prec, "knn", lambda family, features, k=50, relax=None: pack)

    out = rec.design_plate_from_reaction("Brc1ccc(F)cc1.Nc1ccccc1>>Nc1ccc(F)cc1", plate_size=24, relax={}, constraint_rules={})
    rows = out.get("rows") or []
    assert len(rows) == 24
    # Check first three wells diversify across the three cores (in frequency order)
    first_three = [rows[i]["core"] for i in range(3)]
    assert set(first_three) == {"Cu/L1", "Cu/L2", "Cu/L3"}
    # CSV header present
    csv_text = out.get("csv") or ""
    assert csv_text.startswith("well_id,core,base_uid,solvent_uid,additive_uids,T_C,time_h\n")

