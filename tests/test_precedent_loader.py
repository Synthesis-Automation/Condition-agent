from chemtools.precedent import loader


def test_make_row_from_csv_prefers_detected_reaction_type() -> None:
    rec = {
        "reaction_id": "C_N_Coupling",
        "detected_reaction_type": "Azide_coupling",
        "reaction_smiles": "Brc1ccccc1.NCC>>NCCc1ccccc1",
        "catalyst": "Pd(OAc)2",
        "ligand": "SPhos",
        "base": "K3PO4",
        "solvent": "THF",
        "yield": "80",
    }

    row = loader._make_row_from_csv(
        rec,
        row_index=0,
        file_family="C_N_Coupling",
        source_group="literature",
        fast=True,
    )

    assert row is not None
    assert row["rxn_type"] == "Azide_coupling"
    assert row["dataset_reaction_id"] == "Azide_coupling"
    assert row["dataset_reaction_id_raw"] == "C_N_Coupling"


def test_load_literature_cached_filters_on_row_family(monkeypatch) -> None:
    def fake_iter_files():
        return ["fake/C_N_Coupling_canonical.csv"]

    def fake_source_group(path: str) -> str:  # noqa: ARG001
        return "literature"

    def fake_read_records(path: str):  # noqa: ARG001
        return [
            {
                "reaction_id": "C_N_Coupling",
                "detected_reaction_type": "Azide_coupling",
                "reaction_smiles": "Brc1ccccc1.NCC>>NCCc1ccccc1",
                "catalyst": "Pd(OAc)2",
            },
            {
                "reaction_id": "C_N_Coupling",
                "detected_reaction_type": "C_N_Coupling",
                "reaction_smiles": "Brc1ccccc1.NCC>>NCCc1ccccc1",
                "catalyst": "Pd(OAc)2",
            },
        ]

    monkeypatch.setattr(loader, "_iter_literature_files", fake_iter_files)
    monkeypatch.setattr(loader, "_infer_source_group_from_path", fake_source_group)
    monkeypatch.setattr(loader, "_read_csv_records", fake_read_records)
    monkeypatch.setattr(loader, "_load_precedent_disk_cache", lambda family_key: None)
    monkeypatch.setattr(loader, "_save_precedent_disk_cache", lambda family_key, rows: None)

    loader._load_literature_cached.cache_clear()
    rows = loader._load_literature_cached(("Azide_coupling",))

    assert len(rows) == 1
    assert rows[0]["rxn_type"] == "Azide_coupling"
