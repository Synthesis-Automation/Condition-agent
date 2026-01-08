import json

import pandas as pd

import scripts.A_convert_to_hte_format as converter


def test_process_reaction_dataset_emits_canonical_columns(tmp_path, monkeypatch):
    def fake_featurize(smiles: str):
        mapping = {
            "Brc1ccc(C(F)(F)F)cc1": {
                "motifs": [
                    {"compound_id": "Ar-Br", "alt_compound_ids": ["Ar-CF3"]}
                ]
            },
            "Nc1ccccc1": {
                "motifs": [
                    {"compound_id": "Ar-NH2", "alt_compound_ids": ["Any-NH2"]}
                ]
            },
            "FC(F)(F)c1ccc(Nc2ccccc2)cc1": {
                "motifs": [
                    {"compound_id": "Ar-NH", "alt_compound_ids": []}
                ]
            },
        }
        return mapping.get(smiles, {"motifs": []})

    monkeypatch.setattr(converter, "cached_featurize", fake_featurize)

    record = {
        "smiles": "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1",
        "yield": 48.0,
        "reagents": [
            {"name": "K2CO3", "role": "BASE"},
            {"name": "KOPiv", "role": "ADDITIVE"},
        ],
        "solvents": [{"name": "DMAc"}, {"name": "water"}],
        "catalytic_system": "PdCl2, XantPhos",
    }

    input_path = tmp_path / "input.jsonl"
    output_path = tmp_path / "output.csv"
    input_path.write_text(json.dumps(record) + "\n", encoding="utf-8")

    converter.process_reaction_dataset(str(input_path), str(output_path))

    df = pd.read_csv(output_path)
    assert len(df) == 1

    expected_cols = {
        "reaction_type",
        "yield",
        "z_score",
        "reactant_1",
        "reactant_2",
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "additive",
        "reaction_smiles",
        "reaction_key",
        "reacted_motifs",
        "formed_motifs",
        "spectator_motifs",
    }
    assert expected_cols.issubset(set(df.columns))

    row = df.iloc[0]
    assert row["reaction_smiles"] == record["smiles"]
    assert row["reactant_1"] == "Ar-Br"
    assert row["reactant_2"] == "Ar-NH2"
    assert row["catalyst"] == "PdCl2"
    assert row["ligand"] == "XantPhos"
