import json

import pandas as pd

import scripts.A_convert_rule_to_hte_format as converter


def test_convert_rule_db_to_hte_csv(tmp_path, monkeypatch):
    def fake_featurize(smiles: str):
        mapping = {
            "Brc1ccccc1": {"motifs": [{"compound_id": "Ar-Br"}]},
            "Nc1ccccc1": {"motifs": [{"compound_id": "Ar-NH2"}]},
        }
        return mapping.get(smiles, {"motifs": []})

    monkeypatch.setattr(converter, "_cached_featurize", fake_featurize)

    rule_data = {
        "metadata": {"id": "cn_test", "name": "CN Test", "version": "v1"},
        "reaction": {
            "family": "Buchwald_Hartwig_C_N",
            "reference_reactions": ["Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"],
            "notes": "Test notes",
        },
        "applies_if": {"expr": "Ar-Br AND Ar-NH2"},
        "default_rule": {
            "id": "DEF1",
            "description": "Default rule",
            "conditions": {
                "pd_precatalyst": "Pd-XYZ",
                "ligand": "XPhos",
                "base": "K3PO4",
                "solvent": "dioxane",
            },
        },
        "base_rules": [
            {
                "id": "BR1",
                "name": "Base rule 1",
                "description": "Base rule",
                "reactant_features": {"expr": "Ar-NH2"},
                "conditions": {
                    "pd_source": "PdCl2(dtbpf)",
                    "ligand": "dtbpf",
                    "base": "Cs2CO3",
                    "solvent": "dioxane / water",
                    "additives": ["H2O"],
                },
            }
        ],
    }

    input_path = tmp_path / "rule.json"
    output_path = tmp_path / "rule.csv"
    input_path.write_text(json.dumps(rule_data), encoding="utf-8")

    converter.convert_rule_db_to_hte_csv(str(input_path), str(output_path), yield_value=95.0, z_score_value=2.5)

    df = pd.read_csv(output_path)
    assert len(df) == 2
    assert set(
        [
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
        ]
    ).issubset(df.columns)

    row = df.iloc[0]
    assert row["reaction_type"] == "Buchwald_Hartwig_C_N"
    assert row["reactant_1"] == "Ar-Br"
    assert row["reactant_2"] == "Ar-NH2"
