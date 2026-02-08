import argparse
import json

import pandas as pd

import app.HTE_analytics_cli as cli


class _FakeRecommendation:
    catalyst = "Pd"
    ligand = "PPh3"
    base = "K2CO3"
    solvent = "DMF"
    secondary_solvent = ""
    additive = ""
    coupling_reagent = ""


class _FakeResult:
    def __init__(self) -> None:
        self.recommendations = [_FakeRecommendation()]
        self.database_coverage = 12.5
        self.total_matching_experiments = 7


class _FakeRecommender:
    def __init__(self, hte_db_path: str) -> None:
        self.hte_db_path = hte_db_path

    def recommend(self, **kwargs):
        return _FakeResult()


def test_cmd_backtest_reports_hits(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(cli, "HTERecommender", _FakeRecommender)

    df = pd.DataFrame(
        {
            "reaction_smiles": [
                "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
                "Brc1ccc(F)cc1.Nc1ccccc1>>Fc1ccc(Nc2ccccc2)cc1",
                "Brc1ccc(Cl)cc1.Nc1ccccc1>>Clc1ccc(Nc2ccccc2)cc1",
                "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1",
            ],
            "Catalyst": ["Pd", "Pd", "Pd", "Pd"],
            "Ligand": ["PPh3", "PPh3", "PPh3", "PPh3"],
            "Base": ["K2CO3", "K2CO3", "K2CO3", "K2CO3"],
            "Solvent": ["DMF", "DMF", "DMF", "DMF"],
            "Reaction_Type_Standardized": ["C_N_Coupling"] * 4,
            "Source_Group": ["literature"] * 4,
        }
    )
    input_csv = tmp_path / "backtest_input.csv"
    df.to_csv(input_csv, index=False)

    report_json = tmp_path / "report.json"
    per_row_csv = tmp_path / "per_row.csv"
    args = argparse.Namespace(
        input=str(input_csv),
        reaction=None,
        catalyst=None,
        source_group=None,
        test_fraction=0.4,
        test_size=2,
        seed=11,
        top_k=3,
        hit_ks="1,3",
        min_experiments=1,
        reaction_key_only=False,
        use_spectator_groups=True,
        allow_unseen_conditions=False,
        train_output=None,
        output=str(report_json),
        per_row_output=str(per_row_csv),
    )

    cli.cmd_backtest(args)

    report = json.loads(report_json.read_text(encoding="utf-8"))
    assert report["evaluated_rows"] == 2
    assert report["metrics"]["count"] == 2
    assert report["metrics"]["hit@1"] == 1.0
    assert report["metrics"]["hit@3"] == 1.0
    assert report["metrics"]["query_coverage"] == 1.0

    per_row = pd.read_csv(per_row_csv)
    assert len(per_row) == 2
    assert set(per_row["rank"]) == {1}


def test_cmd_backtest_skips_unseen_conditions_by_default(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(cli, "HTERecommender", _FakeRecommender)

    df = pd.DataFrame(
        {
            "reaction_smiles": [
                "A.B>>P1",
                "A.B>>P2",
                "A.B>>P3",
            ],
            "Catalyst": ["Pd", "Ni", "Cu"],
            "Ligand": ["L1", "L2", "L3"],
            "Base": ["B1", "B2", "B3"],
            "Solvent": ["S1", "S2", "S3"],
        }
    )
    input_csv = tmp_path / "backtest_unique.csv"
    df.to_csv(input_csv, index=False)

    report_json = tmp_path / "report_unique.json"
    args = argparse.Namespace(
        input=str(input_csv),
        reaction=None,
        catalyst=None,
        source_group=None,
        test_fraction=0.34,
        test_size=1,
        seed=3,
        top_k=5,
        hit_ks="1,3",
        min_experiments=1,
        reaction_key_only=False,
        use_spectator_groups=True,
        allow_unseen_conditions=False,
        train_output=None,
        output=str(report_json),
        per_row_output=None,
    )

    cli.cmd_backtest(args)

    report = json.loads(report_json.read_text(encoding="utf-8"))
    assert report["test_rows"] == 1
    assert report["evaluated_rows"] == 0
    assert report["skipped_unseen_condition"] == 1
    assert report["metrics"]["count"] == 0
    assert report["metrics"]["hit@1"] == 0.0
