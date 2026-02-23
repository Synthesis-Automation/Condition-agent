from __future__ import annotations

from pathlib import Path

import pandas as pd

from chemtools.recommend.analytics import HTEAnalytics, _collect_hte_files, _infer_source_group


def test_collect_hte_files_recursively_includes_protocols_and_nested_dirs(tmp_path: Path) -> None:
    (tmp_path / "experiments").mkdir()
    (tmp_path / "literature").mkdir()
    (tmp_path / "protocols").mkdir()
    (tmp_path / "nested" / "custom").mkdir(parents=True)

    files = [
        tmp_path / "experiments" / "HTE_canonical.csv",
        tmp_path / "literature" / "test.jsonl",
        tmp_path / "protocols" / "protocols_all_v2_hte.csv",
        tmp_path / "nested" / "custom" / "extra.csv",
    ]
    for path in files:
        path.write_text("col\n1\n", encoding="utf-8")

    found = _collect_hte_files(tmp_path)
    found_set = {p.resolve() for p in found}

    assert {p.resolve() for p in files}.issubset(found_set)


def test_infer_source_group_recognizes_protocols() -> None:
    path = Path("data/HTE_db/protocols/protocols_all_v2_hte.csv")
    assert _infer_source_group(path) == "protocols"


def test_hte_analytics_normalizes_detected_reaction_type_from_literature_exports(tmp_path: Path) -> None:
    (tmp_path / "experiments").mkdir()
    (tmp_path / "literature").mkdir()

    pd.DataFrame(
        [
            {
                "reaction_type": "C_S_Coupling",
                "reactant_1": "Ar-Br",
                "reactant_2": "Ar-SH",
                "yield": 55.0,
                "catalyst": "Pd",
            }
        ]
    ).to_csv(tmp_path / "experiments" / "HTE_canonical.csv", index=False)

    pd.DataFrame(
        [
            {
                "detected_reaction_type": "C_S_Coupling",
                "reactant_1": "Ar-Br",
                "reactant_2": "Ar-SH",
                "yield": 60.0,
                "catalyst": "Pd",
            },
            {
                "detected_reaction_type": "C_S_Coupling",
                "reactant_1": "Ar-I",
                "reactant_2": "Ar-SH",
                "yield": 20.0,
                "catalyst": "Cu",
            },
        ]
    ).to_csv(tmp_path / "literature" / "C_S_Coupling_canonical.csv", index=False)

    analytics = HTEAnalytics(str(tmp_path))
    summary = analytics.get_reaction_type_summary()
    row = summary.loc[summary["Reaction_Type"] == "C_S_Coupling"].iloc[0]

    assert int(row["Num_Experiments"]) == 3
    assert int(row["Num_Reactant_Pairs"]) == 2
