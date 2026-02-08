import pandas as pd

from ml.tuning import (
    evaluate_setting,
    filter_rows_by_min_stage2,
    grid_search_defaults,
    score_composite_objective,
    select_best_by_composite,
)


def _toy_predictions() -> pd.DataFrame:
    # Two folds, two reactions/fold, three conditions each.
    rows = []
    for fold in (1, 2):
        for rxn, grp in (("R1", "G1"), ("R2", "G2")):
            rows.extend(
                [
                    {
                        "fold": fold,
                        "group": grp,
                        "reaction_smiles": rxn,
                        "condition_yield": 80.0,
                        "stage1_score": 0.9,
                        "stage2_score": 0.6,
                    },
                    {
                        "fold": fold,
                        "group": grp,
                        "reaction_smiles": rxn,
                        "condition_yield": 40.0,
                        "stage1_score": 0.5,
                        "stage2_score": 0.7,
                    },
                    {
                        "fold": fold,
                        "group": grp,
                        "reaction_smiles": rxn,
                        "condition_yield": 10.0,
                        "stage1_score": 0.1,
                        "stage2_score": 0.2,
                    },
                ]
            )
    return pd.DataFrame(rows)


def test_evaluate_setting_returns_metrics() -> None:
    df = _toy_predictions()
    result = evaluate_setting(df, shortlist_k=2, w_stage1=0.7, w_stage2=0.3)
    assert result.apyr_mean >= 0.0
    assert result.spearman_mean >= -1.0
    assert result.top1_hit_ge_70 >= 0.0


def test_grid_search_defaults_finds_candidate() -> None:
    df = _toy_predictions()
    best, rows = grid_search_defaults(df, shortlist_grid=[1, 2, 3], w_stage1_grid=[0.0, 0.5, 1.0])
    assert len(rows) == 9
    assert best.shortlist_k in {1, 2, 3}
    assert 0.0 <= best.w_stage1 <= 1.0


def test_composite_selection_prefers_stage2_when_weighted() -> None:
    rows = [
        {
            "shortlist_k": 30,
            "w_stage1": 1.0,
            "w_stage2": 0.0,
            "apyr_mean": 85.0,
            "spearman_mean": 0.40,
            "top1_hit_ge_70": 0.50,
            "top3_hit_ge_70": 0.55,
        },
        {
            "shortlist_k": 20,
            "w_stage1": 0.4,
            "w_stage2": 0.6,
            "apyr_mean": 84.0,
            "spearman_mean": 0.55,
            "top1_hit_ge_70": 0.50,
            "top3_hit_ge_70": 0.55,
        },
    ]
    scored = score_composite_objective(
        rows,
        w_apyr=0.40,
        w_top1=0.10,
        w_spearman=0.35,
        w_top3=0.05,
        w_stage2_pref=0.10,
    )
    best = select_best_by_composite(scored)
    assert best.w_stage2 > 0.0


def test_filter_rows_by_min_stage2_enforces_threshold() -> None:
    rows = [
        {"w_stage2": 0.0, "shortlist_k": 30},
        {"w_stage2": 0.2, "shortlist_k": 20},
    ]
    filtered = filter_rows_by_min_stage2(rows, min_stage2_weight=0.1)
    assert len(filtered) == 1
    assert filtered[0]["shortlist_k"] == 20
