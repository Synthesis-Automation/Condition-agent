import pandas as pd

from ml.tune_and_holdout_chanlam import _derive_external_holdout_groups


def test_derive_external_holdout_groups_partial_coverage() -> None:
    df = pd.DataFrame(
        {
            "sulfonamide_smiles": ["A", "A", "B", "B", "C"],
            "catalyst": ["c1", "c2", "c1", "c2", "c1"],
            "ligand": ["NA"] * 5,
            "base": ["b1", "b2", "b1", "b2", "b1"],
            "acid": ["NA"] * 5,
            "oxidant": ["NA"] * 5,
            "reductant": ["NA"] * 5,
            "additive": ["NA"] * 5,
            "condensation_agent": ["NA"] * 5,
            "other_reagent": ["NA"] * 5,
            "solvent": ["s1", "s1", "s1", "s1", "s1"],
        }
    )
    # A and B have 2 unique conditions each (max coverage), C has 1 (holdout)
    holdout = _derive_external_holdout_groups(df)
    assert holdout == ["C"]
