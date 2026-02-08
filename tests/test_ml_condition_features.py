import pandas as pd

from ml.condition_features import add_condition_property_features


def test_add_condition_property_features_outputs_expected_columns() -> None:
    df = pd.DataFrame(
        {
            "catalyst": ["Cu(OAc)2"],
            "base": ["K2CO3"],
            "solvent": ["DCE"],
        }
    )
    out = add_condition_property_features(df)
    for col in (
        "catalyst_metal",
        "base_strength_band",
        "base_inorganic_band",
        "solvent_class",
        "is_cu_catalyst",
        "is_strong_base",
        "is_polar_solvent",
    ):
        assert col in out.columns
    assert out.loc[0, "catalyst_metal"] == "Cu"
