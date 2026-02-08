import pandas as pd

from ml.simple_chanlam_model import (
    SIMPLE_CONDITION_COLUMNS,
    SIMPLE_NUMERIC_COLUMNS,
    SIMPLE_TEXT_COLUMNS,
    _feature_schema,
    _normalize_condition_space,
)


def test_normalize_condition_space_unique_rows() -> None:
    df = pd.DataFrame(
        {
            "catalyst": ["Cu(OAc)2", "Cu(OAc)2", None],
            "base": ["K2CO3", "K2CO3", "CsF"],
            "solvent": ["DCE", "DCE", "MeCN"],
        }
    )
    out = _normalize_condition_space(df, SIMPLE_CONDITION_COLUMNS)
    assert list(out.columns) == SIMPLE_CONDITION_COLUMNS
    assert len(out) == 2
    assert "NA" in out["catalyst"].values


def test_simple_feature_lists_non_empty() -> None:
    assert set(SIMPLE_CONDITION_COLUMNS) == {"catalyst", "base", "solvent"}
    assert len(SIMPLE_TEXT_COLUMNS) >= 2
    assert len(SIMPLE_NUMERIC_COLUMNS) >= 6


def test_feature_schema_adds_condition_props() -> None:
    cat_base, num_base, text_base = _feature_schema(False, "core_full")
    cat_props, num_props, text_props = _feature_schema(True, "core_full")
    assert text_base == text_props
    assert len(cat_props) > len(cat_base)
    assert len(num_props) > len(num_base)


def test_base_profile_uses_base_and_reactant_motif_spectator() -> None:
    cat, num, text = _feature_schema(False, "base_motif_spectator")
    assert cat == ["base"]
    assert set(num) == {"sulf_motif_count", "bor_motif_count"}
    assert "spectator_groups_tokens" in text
