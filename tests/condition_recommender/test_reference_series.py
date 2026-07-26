from condition_recommender.conversion.reference_series import (
    reference_condition_series_id,
)


def test_reference_series_groups_variants_with_one_core_and_stage_count() -> None:
    first = reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
    )
    repeated = reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
    )

    assert first.startswith("RCS1:")
    assert first == repeated


def test_reference_series_separates_recipe_stage_and_reference() -> None:
    base = reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
    )

    assert base != reference_condition_series_id(
        reference_id="REF1:other",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
    )
    assert base != reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:other",
        chemistry_key="L2:chemistry",
        stages="1",
    )
    assert base != reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="2",
    )
    assert base != reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
        temperature_c=80.0,
    )
    assert base != reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:chemistry",
        stages="1",
        steps="2",
    )
    assert base != reference_condition_series_id(
        reference_id="REF1:paper",
        recipe_core_id="RCORE1:recipe",
        chemistry_key="L2:other",
        stages="1",
    )


def test_reference_series_requires_reference_and_recipe_core() -> None:
    assert (
        reference_condition_series_id(
            reference_id="",
            recipe_core_id="RCORE1:recipe",
            chemistry_key="L2:chemistry",
            stages="1",
        )
        == ""
    )
    assert (
        reference_condition_series_id(
            reference_id="REF1:paper",
            recipe_core_id="RCORE1:recipe",
            chemistry_key="",
            stages="1",
        )
        == ""
    )
