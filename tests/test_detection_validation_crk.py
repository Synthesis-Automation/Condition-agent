from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)


def test_validate_detection_with_crk_key_amide() -> None:
    crk_raw = (
        "|Bn-NH2|HeteroAr-CO2H -> Ar-Alkyl|Bn-NHCOR|HeteroAr-CONHR "
        "| bond_formed: C-N | bond_broken: C-O | spectators: Ar-F|HeteroAr-H|Pyrimidine"
    )
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Amide_formation"
    assert result["validation_method"] == "crk_pattern"
    assert result["corrected_from"] == "Unknown"


def test_validate_detection_with_crk_key_normalizes_double_dash_tokens() -> None:
    crk_raw = (
        "|Bn-NH2|HeteroAr-CO2H -> Ar-Alkyl|Bn--NHCOR|HeteroAr--CONHR "
        "| bond_formed: C-N | bond_broken: C-O"
    )
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Amide_formation"

