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


def test_validate_detection_with_crk_key_avoids_snar_for_electron_rich_aryl_halide() -> None:
    crk_raw = (
        "|Ar-Br|R2CH-NH2 -> Ar-NHR "
        "| bond_formed: C(ar)-N | bond_broken: Br-C(ar) "
        "| spectators: Ar-OR"
    )
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "C_N_Coupling"


def test_validate_detection_with_crk_key_keeps_snar_when_heteroaryl_activation_present() -> None:
    crk_raw = (
        "|Ar-Cl|HeteroAr-NH2 -> Ar-NHR "
        "| bond_formed: C(ar)-N | bond_broken: C(ar)-Cl "
        "| spectators: HeteroAr-H|Pyrimidine"
    )
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "SNAr_CN"
