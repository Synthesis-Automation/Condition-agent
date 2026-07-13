from reactive_taxonomy import resolve_source_label, validate_source_label_mappings
from scripts.clean_reaction_label_csv import clean_rows


def test_reactive_source_label_crosswalk_matches_taxonomy_output() -> None:
    assert validate_source_label_mappings() == []


def test_qualified_labels_separate_handle_from_attachment_sterics() -> None:
    primary_amine = resolve_source_label("RNH2 a-branch")
    secondary_amine = resolve_source_label("R2NH a-branch")
    primary_alcohol = resolve_source_label("Alkyl-OH primary")

    assert primary_amine.base_label == "R-NH2"
    assert primary_amine.display_label == "R–NH₂ (α-C branched)"
    assert primary_amine.center_substitution_class == "primary"
    assert primary_amine.alpha_branched is True

    assert secondary_amine.base_label == "R1R2-NH"
    assert secondary_amine.display_label == "R¹R²–NH (≥1 α-C branched)"
    assert secondary_amine.center_substitution_class == "secondary"

    assert primary_alcohol.base_label == "R-OH"
    assert primary_alcohol.display_label == "RCH₂–OH"
    assert primary_alcohol.attachment_carbon_class == "primary"


def test_cleanup_maps_supported_labels_and_preserves_unsolved_labels() -> None:
    rows = [
        {"FG A": "ArBr", "FG B": "RNH2"},
        {"FG A": "RCO2H or M", "FG B": "RNH2 a-branch"},
        {"FG A": "ArBr", "FG B": "ArBr"},
        {"FG A": "", "FG B": ""},
        {"FG A": "Protecting Group", "FG B": "ArNH2"},
    ]

    cleaned, stats = clean_rows(rows)

    assert len(cleaned) == 2
    assert cleaned[0]["FG A Source"] == "ArBr"
    assert cleaned[0]["FG A"] == "Ar-Br"
    assert cleaned[0]["FG A Display"] == "Ar–Br"
    assert cleaned[0]["FG B"] == "R-NH2"
    assert cleaned[0]["FG B Center Class"] == "primary"

    assert cleaned[1]["FG A"] == "RCO2H or M"
    assert cleaned[1]["FG A Mapping Status"] == "unresolved"
    assert cleaned[1]["FG B"] == "R-NH2"
    assert cleaned[1]["FG B Alpha Branched"] == "true"
    assert cleaned[1]["FG B Mapping Status"] == "qualified"

    assert stats["matched_both_blank"] == 1
    assert stats["matched_identical"] == 2
    assert stats["matched_protecting_group"] == 1
    assert stats["removed_union"] == 3
