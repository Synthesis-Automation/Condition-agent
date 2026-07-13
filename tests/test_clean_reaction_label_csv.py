from scripts.clean_reaction_label_csv import (
    FG_LABEL_RULES,
    clean_rows,
    validate_label_rules,
)


def test_reactive_label_rules_match_taxonomy_output() -> None:
    validate_label_rules(FG_LABEL_RULES)


def test_cleanup_maps_supported_labels_and_preserves_unsolved_labels() -> None:
    rows = [
        {"FG A": "ArBr", "FG B": "RNH2"},
        {"FG A": "RCO2H or M", "FG B": "RNH2 a-branch"},
        {"FG A": "ArBr", "FG B": "ArBr"},
        {"FG A": "", "FG B": ""},
        {"FG A": "Protecting Group", "FG B": "ArNH2"},
    ]

    cleaned, stats = clean_rows(rows)

    assert cleaned == [
        {"FG A": "Ar-Br", "FG B": "R-NH2"},
        {"FG A": "RCO2H or M", "FG B": "RNH2 a-branch"},
    ]
    assert stats["matched_both_blank"] == 1
    assert stats["matched_identical"] == 2
    assert stats["matched_protecting_group"] == 1
    assert stats["removed_union"] == 3
