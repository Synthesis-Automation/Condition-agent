from reactive_taxonomy import resolve_source_label, validate_source_label_mappings
from scripts.clean_reaction_label_csv import OUTPUT_FIELDNAMES, clean_rows


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


def test_alkene_source_label_maps_to_generic_pi_handle() -> None:
    alkene = resolve_source_label("alkene")

    assert alkene.base_label == "H2C=CH2"
    assert alkene.canonical_signature == "PI|Alkene"
    assert alkene.mapping_status == "canonical"


def test_acid_or_carboxylate_source_label_maps_to_latent_acyl_donor() -> None:
    acid = resolve_source_label("RCO2H or M")

    assert acid.base_label == "R-C(O)OH"
    assert acid.display_label == "R–C(O)OH"
    assert acid.canonical_signature == "EC|Acyl|Alkyl|OH|latent"
    assert acid.mapping_status == "qualified"
    assert acid.qualifier_scope == "acid_or_carboxylate_collapsed_to_acid"


def test_acidic_alkyl_h_maps_to_qualified_activated_carbon_family() -> None:
    acidic_carbon = resolve_source_label("Alkyl-H acidic")

    assert acidic_carbon.base_label == "activated C-H"
    assert acidic_carbon.display_label == "Activated alkyl C–H"
    assert acidic_carbon.canonical_signature == "XH|Csp3"
    assert acidic_carbon.signature_scope == "family"
    assert acidic_carbon.mapping_status == "qualified"


def test_cleanup_maps_supported_labels_and_filters_invalid_rows() -> None:
    rows = [
        {"FG A": "ArBr", "FG B": "RNH2"},
        {"FG A": "RCO2H or M", "FG B": "RNH2 a-branch"},
        {"FG A": "ArBr", "FG B": "ArBr"},
        {"FG A": "", "FG B": ""},
        {"FG A": "Protecting Group", "FG B": "ArNH2"},
    ]

    cleaned, stats = clean_rows(rows)

    assert len(cleaned) == 2
    assert set(cleaned[0]) == set(OUTPUT_FIELDNAMES)
    assert cleaned[0]["reactive_site_1_display_label"] == "Ar–Br"
    assert cleaned[0]["reactive_site_2_display_label"] == "R–NH₂"
    assert cleaned[0]["reactive_site_1_signature"] == "LG|Ar|Br"
    assert cleaned[0]["reactive_site_2_signature"] == "XH|N|H2|Alkyl"
    assert cleaned[0]["reactive_site_2_center_class"] == "primary"

    assert cleaned[1]["reactive_site_1_signature"] == "EC|Acyl|Alkyl|OH|latent"
    assert cleaned[1]["reactive_site_2_alpha_branched"] == "true"

    assert stats["matched_both_blank"] == 1
    assert stats["matched_identical"] == 2
    assert stats["matched_protecting_group"] == 1
    assert stats["removed_union"] == 3


def test_output_schema_keeps_display_labels_but_excludes_dead_metadata() -> None:
    excluded = {
        f"reactive_site_{index}_{suffix}"
        for index in (1, 2)
        for suffix in (
            "source_label",
            "normalized_label",
            "qualifier_scope",
            "mapping_status",
        )
    }

    assert len(OUTPUT_FIELDNAMES) == 22
    assert "reactive_site_1_display_label" in OUTPUT_FIELDNAMES
    assert "reactive_site_2_display_label" in OUTPUT_FIELDNAMES
    assert excluded.isdisjoint(OUTPUT_FIELDNAMES)
