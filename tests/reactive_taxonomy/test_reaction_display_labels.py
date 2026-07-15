from reactive_taxonomy import featurize_reaction, load_reaction_label_rendering


def test_mapped_unknown_reaction_receives_observed_edit_label() -> None:
    result = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")

    assert result.named_family is None
    assert result.reaction_label == "C–N bond formation"
    assert result.reaction_label_status == "mapped_edit_summary"
    assert result.display_label is not None
    assert result.display_label.status == "observed_edits"
    assert result.display_label.detailed == "C(map 1)–N(map 2) bond formation"
    assert len(result.display_label.clauses) == 1
    assert result.display_label.clauses[0].evidence == "supplied_atom_mapping"


def test_multiple_edits_compose_and_collapse_repeated_generic_clauses() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")

    assert result.reaction_label == "C=C → C–C; 2 × H gain at C"
    assert result.display_label is not None
    assert len(result.display_label.clauses) == 3
    assert "H gain at C(map 1)" in result.display_label.detailed
    assert "H gain at C(map 2)" in result.display_label.detailed


def test_generic_label_is_invariant_to_reactant_component_order() -> None:
    forward = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")
    reversed_order = featurize_reaction("[NH2:2].[CH3:1]>>[CH3:1][NH2:2]")

    assert forward.reaction_label == reversed_order.reaction_label
    assert forward.display_label is not None
    assert reversed_order.display_label is not None
    assert forward.display_label.detailed == reversed_order.display_label.detailed


def test_ascii_rendering_changes_display_only() -> None:
    unicode_result = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")
    ascii_result = featurize_reaction(
        "[CH3:1].[NH2:2]>>[CH3:1][NH2:2]", label_style="ascii"
    )

    assert unicode_result.reaction_label == "C–N bond formation"
    assert ascii_result.reaction_label == "C-N bond formation"
    assert unicode_result.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert (
        unicode_result.reaction_signature.signature_id
        == ascii_result.reaction_signature.signature_id
    )


def test_exact_grammar_label_is_preserved_as_overlay() -> None:
    result = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")

    assert result.reaction_label == "Ar–Br + R–NH2 → Ar–NH–R"
    assert result.reaction_label_status == "exact_product"
    assert result.display_label is not None
    assert result.display_label.status == "exact_reconstruction"
    assert "C–Br bond cleavage" in result.display_label.detailed
    assert "C–N bond formation" in result.display_label.detailed
    assert "N–H loss" in result.display_label.detailed


def test_conflicting_evidence_is_visible_in_label_contract() -> None:
    reaction = (
        "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[NH2:8][CH3:9]>>"
        "[cH:2]1[c:3]([NH:8][CH3:9])[cH:4][cH:5][cH:6][cH:7]1"
    )
    result = featurize_reaction(reaction)

    assert result.reaction_label_status == "conflicting_edit_summary"
    assert result.reaction_label is not None
    assert result.reaction_label.startswith("Conflicting evidence:")
    assert result.display_label is not None
    assert result.display_label.status == "conflicting_evidence"
    assert result.display_label.confidence == 0.5
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.display_label.warnings


def test_display_label_serializes_as_nested_evidence() -> None:
    payload = featurize_reaction(
        "[CH3:1].[NH2:2]>>[CH3:1][NH2:2]"
    ).to_dict()

    assert payload["display_label"]["schema_version"] == "1.0"
    assert payload["display_label"]["clauses"][0]["edit_type"] == "formed"
    assert payload["display_label"]["clauses"][0]["atom_map_numbers"] == (1, 2)


def test_reaction_label_definition_is_versioned() -> None:
    rendering = load_reaction_label_rendering()

    assert rendering["schema_version"] == "1.0"
    assert rendering["label_schema_version"] == "1.0"
    assert set(rendering["clause_order"]) == {
        "formed",
        "broken",
        "order_changed",
        "hydrogen_change",
    }
