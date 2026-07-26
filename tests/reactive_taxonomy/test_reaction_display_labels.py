from reactive_taxonomy import (
    featurize_reaction,
    load_reaction_label_patterns,
    load_reaction_label_rendering,
)


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

    assert result.reaction_label == "H2C=CH2 → H3C–CH3"
    assert result.reaction_label_status == "mapped_generic_pattern"
    assert result.display_label is not None
    assert result.display_label.status == "generic_pattern"
    assert result.display_label.pattern_id == "hydrogenation"
    assert result.display_label.structural_label == "C=C → C–C; 2 × H gain at C"
    assert result.display_label.transformation_label == "C=C hydrogenation"
    assert result.display_label.contextual_label == "H2C=CH2 → H3C–CH3"
    assert result.display_label.product_context_label == "H3C–CH3"
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


def test_detailed_label_uses_superscript_fragment_indices() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1"
    )

    assert result.reaction_label == "Ar1–Br + Ar2–OH → Ar1–O–Ar2"
    assert result.display_label is not None
    assert result.display_label.detailed == (
        "Ar¹–Br + Ar²–OH → Ar¹–O–Ar²; edits: "
        "C–Br bond cleavage; C–O bond formation; O–H loss"
    )


def test_detailed_label_uses_subscripts_for_molecular_formula_counts() -> None:
    result = featurize_reaction("CC(=O)O.CN>>CC(=O)NC")

    assert result.reaction_label == (
        "R1–C(O)OH + R2–NH2 → R1–C(O)–NH–R2"
    )
    assert result.display_label is not None
    assert result.display_label.detailed == (
        "R¹–C(O)OH + R²–NH₂ → R¹–C(O)–NH–R²; edits: "
        "C–O bond cleavage; C–N bond formation; N–H loss"
    )


def test_formula_subscripts_do_not_change_maps_ring_sizes_or_edit_counts() -> None:
    mapped = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    ring = featurize_reaction(
        "OCCCCc1ccccc1Br>>c1ccc2c(c1)CCCCO2"
    )

    assert mapped.display_label is not None
    assert mapped.display_label.detailed.startswith("H₂C=CH₂ → H₃C–CH₃")
    assert "C(map 1)=C(map 2)" in mapped.display_label.detailed
    assert ring.display_label is not None
    assert "(7-membered ring)" in ring.display_label.detailed


def test_detailed_label_puts_seven_membered_ring_context_after_transformation() -> None:
    result = featurize_reaction(
        "OCCCCc1ccccc1Br>>c1ccc2c(c1)CCCCO2"
    )

    assert result.reaction_label == (
        "intramolecular (7-membered ring) Ar–Br / R–OH → Ar–O–R"
    )
    assert result.display_label is not None
    assert result.display_label.detailed == (
        "Ar–Br / R–OH → Ar–O–R; intramolecular (7-membered ring); "
        "edits: C–Br bond cleavage; C–O bond formation; O–H loss"
    )


def test_ascii_detailed_label_keeps_plain_fragment_indices() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
        label_style="ascii",
    )

    assert result.display_label is not None
    assert result.display_label.detailed.startswith(
        "Ar1-Br + Ar2-OH -> Ar1-O-Ar2; edits:"
    )


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

    assert payload["display_label"]["schema_version"] == "1.3"
    assert payload["display_label"]["clauses"][0]["edit_type"] == "formed"
    assert payload["display_label"]["clauses"][0]["atom_map_numbers"] == (1, 2)


def test_unknown_mapped_substitution_receives_generic_pattern_label() -> None:
    result = featurize_reaction(
        "[CH3:1][O:2][CH3:5].[NH2:3]>>"
        "[CH3:1][NH:3].[O-:2][CH3:5]"
    )

    assert result.named_family is None
    assert result.selected_candidate is None
    assert result.reaction_label == "C–N substitution"
    assert result.reaction_label_status == "mapped_generic_pattern"
    assert result.display_label is not None
    assert result.display_label.pattern_id == "substitution"
    assert result.display_label.grammar_id is None
    assert result.display_label.contextual_label is None
    assert result.display_label.structural_label == (
        "C–O bond cleavage; C–N bond formation; N–H loss"
    )


def test_unknown_mapped_dehydrogenation_receives_generic_pattern_label() -> None:
    result = featurize_reaction("[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]")

    assert result.named_family is None
    assert result.reaction_label == "H3C–CH3 → H2C=CH2"
    assert result.display_label is not None
    assert result.display_label.pattern_id == "dehydrogenation"
    assert result.display_label.transformation_label == (
        "C=C formation by dehydrogenation"
    )


def test_unknown_intramolecular_formation_receives_ring_closure_pattern() -> None:
    result = featurize_reaction(
        "[CH3:1][CH2:2][NH2:3]>>[CH2:1]1[CH2:2][NH:3]1"
    )

    assert result.named_family is None
    assert result.reaction_label == (
        "intramolecular (3-membered ring) C–N bond formation"
    )
    assert result.display_label is not None
    assert result.display_label.pattern_id == "intramolecular_bond_formation"
    assert result.display_label.structural_label == (
        "C–N bond formation; C–H loss; N–H loss"
    )


def test_exact_grammar_label_records_grammar_overlay_provenance() -> None:
    result = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")

    assert result.display_label is not None
    assert result.display_label.grammar_id == "sp2_c_n_substitution"
    assert result.display_label.grammar_label == result.reaction_label
    assert result.display_label.pattern_id == "substitution"


def test_reaction_label_definition_is_versioned() -> None:
    rendering = load_reaction_label_rendering()
    patterns = load_reaction_label_patterns()

    assert rendering["schema_version"] == "1.3"
    assert rendering["label_schema_version"] == "1.3"
    assert patterns["schema_version"] == "1.0"
    assert {pattern["id"] for pattern in patterns["patterns"]} >= {
        "substitution",
        "hydrogenation",
        "dehydrogenation",
        "reductive_bond_cleavage",
        "intramolecular_bond_formation",
    }
    assert set(rendering["clause_order"]) == {
        "formed",
        "broken",
        "order_changed",
        "hydrogen_change",
    }


def test_mapped_heteroatom_bond_reductions_receive_generic_labels() -> None:
    examples = {
        "c1ccc([N:1]=[N:2]c2ccccc2)cc1>>"
        "[NH2:1]c1ccccc1.[NH2:2]c1ccccc1": (
            "N=N reductive cleavage",
            "HB|Azo",
        ),
        "C[S:1][S:2]C>>C[SH:1].C[SH:2]": (
            "S–S reductive cleavage",
            "HB|Disulfide",
        ),
        "C[O:1][O:2]C>>C[OH:1].C[OH:2]": (
            "O–O reductive cleavage",
            "HB|Peroxide",
        ),
    }

    for reaction, (label, handle_signature) in examples.items():
        result = featurize_reaction(reaction)
        assert result.valid
        assert result.reaction_label == label
        assert result.reaction_label_status == "mapped_generic_pattern"
        assert result.display_label is not None
        assert result.display_label.pattern_id == "reductive_bond_cleavage"
        assert any(
            site.canonical_signature == handle_signature
            for component in result.reactants
            for site in component.compound_analysis.sites
        )
