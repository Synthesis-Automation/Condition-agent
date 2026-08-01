from reactive_taxonomy import (
    featurize_reaction,
    load_reaction_label_patterns,
    load_reaction_label_rendering,
)


def test_mapped_unknown_reaction_receives_observed_edit_label() -> None:
    result = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")

    assert result.named_family is None
    assert result.reaction_label.concise == "C–N bond formation"
    assert result.reaction_label.status == "observed_edits"
    assert result.reaction_label is not None
    assert result.reaction_label.status == "observed_edits"
    assert result.reaction_label.detailed == (
        "C–N bond formation; edits: C(map 1)–N(map 2) bond formation"
    )
    assert len(result.reaction_label.clauses) == 1
    assert result.reaction_label.clauses[0].evidence == "supplied_atom_mapping"


def test_verified_core_enriches_generic_observed_edit_label() -> None:
    result = featurize_reaction(
        "[CH3:1][S:2][CH3:3].[O:4]>>[CH3:1][S:2](=[O:4])[CH3:3]"
    )

    assert result.reaction_label.concise == "S(R)2 + O → S(R)2(=O)"
    assert result.reaction_label.status == "core_projection"
    assert result.reaction_label.source == "reaction_core"
    assert result.reaction_label.structural_label == "S(R)2 + O → S(R)2(=O)"
    assert "retained context: R (C)" in result.reaction_label.detailed
    assert "literal edits: O(map 4)=S(map 2) bond formation" in (
        result.reaction_label.detailed
    )


def test_simple_mapped_formation_keeps_clearer_literal_edit_label() -> None:
    result = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")

    assert result.reaction_core is not None
    assert result.reaction_core.generic_label == "∅ → C–N"
    assert result.reaction_label.concise == "C–N bond formation"
    assert result.reaction_label.source == "literal_edits"


def test_multiple_edits_compose_and_collapse_repeated_generic_clauses() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")

    assert result.reaction_label.concise == "H2C=CH2 → H3C–CH3"
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_label is not None
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_label.pattern_id == "hydrogenation"
    assert result.reaction_label.structural_label == "C=C → C–C; 2 × H gain at C"
    assert result.reaction_label.transformation_label == "C=C hydrogenation"
    assert result.reaction_label.contextual_label == "H2C=CH2 → H3C–CH3"
    assert result.reaction_label.product_context_label == "H3C–CH3"
    assert len(result.reaction_label.clauses) == 3
    assert "H gain at C(map 1)" in result.reaction_label.detailed
    assert "H gain at C(map 2)" in result.reaction_label.detailed


def test_generic_label_is_invariant_to_reactant_component_order() -> None:
    forward = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")
    reversed_order = featurize_reaction("[NH2:2].[CH3:1]>>[CH3:1][NH2:2]")

    assert forward.reaction_label == reversed_order.reaction_label
    assert forward.reaction_label is not None
    assert reversed_order.reaction_label is not None
    assert forward.reaction_label.detailed == reversed_order.reaction_label.detailed


def test_ascii_rendering_changes_display_only() -> None:
    unicode_result = featurize_reaction("[CH3:1].[NH2:2]>>[CH3:1][NH2:2]")
    ascii_result = featurize_reaction(
        "[CH3:1].[NH2:2]>>[CH3:1][NH2:2]", label_style="ascii"
    )

    assert unicode_result.reaction_label.concise == "C–N bond formation"
    assert ascii_result.reaction_label.concise == "C-N bond formation"
    assert unicode_result.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert (
        unicode_result.reaction_signature.signature_id
        == ascii_result.reaction_signature.signature_id
    )


def test_exact_grammar_label_is_preserved_as_overlay() -> None:
    result = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")

    assert result.reaction_label.concise == "Ar–Br + R–NH2 → Ar–NH–R"
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_label is not None
    assert result.reaction_label.status == "exact_reconstruction"
    assert "C–Br bond cleavage" in result.reaction_label.detailed
    assert "C–N bond formation" in result.reaction_label.detailed
    assert "N–H loss" in result.reaction_label.detailed


def test_detailed_label_uses_superscript_fragment_indices() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1"
    )

    assert result.reaction_label.concise == "Ar1–Br + Ar2–OH → Ar1–O–Ar2"
    assert result.reaction_label is not None
    assert result.reaction_label.detailed == (
        "Ar¹–Br + Ar²–OH → Ar¹–O–Ar²; edits: "
        "C–Br bond cleavage; C–O bond formation; O–H loss"
    )


def test_detailed_label_uses_subscripts_for_molecular_formula_counts() -> None:
    result = featurize_reaction("CC(=O)O.CN>>CC(=O)NC")

    assert result.reaction_label.concise == (
        "R1–C(O)OH + R2–NH2 → R1–C(O)–NH–R2"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.detailed == (
        "R¹–C(O)OH + R²–NH₂ → R¹–C(O)–NH–R²; edits: "
        "C–O bond cleavage; C–N bond formation; N–H loss"
    )


def test_formula_subscripts_do_not_change_maps_ring_sizes_or_edit_counts() -> None:
    mapped = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    ring = featurize_reaction(
        "OCCCCc1ccccc1Br>>c1ccc2c(c1)CCCCO2"
    )

    assert mapped.reaction_label is not None
    assert mapped.reaction_label.detailed.startswith("H₂C=CH₂ → H₃C–CH₃")
    assert "C(map 1)=C(map 2)" in mapped.reaction_label.detailed
    assert ring.reaction_label is not None
    assert "(7-membered ring)" in ring.reaction_label.detailed


def test_detailed_label_puts_seven_membered_ring_context_after_transformation() -> None:
    result = featurize_reaction(
        "OCCCCc1ccccc1Br>>c1ccc2c(c1)CCCCO2"
    )

    assert result.reaction_label.concise == (
        "intramolecular (7-membered ring) Ar–Br / R–OH → Ar–O–R"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.detailed == (
        "Ar–Br / R–OH → Ar–O–R; intramolecular (7-membered ring); "
        "edits: C–Br bond cleavage; C–O bond formation; O–H loss"
    )


def test_ascii_detailed_label_keeps_plain_fragment_indices() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
        label_style="ascii",
    )

    assert result.reaction_label is not None
    assert result.reaction_label.detailed.startswith(
        "Ar1-Br + Ar2-OH -> Ar1-O-Ar2; edits:"
    )


def test_conflicting_evidence_is_visible_in_label_contract() -> None:
    reaction = (
        "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[NH2:8][CH3:9]>>"
        "[cH:2]1[c:3]([NH:8][CH3:9])[cH:4][cH:5][cH:6][cH:7]1"
    )
    result = featurize_reaction(reaction)

    assert result.reaction_label.status == "conflicting_evidence"
    assert result.reaction_label is not None
    assert result.reaction_label.concise.startswith("Conflicting evidence:")
    assert result.reaction_label is not None
    assert result.reaction_label.status == "conflicting_evidence"
    assert result.reaction_label.confidence == 0.5
    assert result.reaction_label.detailed.startswith(
        result.reaction_label.structural_label + "; conflicting evidence;"
    )
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.reaction_label.warnings


def test_display_label_serializes_as_nested_evidence() -> None:
    payload = featurize_reaction(
        "[CH3:1].[NH2:2]>>[CH3:1][NH2:2]"
    ).to_dict()

    assert payload["reaction_label"]["schema_version"] == "2.0"
    assert payload["reaction_label"]["source"] == "literal_edits"
    assert payload["reaction_label"]["clauses"][0]["edit_type"] == "formed"
    assert payload["reaction_label"]["clauses"][0]["atom_map_numbers"] == (1, 2)


def test_unknown_mapped_substitution_combines_core_and_pattern_label() -> None:
    result = featurize_reaction(
        "[CH3:1][O:2][CH3:5].[NH2:3]>>"
        "[CH3:1][NH:3].[O-:2][CH3:5]"
    )

    assert result.named_family is None
    assert result.selected_candidate is None
    assert result.reaction_label.concise == (
        "C(H)3(O-R) + N(H)2 → C(H)3(N-H)"
    )
    assert result.reaction_label.status == "core_projection"
    assert result.reaction_label.source == "reaction_core"
    assert result.reaction_label is not None
    assert result.reaction_label.pattern_id == "substitution"
    assert result.reaction_label.grammar_id is None
    assert result.reaction_label.contextual_label is None
    assert result.reaction_label.structural_label == (
        "C(H)3(O-R) + N(H)2 → C(H)3(N-H)"
    )
    assert "transformation pattern: C–N substitution" in (
        result.reaction_label.detailed
    )


def test_exact_grammar_keeps_concise_label_and_adds_core_context() -> None:
    result = featurize_reaction(
        "[O:1]=[C:2]1[CH2:3][CH2:4][CH2:5][CH2:6][CH2:7]1>>"
        "[OH:1][CH:2]1[CH2:3][CH2:4][CH2:5][CH2:6][CH2:7]1"
    )

    assert result.reaction_label.concise == "R–C(R)=O → R–CH(R)–OH"
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_label.source == "verified_grammar"
    assert "core projection: C(Cycloalkyl)₂(=O)" in (
        result.reaction_label.detailed
    )
    assert "retained context: Cycloalkyl" in result.reaction_label.detailed


def test_contradicted_reactant_fallback_never_ends_in_arrow() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"
    )

    assert result.reaction_label.status == "product_contradicted_reactants"
    assert not result.reaction_label.concise.endswith("→")
    assert result.reaction_label.concise.endswith("[product contradicted]")
    assert "no candidate product transformation is asserted" in (
        result.reaction_label.detailed
    )


def test_unknown_mapped_dehydrogenation_receives_generic_pattern_label() -> None:
    result = featurize_reaction("[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]")

    assert result.named_family is None
    assert result.reaction_label.concise == "H3C–CH3 → H2C=CH2"
    assert result.reaction_label is not None
    assert result.reaction_label.pattern_id == "dehydrogenation"
    assert result.reaction_label.transformation_label == (
        "C=C formation by dehydrogenation"
    )


def test_unknown_intramolecular_formation_receives_ring_closure_pattern() -> None:
    result = featurize_reaction(
        "[CH3:1][CH2:2][NH2:3]>>[CH2:1]1[CH2:2][NH:3]1"
    )

    assert result.named_family is None
    assert result.reaction_label.concise == (
        "intramolecular (3-membered ring) C–N bond formation"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.pattern_id == "intramolecular_bond_formation"
    assert result.reaction_label.structural_label == (
        "C–N bond formation; C–H loss; N–H loss"
    )


def test_exact_grammar_label_records_grammar_overlay_provenance() -> None:
    result = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")

    assert result.reaction_label is not None
    assert result.reaction_label.grammar_id == "sp2_c_n_substitution"
    assert result.reaction_label.grammar_label == result.reaction_label.concise
    assert result.reaction_label.pattern_id == "substitution"


def test_reaction_label_definition_is_versioned() -> None:
    rendering = load_reaction_label_rendering()
    patterns = load_reaction_label_patterns()

    assert rendering["schema_version"] == "2.0"
    assert rendering["label_schema_version"] == "2.1"
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
        assert result.reaction_label.concise == label
        assert result.reaction_label.status == "generic_pattern"
        assert result.reaction_label is not None
        assert result.reaction_label.pattern_id == "reductive_bond_cleavage"
        assert any(
            site.canonical_signature == handle_signature
            for component in result.reactants
            for site in component.compound_analysis.sites
        )
