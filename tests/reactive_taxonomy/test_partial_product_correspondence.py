"""Regressions for incomplete but structurally informative reactions."""

from reactive_taxonomy import featurize_reaction


def test_missing_chlorine_source_yields_partial_acyl_substitution() -> None:
    result = featurize_reaction(
        "O=C(O)c1cccc(I)c1>>O=C(Cl)c1cccc(I)c1"
    )

    assert result.valid
    assert result.reaction_label == (
        "R–C(=O)–OH → R–C(=O)–Cl [Cl source missing]"
    )
    assert result.reaction_label_status == "partial_product_correspondence"
    assert result.evidence_quality == "partial_product_correspondence"
    assert result.transformation_class == "acyl_heteroatom_substitution"
    assert result.named_family is None
    assert result.reaction_signature is None
    assert result.reaction_completeness.status == "incomplete"
    assert result.display_label.status == "partial_product_correspondence"
    assert result.display_label.detailed == (
        "R–C(=O)–OH → R–C(=O)–Cl; "
        "partial conserved-scaffold observation; "
        "the reactants do not account for Cl in the product."
    )
    assert "Ar–I" not in result.reaction_label

    observation = result.partial_product_transformation
    assert observation is not None
    assert observation.transformation_type == "attachment_replacement"
    assert observation.reactant_center.element == "C"
    assert observation.removed_attachment.element == "O"
    assert observation.added_attachment.element == "Cl"
    assert observation.missing_product_atom_elements == ("Cl",)
    assert observation.product_heavy_atom_coverage == 0.9
    assert "PRODUCT_ATOM_SOURCE_UNRESOLVED:Cl" in result.warnings
    assert "PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES" in result.warnings
    assert all(
        candidate.verification
        in {"construction_failed", "product_mismatch"}
        for candidate in result.candidates
    )


def test_partial_attachment_replacement_is_not_acyl_specific() -> None:
    result = featurize_reaction("CC(C)O>>CC(C)Cl")

    assert result.reaction_label == "C–OH → C–Cl [Cl source missing]"
    assert result.transformation_class == "attachment_replacement"
    assert result.partial_product_transformation is not None
    assert result.reaction_signature is None
    assert result.display_label.detailed.startswith(
        "C–OH → C–Cl; partial conserved-scaffold observation;"
    )


def test_partial_replacement_ascii_label_is_deterministic() -> None:
    result = featurize_reaction(
        "CC(=O)O>>CC(=O)F",
        label_style="ascii",
    )

    assert result.reaction_label == (
        "R-C(=O)-OH -> R-C(=O)-F [F source missing]"
    )
    assert result.partial_product_transformation is not None
    assert result.display_label.detailed.startswith(
        "R-C(=O)-OH -> R-C(=O)-F; partial conserved-scaffold observation;"
    )


def test_product_growth_without_attachment_exchange_is_not_forced() -> None:
    result = featurize_reaction("CC>>CCC")

    assert result.partial_product_transformation is None
    assert result.reaction_signature is None
