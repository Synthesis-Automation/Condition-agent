from reactive_taxonomy import featurize_reaction


def test_valid_mapping_produces_typed_order_change() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")

    assert result.evidence_quality == "validated_atom_mapping"
    assert result.reaction_signature is not None
    assert result.reaction_signature.named_family is None
    assert len(result.reaction_signature.edits) == 3
    edit = next(
        edit
        for edit in result.reaction_signature.edits
        if edit.edit_type == "order_changed"
    )
    assert edit.edit_type == "order_changed"
    assert edit.old_order == "DOUBLE"
    assert edit.new_order == "SINGLE"
    assert edit.atom_1.atom_map_number == 1
    assert edit.atom_2 is not None
    assert edit.atom_2.atom_map_number == 2
    hydrogen_edits = [
        edit
        for edit in result.reaction_signature.edits
        if edit.edit_type == "hydrogen_change"
    ]
    assert len(hydrogen_edits) == 2
    assert {edit.atom_1.atom_map_number for edit in hydrogen_edits} == {1, 2}
    assert all(edit.old_order is None for edit in hydrogen_edits)
    assert all(edit.new_order == "SINGLE" for edit in hydrogen_edits)


def test_mapped_dehydrogenation_produces_hydrogen_changes() -> None:
    result = featurize_reaction("[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]")

    assert result.reaction_signature is not None
    hydrogen_edits = [
        edit
        for edit in result.reaction_signature.edits
        if edit.edit_type == "hydrogen_change"
    ]
    assert len(hydrogen_edits) == 2
    assert all(edit.old_order == "SINGLE" for edit in hydrogen_edits)
    assert all(edit.new_order is None for edit in hydrogen_edits)


def test_duplicate_atom_map_does_not_upgrade_evidence() -> None:
    result = featurize_reaction(
        "[CH3:1][Br:1].[NH2:2]C>>[CH3:1][NH:2]C"
    )

    assert result.reaction_signature is None
    assert result.evidence_quality == "unresolved"
    assert "DUPLICATE_ATOM_MAP:reactant:1" in result.warnings


def test_atom_map_element_mismatch_is_rejected() -> None:
    result = featurize_reaction("[CH3:1]Br>>[NH3:1]")

    assert result.reaction_signature is None
    assert result.evidence_quality == "unresolved"
    assert "ATOM_MAP_ELEMENT_MISMATCH:1:C:N" in result.warnings


def test_mapping_and_reconstruction_conflict_is_preserved() -> None:
    reaction = (
        "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[NH2:8][CH3:9]>>"
        "[cH:2]1[c:3]([NH:8][CH3:9])[cH:4][cH:5][cH:6][cH:7]1"
    )
    result = featurize_reaction(reaction)

    assert result.selected_candidate is not None
    assert result.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in result.warnings
    assert result.reaction_signature is not None
    assert result.reaction_signature.evidence_quality == "conflicting_edit_evidence"
    assert result.reaction_signature.transformation_confidence == 0.5


def test_mapping_and_reconstruction_agreement_keeps_exact_compatibility_view() -> None:
    reaction = (
        "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1."
        "[NH2:8][CH3:9]>>"
        "[c:2]1([NH:8][CH3:9])[cH:3][cH:4][cH:5][cH:6][cH:7]1"
    )
    result = featurize_reaction(reaction)

    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.reaction_signature is not None
    assert (
        result.reaction_signature.evidence_quality
        == "validated_mapping_and_exact_reconstruction"
    )
    assert result.product_connection is not None


def test_scaffold_fallback_sorts_mixed_edit_orders_deterministically() -> None:
    reaction = (
        "C[C@@H]1CN2C(=O)[C@H]3[C@@H]4C=C5CCc6cccc(c65)"
        "[C@@]43[C@H]2CC1=O"
        ">>C[C@@H]1CN2C(=O)C3=C(c4cccc5c4[C@@H](CC3)CC5)"
        "[C@H]2C[C@H]1O"
    )

    first = featurize_reaction(reaction)
    second = featurize_reaction(reaction)

    assert first.valid
    assert first.to_dict() == second.to_dict()
