"""Regression tests for conservative unmapped scaffold correspondence."""

from reactive_taxonomy import featurize_reaction


def test_unmapped_alkene_hydrogenation_receives_inferred_label() -> None:
    result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )

    assert result.valid
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.reaction_label.concise == "Ar–CH=CH2 → Ar–CH2–CH3"
    assert result.reaction_label.status == "exact_reconstruction"
    assert result.reaction_label is not None
    assert result.reaction_label.pattern_id == "hydrogenation"
    assert result.reaction_label.transformation_label == "C=C hydrogenation"
    assert result.reaction_label.reactant_context_label == "Ar–CH=CH2"
    assert result.reaction_label.product_context_label == "Ar–CH2–CH3"
    assert result.reaction_label.confidence == 1.0
    assert result.reaction_signature is not None
    assert result.reaction_signature.product_transformation is not None
    assert (
        result.reaction_signature.product_transformation.concise_label
        == "Ar–CH2–CH3"
    )
    assert "INFERRED_ATOM_CORRESPONDENCE" not in result.warnings


def test_unmapped_carbonyl_redox_receives_generic_labels() -> None:
    reduction = featurize_reaction(
        "c1ccc(C=O)cc1.[H][H]>>c1ccc(CO)cc1"
    )
    oxidation = featurize_reaction(
        "CC(O)c1ccccc1.[O]>>CC(=O)c1ccccc1"
    )

    assert reduction.reaction_label.concise == "Ar–CH=O → Ar–CH2–OH"
    assert reduction.reaction_label is not None
    assert reduction.reaction_label.pattern_id == "heteroatom_bond_reduction"
    assert reduction.reaction_label.transformation_label == "C=O reduction"
    assert oxidation.reaction_label.concise == "Ar–CH(R)–OH → Ar–C(R)=O"
    assert oxidation.reaction_label is not None
    assert oxidation.reaction_label.pattern_id == "heteroatom_bond_oxidation"


def test_unmapped_complete_alkyne_hydrogenation_receives_generic_label() -> None:
    result = featurize_reaction(
        "C#Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )

    assert result.reaction_label.concise == "Ar–C≡CH → Ar–CH2–CH3"
    assert result.reaction_label is not None
    assert result.reaction_label.pattern_id == "complete_alkyne_hydrogenation"
    assert result.reaction_label.structural_label == (
        "C≡C → C–C; 4 × H gain at C"
    )


def test_symmetry_equivalent_correspondence_is_order_invariant() -> None:
    forward = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )
    reversed_order = featurize_reaction(
        "[H][H].C=Cc1ccccc1>>CCc1ccccc1"
    )

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )
    assert forward.reaction_label == reversed_order.reaction_label


def test_contextual_label_style_changes_display_only() -> None:
    unicode_result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1"
    )
    ascii_result = featurize_reaction(
        "C=Cc1ccccc1.[H][H]>>CCc1ccccc1", label_style="ascii"
    )

    assert unicode_result.reaction_label.concise == "Ar–CH=CH2 → Ar–CH2–CH3"
    assert ascii_result.reaction_label.concise == "Ar-CH=CH2 -> Ar-CH2-CH3"
    assert unicode_result.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert unicode_result.reaction_signature.signature_id == (
        ascii_result.reaction_signature.signature_id
    )


def test_multiple_order_changes_render_as_repeated_events() -> None:
    result = featurize_reaction(
        "[CH2:1]=[CH:2][CH:3]=[CH2:4]>>"
        "[CH3:1][CH2:2][CH2:3][CH3:4]"
    )

    assert result.reaction_label is not None
    assert result.reaction_label.contextual_label is None
    assert result.reaction_label.pattern_id is None
    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert result.reaction_label.concise == "2 × C=C hydrogenation"
    assert result.reaction_label.status == "multi_event"


def test_beta_elimination_resolves_previously_ambiguous_correspondence() -> None:
    result = featurize_reaction("CCCBr>>CC=C")

    assert result.valid
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate is not None
    assert result.selected_candidate.annotation_id == "beta_halo_elimination"
    assert result.edit_archetype == "elimination"
    assert result.reaction_signature is not None
    assert result.reaction_signature.edit_archetype == "elimination"


def test_unresolvable_multi_substrate_assembly_remains_unresolved() -> None:
    result = featurize_reaction("CC.CN>>CCN")

    assert result.valid
    assert result.evidence_quality == "unresolved"
    assert result.reaction_label.status == "unavailable"
    assert result.reaction_signature is None
    assert "GLOBAL_CORRESPONDENCE_NOT_FOUND" in result.warnings


def test_global_correspondence_recovers_aldol_graph_edits() -> None:
    result = featurize_reaction(
        "CC=O.CC=O>>CC(O)CC=O",
        label_style="ascii",
    )

    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_signature is not None
    assert result.transformation_class == "generic_graph_transformation"
    assert result.reaction_completeness is not None
    assert result.reaction_completeness.status == "verified"
    assert all(
        edit.confidence == 0.8 for edit in result.reaction_signature.edits
    )
    assert result.reaction_label.concise == "C-H + C=O -> C-C + C-O + O-H"
    assert result.reaction_label.status == "observed_edits"
    assert "INFERRED_GLOBAL_ATOM_CORRESPONDENCE" in result.warnings


def test_global_correspondence_recovers_cycloaddition_graph_edits() -> None:
    result = featurize_reaction(
        "CC#C.CN=[N+]=[N-]>>Cc1nnn(C)c1",
        label_style="ascii",
    )

    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_signature is not None
    assert result.transformation_class == "generic_graph_transformation"
    assert result.reaction_signature.topology.reaction_scope == "intermolecular"
    assert result.reaction_label is not None
    assert "->" in result.reaction_label.concise
    assert result.reaction_label.concise == "C#C + N=N=N -> aromatic 5-membered C2N3 ring"
    assert result.reaction_label.status == "ring_formation"
    assert result.reaction_signature.topology.formed_ring_sizes == (5,)


def test_global_correspondence_recovers_condensation_graph_edits() -> None:
    result = featurize_reaction(
        "CC=O.CN>>CC=NC",
        label_style="ascii",
    )

    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_signature is not None
    assert result.transformation_class == "substitution"
    assert result.reaction_label.concise == "C=O + 2 x N-H -> C=N"
    assert result.reaction_signature.formed_bond_types == ("C-N:DOUBLE",)
    assert result.reaction_signature.broken_bond_types == ("C-O:DOUBLE",)
    assert result.reaction_signature.order_changes == ()
    assert {
        edit.edit_type for edit in result.reaction_signature.edits
    } == {"broken", "formed", "hydrogen_change"}


def test_unique_scaffold_edits_receive_generic_transformation_class() -> None:
    result = featurize_reaction(
        "CCCCCCCCCCCCOC(=S)SC>>CCCCCCCCCCCC",
        label_style="ascii",
    )

    assert result.evidence_quality == "unique_scaffold_correspondence"
    assert result.reaction_signature is not None
    assert result.transformation_class == "generic_graph_transformation"


def test_explicit_atom_stereo_is_retained_in_signature_identity() -> None:
    retained = featurize_reaction(
        "C[C@H](O)C=O.CN>>C[C@H](O)C=NC",
        label_style="ascii",
    )
    unspecified = featurize_reaction(
        "C[C@H](O)C=O.CN>>CC(O)C=NC",
        label_style="ascii",
    )

    assert retained.evidence_quality == "global_atom_correspondence"
    assert retained.reaction_signature is not None
    assert unspecified.reaction_signature is not None
    assert any(
        change.stereo_type == "atom"
        and change.change_type == "retained"
        and change.old_descriptor == change.new_descriptor
        for change in retained.reaction_signature.stereo_changes
    )
    assert retained.reaction_signature.product_transformation is not None
    assert (
        retained.reaction_signature.product_transformation.stereo_changes
        == retained.reaction_signature.stereo_changes
    )
    assert (
        retained.reaction_signature.signature_id
        != unspecified.reaction_signature.signature_id
    )


def test_global_correspondence_is_reactant_order_invariant() -> None:
    forward = featurize_reaction("CC=O.CN>>CC=NC")
    reversed_order = featurize_reaction("CN.CC=O>>CC=NC")

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert forward.reaction_signature.signature_id == (
        reversed_order.reaction_signature.signature_id
    )
    assert forward.reaction_label == reversed_order.reaction_label


def test_fragmented_correspondence_recovers_topology_changing_scaffold() -> None:
    result = featurize_reaction("CNCCCN>>CN1CCC1", label_style="ascii")

    assert result.evidence_quality == "fragmented_scaffold_correspondence"
    assert result.reaction_signature is not None
    assert result.reaction_signature.named_family is None
    assert result.reaction_signature.topology.ring_count_delta == 1
    assert result.reaction_signature.formed_bond_types == ("C-N:SINGLE",)
    assert result.reaction_signature.broken_bond_types == ("C-N:SINGLE",)
    assert "INFERRED_FRAGMENTED_SCAFFOLD_CORRESPONDENCE" in result.warnings


def test_fragmented_correspondence_is_smiles_serialization_invariant() -> None:
    forward = featurize_reaction("CNCCCN>>CN1CCC1")
    reversed_smiles = featurize_reaction("NCCCNC>>C1CCN(C)1")

    assert forward.reaction_signature is not None
    assert reversed_smiles.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_smiles.reaction_signature.signature_id
    )


def test_fragmented_correspondence_abstains_for_competing_nitrogen_origins() -> None:
    result = featurize_reaction(
        "CCOC(=O)/C(CC)=N\\Nc1c(Cl)cccc1Cl"
        ">>CCOC(=O)c1[nH]c2c(Cl)cc(Cl)cc2c1C"
    )

    assert result.evidence_quality == "ambiguous_atom_correspondence"
    assert result.reaction_signature is None
    assert any(
        warning.startswith("AMBIGUOUS_SCAFFOLD_CORRESPONDENCE:")
        for warning in result.warnings
    )
