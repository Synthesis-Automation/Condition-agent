from reactive_taxonomy import featurize_reaction


def test_supported_exact_families_produce_signatures() -> None:
    cases = {
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1": "C-C:SINGLE",
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1": "C-N:SINGLE",
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1": "C-O:SINGLE",
        "Brc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1": "C-S:SINGLE",
        "CC(=O)Cl.CN>>CC(=O)NC": "C-N:SINGLE",
    }

    for reaction, formed_bond in cases.items():
        result = featurize_reaction(reaction)
        assert result.reaction_signature is not None, reaction
        assert formed_bond in result.reaction_signature.formed_bond_types
        assert result.reaction_signature.product_transformation is not None
        assert result.reaction_signature.product_transformation.exact_product_verified


def test_signature_is_stable_across_reactant_order() -> None:
    forward = featurize_reaction("Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1")
    reversed_order = featurize_reaction(
        "c1ccc(B(O)O)cc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_order.reaction_signature.signature_id
    )
    assert (
        forward.reaction_signature.exact_signature_key
        == reversed_order.reaction_signature.exact_signature_key
    )
    assert forward.fallback_descriptor is not None
    assert reversed_order.fallback_descriptor is not None
    assert (
        forward.fallback_descriptor.descriptor_id
        == reversed_order.fallback_descriptor.descriptor_id
    )


def test_signature_identity_is_display_style_independent() -> None:
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    unicode_result = featurize_reaction(reaction, label_style="unicode")
    ascii_result = featurize_reaction(reaction, label_style="ascii")

    assert unicode_result.reaction_label != ascii_result.reaction_label
    assert unicode_result.reaction_signature is not None
    assert ascii_result.reaction_signature is not None
    assert (
        unicode_result.reaction_signature.signature_id
        == ascii_result.reaction_signature.signature_id
    )


def test_signature_contains_versioned_hierarchical_keys() -> None:
    result = featurize_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")
    signature = result.reaction_signature

    assert signature is not None
    assert signature.signature_id.startswith("RS3:")
    assert signature.exact_signature_key.startswith("L0:")
    assert signature.handle_signature_key.startswith("L1:")
    assert signature.transformation_signature_key.startswith("L2:")
    assert signature.bond_edit_signature_key.startswith("L3:")
    assert signature.environment_signature_key.startswith("L4:")
    assert "signature_features.v3.json" in signature.definition_versions
    assert "reactivity_descriptor_rules.v1.json" in signature.definition_versions


def test_signature_partners_retain_typed_reactivity_profiles() -> None:
    cases = {
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1": (
            "Ar",
            "aromatic",
            "open",
            "balanced",
        ),
        "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1": (
            "HeteroAr",
            "aromatic",
            "open",
            "electron_poor",
        ),
        "CCBr.Nc1ccccc1>>CCNc1ccccc1": ("Alkyl", "alkyl", "open", "balanced"),
        "CC(C)Br.N>>CC(C)N": ("Alkyl", "alkyl", "hindered", "balanced"),
    }

    for reaction, expected in cases.items():
        result = featurize_reaction(reaction)
        assert result.reaction_signature is not None
        electrophile = next(
            partner
            for partner in result.reaction_signature.partners
            if partner.role == "electrophile"
        )
        assert electrophile.anchor_contexts == (expected[0],)
        assert electrophile.reactivity_profile is not None
        assert electrophile.reactivity_profile.context_kind == expected[1]
        assert electrophile.reactivity_profile.steric.accessibility_class == expected[2]
        assert (
            electrophile.reactivity_profile.electronic.activation_class == expected[3]
        )


def test_signature_serializes_with_analysis() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    payload = result.to_dict()

    assert payload["reaction_signature"]["named_family"] is None
    assert payload["reaction_signature"]["order_changes"] == ("C-C:DOUBLE>SINGLE",)
    assert payload["reaction_signature"]["hydrogen_changes"] == (
        "C-H:NONE>SINGLE",
        "C-H:NONE>SINGLE",
    )
    assert payload["schema_version"] == "3.3"
    assert payload["reaction_signature"]["schema_version"] == "3.0"
    assert payload["reaction_signature"]["topology"]["reaction_scope"] == (
        "unimolecular"
    )


def test_mixed_c_o_c_s_reaction_is_partitioned_into_two_events() -> None:
    reaction = (
        "[OH:1][CH3:2].[SH:3][CH3:4]."
        "[c:5]1([F:6])[cH:7][cH:8][c:9]([F:10])[cH:11][cH:12]1"
        ">>[c:5]1([O:1][CH3:2])[cH:7][cH:8]"
        "[c:9]([S:3][CH3:4])[cH:11][cH:12]1"
    )
    result = featurize_reaction(reaction, label_style="ascii")
    signature = result.reaction_signature

    assert signature is not None
    assert signature.event_count == 2
    assert signature.event_scope == "multi_event"
    assert [event.formed_bond_types for event in signature.events] == [
        ("C-O:SINGLE",),
        ("C-S:SINGLE",),
    ]
    assert all(event.broken_bond_types == ("C-F:SINGLE",) for event in signature.events)
    assert [relation.relation_type for relation in signature.event_relations] == [
        "shared_component"
    ]
    assert result.reaction_label == "C-O substitution + C-S substitution"
    assert result.reaction_label_status == "multi_event_edit_summary"
    assert result.display_label is not None
    assert result.display_label.event_labels == (
        "C-O substitution",
        "C-S substitution",
    )
    assert result.display_label.detailed.startswith(
        "C-O substitution + C-S substitution; events:"
    )
    assert "Event 1:" in result.display_label.detailed


def test_multi_event_signature_is_reactant_order_invariant() -> None:
    substrate = "[c:5]1([F:6])[cH:7][cH:8][c:9]([F:10])[cH:11][cH:12]1"
    product = "[c:5]1([O:1][CH3:2])[cH:7][cH:8][c:9]([S:3][CH3:4])[cH:11][cH:12]1"
    forward = featurize_reaction(f"[OH:1][CH3:2].[SH:3][CH3:4].{substrate}>>{product}")
    reversed_order = featurize_reaction(
        f"{substrate}.[SH:3][CH3:4].[OH:1][CH3:2]>>{product}"
    )

    assert forward.reaction_signature is not None
    assert reversed_order.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_order.reaction_signature.signature_id
    )


def test_repeated_events_are_counted_in_display_label() -> None:
    reaction = (
        "[SH:1][CH3:2].[SH:3][CH3:4]."
        "[c:5]1([F:6])[cH:7][cH:8][c:9]([F:10])[cH:11][cH:12]1"
        ">>[c:5]1([S:1][CH3:2])[cH:7][cH:8]"
        "[c:9]([S:3][CH3:4])[cH:11][cH:12]1"
    )
    result = featurize_reaction(reaction, label_style="ascii")

    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert result.reaction_label == "2 x C-S substitution"
    assert result.display_label is not None
    assert result.display_label.detailed.startswith("2 x C-S substitution; events:")


def test_balanced_unmapped_multi_event_reaction_is_exactly_reconstructed() -> None:
    reaction = "CO.CS.Fc1ccc(F)cc1>>COc1ccc(SC)cc1"
    result = featurize_reaction(reaction, label_style="ascii")

    assert result.evidence_quality == "exact_multi_event_reconstruction"
    assert len(result.selected_events) == 2
    assert result.selected_candidate is None
    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert result.reaction_signature.product_transformation is not None
    assert result.reaction_signature.product_transformation.exact_product_verified
    assert result.transformation_class == "generic_multi_event_graph_transformation"
    assert result.reaction_label == "C-O substitution + C-S substitution"
    assert [
        event.transformation_class for event in result.reaction_signature.events
    ] == [
        "sp2_c_o_substitution",
        "sp2_c_s_substitution",
    ]


def test_overlapping_multi_event_candidates_do_not_use_removed_join_atoms() -> None:
    reaction = (
        "O=C(O)c1ccccc1."
        "CN(C)c1ccncc1[C@H](O)[C@@H](N)C(C)(C)C."
        "CC(=O)OC(C)=O"
        ">>CC(=O)O[C@@H](c1cnccc1N(C)C)"
        "[C@@H](NC(=O)c1ccccc1)C(C)(C)C"
    )

    result = featurize_reaction(reaction, label_style="ascii")

    assert result.valid
    assert result.evidence_quality == "exact_multi_event_reconstruction"
    assert result.reaction_label == "C-O substitution + C-N substitution"
    assert len(result.selected_events) == 2


def test_unbalanced_multi_event_reaction_does_not_invent_partner_copy() -> None:
    reaction = (
        "Sc1ccccc1.O=[N+]([O-])c1c(F)cccc1F>>O=[N+]([O-])c1c(Sc2ccccc2)cccc1Sc1ccccc1"
    )
    result = featurize_reaction(reaction, label_style="ascii")

    assert result.evidence_quality == "reactant_grammar_only"
    assert result.selected_events == ()
    assert result.reaction_signature is None
    assert result.reaction_label is None
    assert result.reaction_label_status == "product_contradicted_candidates"
    assert "PRODUCT_CONTRADICTED_GRAMMAR_CANDIDATES" in result.warnings


def test_l3_identity_includes_schema_level_hydrogen_changes() -> None:
    hydrogenation = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    order_change_only = featurize_reaction("[CH2:1]=[CH2:2]>>[CH2:1][CH2:2]")

    assert hydrogenation.reaction_signature is not None
    assert order_change_only.reaction_signature is not None
    assert hydrogenation.reaction_signature.order_changes == (
        order_change_only.reaction_signature.order_changes
    )
    assert hydrogenation.reaction_signature.hydrogen_changes
    assert not order_change_only.reaction_signature.hydrogen_changes
    assert (
        hydrogenation.reaction_signature.bond_edit_signature_key
        != order_change_only.reaction_signature.bond_edit_signature_key
    )
