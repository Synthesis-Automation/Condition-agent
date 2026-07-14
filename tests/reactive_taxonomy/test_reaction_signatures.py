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
    forward = featurize_reaction(
        "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    )
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
    result = featurize_reaction(
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )
    signature = result.reaction_signature

    assert signature is not None
    assert signature.signature_id.startswith("RS1:")
    assert signature.exact_signature_key.startswith("L0:")
    assert signature.handle_signature_key.startswith("L1:")
    assert signature.transformation_signature_key.startswith("L2:")
    assert signature.bond_edit_signature_key.startswith("L3:")
    assert signature.environment_signature_key.startswith("L4:")
    assert "signature_features.v1.json" in signature.definition_versions


def test_signature_serializes_with_analysis() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    payload = result.to_dict()

    assert payload["reaction_signature"]["named_family"] is None
    assert payload["reaction_signature"]["order_changes"] == (
        "C-C:DOUBLE>SINGLE",
    )
    assert payload["schema_version"] == "1.3"

