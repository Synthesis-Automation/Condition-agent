from reactive_taxonomy import featurize_reaction, validate_taxonomy


def test_reaction_grammar_taxonomy_validates() -> None:
    assert validate_taxonomy() == []


def test_suzuki_exact_product_reconstruction() -> None:
    reaction = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.valid
    assert result.transformation_class == "c_c_transfer_coupling"
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.reaction_label == "Ar–Br + Ar–B(OH)2 → Ar–Ar"
    assert [change.change_type for change in result.selected_candidate.predicted_bond_changes] == ["broken", "broken", "formed"]
    assert result.product_connection is not None
    assert result.product_connection.endpoint_1.role == "electrophile"
    assert result.product_connection.endpoint_2.role == "transfer_partner"


def test_product_connection_canonicalizes_display_without_swapping_provenance() -> None:
    result = featurize_reaction("ClC=C.OB(O)c1ccccc1>>C=Cc1ccccc1")
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.reaction_label == "Alkenyl–Cl + Ar–B(OH)2 → Ar–Alkenyl"
    connection = result.product_connection
    assert connection.concise_label == "Ar–Alkenyl"
    assert connection.endpoint_1.role == "electrophile"
    assert connection.endpoint_1.context == "Alkenyl"
    assert connection.endpoint_2.role == "transfer_partner"
    assert connection.endpoint_2.context == "Ar"


def test_cn_exact_product_reconstruction() -> None:
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.transformation_class == "sp2_c_n_substitution"
    assert result.reaction_label == "Ar–Br + Ar–NH2 → Ar–NH–Ar"
    assert result.product_connection.connection_type == "C_N"
    assert result.family_environment.family_id == "c_n_coupling"


def test_secondary_amine_product_consumes_nitrogen_hydrogen() -> None:
    result = featurize_reaction("Brc1ccccc1.CNC>>CN(C)c1ccccc1")
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.product_connection.connection_type == "C_N"
    assert "NH" not in result.product_connection.concise_label


def test_amide_nitrogen_context_is_rendered_from_attachment_atom() -> None:
    result = featurize_reaction("Brc1ccccc1.CC(=O)NC>>CC(=O)N(C)c1ccccc1")
    assert result.evidence_quality == "exact_product_reconstruction"
    assert "N(C(O)–R)(R)" in result.product_connection.concise_label


def test_amide_exact_product_reconstruction() -> None:
    reaction = "CC(=O)Cl.CN>>CC(=O)NC"
    result = featurize_reaction(reaction)
    assert result.transformation_class == "amide_formation"
    assert result.reaction_label == "R–C(O)Cl + R–NH2 → R–C(O)–NHR"


def test_wrong_product_is_not_confirmed() -> None:
    result = featurize_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccccc1")
    assert result.valid
    assert result.selected_candidate is None
    assert result.evidence_quality == "reactant_grammar_only"
    assert result.candidates[0].verification == "product_mismatch"
    assert result.reaction_label == "Ar–Br + Ar–NH2 →"
    assert result.reaction_label_status == "reactant_only"


def test_supplied_mapping_yields_exact_bond_differences() -> None:
    reaction = "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"
    result = featurize_reaction(reaction)
    changes = {(change["change_type"], tuple(change["atom_maps"])) for change in result.mapped_bond_changes}
    assert ("broken", (1, 2)) in changes
    assert ("formed", (1, 3)) in changes


def test_three_part_reaction_smiles_preserves_agents() -> None:
    result = featurize_reaction("Brc1ccccc1.Nc1ccccc1>CCN>c1ccc(Nc2ccccc2)cc1")
    assert result.valid
    assert len(result.agents) == 1
    assert result.transformation_class == "sp2_c_n_substitution"


def test_additional_v1_grammars_reconstruct_products() -> None:
    cases = {
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1": "sp2_c_o_substitution",
        "Brc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1": "sp2_c_s_substitution",
        "Brc1ccccc1.C#Cc1ccccc1>>c1ccc(C#Cc2ccccc2)cc1": "sonogashira_coupling",
        "Brc1ccccc1.Cl[Zn]c1ccccc1>>c1ccc(-c2ccccc2)cc1": "other_metal_transfer_coupling",
        "OB(O)c1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1": "chan_lam_heteroatom_coupling",
        "CC(=O)Cl.CCO>>CCOC(C)=O": "ester_formation",
        "CC(=O)Cl.CCS>>CCSC(C)=O": "thioester_formation",
        "CS(=O)(=O)Cl.CN>>CS(=O)(=O)NC": "sulfonamide_formation",
    }
    for reaction, grammar_id in cases.items():
        result = featurize_reaction(reaction)
        assert result.evidence_quality == "exact_product_reconstruction", reaction
        assert result.selected_candidate.grammar_id == grammar_id


def test_oxygen_and_sulfur_products_have_typed_connections_and_environments() -> None:
    cases = {
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1": ("C_O", "c_o_coupling", "Ar–O–Ar"),
        "Brc1ccccc1.CCS>>CCSc1ccccc1": ("C_S", "c_s_coupling", "Ar–S–R"),
    }
    for reaction, (connection_type, family_id, label) in cases.items():
        result = featurize_reaction(reaction)
        assert result.evidence_quality == "exact_product_reconstruction"
        assert result.product_connection.connection_type == connection_type
        assert result.product_connection.concise_label == label
        assert result.family_environment.family_id == family_id


def test_composite_triflate_handle_is_removed_as_complete_fragment() -> None:
    reaction = "Fc1ccc(OS(=O)(=O)C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>Fc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.evidence_quality == "exact_product_reconstruction"
    assert result.selected_candidate.grammar_id == "boron_transfer_coupling"
    assert result.selected_candidate.predicted_product_smiles == "Fc1ccc(-c2ccccc2)cc1"
