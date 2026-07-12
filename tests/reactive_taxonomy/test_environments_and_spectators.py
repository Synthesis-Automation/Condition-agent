from reactive_taxonomy import featurize_molecule, featurize_reaction, validate_taxonomy


def test_functional_groups_and_site_environments_are_public() -> None:
    result = featurize_molecule("Brc1ccccc1C#N")
    assert result.valid
    assert {group.group_id for group in result.functional_groups} >= {"aryl_halide", "nitrile"}
    assert result.site_environments
    leaving_environment = next(
        environment for environment in result.site_environments
        if environment.site_id == next(site.site_id for site in result.sites if site.site_type == "leaving_group")
    )
    assert leaving_environment.steric["method"] == "graph_local_v1"
    assert any(group["group_id"] == "nitrile" for group in leaving_environment.nearby_groups)


def test_reaction_spectators_exclude_consumed_handle_and_keep_nitrile() -> None:
    reaction = "Brc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.evidence_quality == "exact_product_reconstruction"
    ids = {group.group_id for group in result.spectator_groups}
    assert "nitrile" in ids
    assert "aryl_halide" not in ids
    nitrile = next(group for group in result.spectator_groups if group.group_id == "nitrile")
    assert nitrile.graph_distance is not None
    assert nitrile.unchanged_evidence == "exact_product_reconstruction_event_exclusion"


def test_functional_group_registry_is_validated() -> None:
    assert validate_taxonomy() == []


def test_functional_group_ownership_suppresses_generic_sulfonamide_amine() -> None:
    result = featurize_molecule("CS(=O)(=O)NC")
    ids = [group.group_id for group in result.functional_groups]
    assert "sulfonamide" in ids
    assert "secondary_amine" not in ids


def test_suzuki_family_environment_separates_partner_roles() -> None:
    reaction = "Brc1ccccc1.OB(O)c1ccccn1>>c1ccc(-c2ccccn2)cc1"
    result = featurize_reaction(reaction)
    assert result.named_family == "suzuki_miyaura"
    assert result.family_environment is not None
    partners = {partner.role: partner for partner in result.family_environment.partners}
    assert partners["electrophile"].handle_token == "Br"
    assert partners["electrophile"].anchor_context == "Ar"
    assert partners["transfer_partner"].handle_token == "B(OH)2"
    assert partners["transfer_partner"].anchor_context == "HeteroAr"
    assert "protodeboronation_risk" in partners["transfer_partner"].flags


def test_suzuki_environment_reports_competing_halide() -> None:
    reaction = "Brc1ccc(Br)cc1.OB(O)c1ccccc1>>Brc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.family_environment is not None
    electrophile = next(
        partner for partner in result.family_environment.partners
        if partner.role == "electrophile"
    )
    assert electrophile.competing_site_labels == ("Ar–Br",)
    assert "competing_reactive_handle" in electrophile.flags


def test_suzuki_site_descriptors_distinguish_ortho_bulk_and_electronics() -> None:
    result = featurize_molecule("Brc1c(C)cccc1C#N")
    leaving_site = next(site for site in result.sites if site.site_type == "leaving_group")
    environment = next(item for item in result.site_environments if item.site_id == leaving_site.site_id)
    assert environment.steric["ortho_substituent_count"] >= 1
    assert environment.electronic["class"] == "electron_poor"


def test_alkyl_attachment_sterics_are_distinct_from_amine_substitution() -> None:
    examples = {
        "CN": "methyl",
        "CCN": "primary",
        "CC(C)N": "secondary",
        "CC(C)(C)N": "tertiary",
    }
    for smiles, expected in examples.items():
        result = featurize_molecule(smiles)
        site = next(item for item in result.sites if item.site_type == "pronucleophile_XH")
        environment = next(item for item in result.site_environments if item.site_id == site.site_id)
        assert environment.steric["center_substitution_class"] == "primary"
        alkyl = next(group for group in environment.steric["attached_groups"] if group["context"] == "Alkyl")
        assert alkyl["attachment_carbon_class"] == expected

    tert_butyl = featurize_molecule("CC(C)(C)N").site_environments[0].steric["attached_groups"][0]
    assert tert_butyl["alpha_branched"] is True
