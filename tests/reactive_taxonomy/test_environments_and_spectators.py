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


def test_aldehyde_definition_includes_formaldehyde_but_excludes_formamides() -> None:
    assert [group.group_id for group in featurize_molecule("C=O").functional_groups] == ["aldehyde"]
    assert [group.group_id for group in featurize_molecule("CC=O").functional_groups] == ["aldehyde"]
    assert [group.group_id for group in featurize_molecule("CN(C)C=O").functional_groups] == ["amide"]
    assert [group.group_id for group in featurize_molecule("NC=O").functional_groups] == ["amide"]


def test_specific_carbonyl_groups_own_generic_overlaps() -> None:
    examples = {
        "CC(=O)OC(=O)C": {"acid_anhydride"},
        "COC(=O)NC": {"carbamate"},
        "O=C1C=CC(=O)N1c1ccccc1": {"imide", "alkene"},
    }
    for smiles, expected in examples.items():
        ids = {group.group_id for group in featurize_molecule(smiles).functional_groups}
        assert ids == expected
    anhydride_sites = [
        site for site in featurize_molecule("CC(=O)OC(=O)C").sites
        if site.site_type == "electrophilic_center"
    ]
    assert len(anhydride_sites) == 2
    assert {site.availability for site in anhydride_sites} == {"activated"}
    assert {site.details["acyl_subtype"] for site in anhydride_sites} == {"acid_anhydride"}
    protected_aniline = featurize_molecule("Brc1ccc(NC(=O)OC(C)(C)C)cc1")
    assert not [
        site for site in protected_aniline.sites
        if site.site_type == "electrophilic_center"
    ]


def test_common_synthesis_functional_group_pack() -> None:
    examples = {
        "c1ccccc1C1CO1": "epoxide",
        "C(=Nc1ccccc1)c1ccccc1": "imine",
        "[N-]=[N+]=Nc1ccccc1": "azide",
        "O=C=Nc1ccccc1": "isocyanate",
        "S=C=Nc1ccccc1": "isothiocyanate",
        "CS(C)=O": "sulfoxide",
        "CSSC": "disulfide",
        "CCOP(=O)(CC(=O)OCC)OCC": "phosphonate",
    }
    for smiles, expected in examples.items():
        assert expected in {group.group_id for group in featurize_molecule(smiles).functional_groups}


def test_site_ids_use_the_reactive_locus() -> None:
    aldehyde = next(
        site for site in featurize_molecule("CC=O").sites
        if site.site_type == "electrophilic_center"
    )
    assert aldehyde.details["center_atom_index"] == 1
    assert aldehyde.site_id.startswith("mol0:atom1:")

    alkene = next(
        site for site in featurize_molecule("CC=C").sites
        if site.site_type == "unsaturated_bond"
    )
    assert alkene.bond_indices == [1]
    assert alkene.site_id.startswith("mol0:bond1:")


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


def test_nitrogen_substitution_class_does_not_count_sulfur_as_carbon() -> None:
    environment = featurize_molecule("CS(=O)(=O)NC").site_environments[0]
    assert environment.steric["center_substitution_class"] == "secondary"
