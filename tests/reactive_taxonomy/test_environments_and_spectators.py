from reactive_taxonomy import analyze_molecule, featurize_reaction, validate_taxonomy


def test_functional_groups_and_site_environments_are_public() -> None:
    result = analyze_molecule("Brc1ccccc1C#N")
    assert result.valid
    assert {group.motif_id for group in result.motifs} >= {"aryl_halide", "nitrile"}
    assert result.reactive_site_environments
    leaving_environment = next(
        environment for environment in result.reactive_site_environments
        if environment.hypothesis_id == next(site.hypothesis_id for site in result.reactive_site_hypotheses if site.site_type == "leaving_group")
    )
    assert leaving_environment.reactivity_profile is not None
    assert (
        leaving_environment.reactivity_profile.steric.evidence.method
        == "aromatic_steric_graph_v1"
    )
    assert any(group["motif_id"] == "nitrile" for group in leaving_environment.nearby_motifs)


def test_reaction_spectators_exclude_consumed_handle_and_keep_nitrile() -> None:
    reaction = "Brc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.evidence_quality == "global_atom_correspondence"
    ids = {group.group_id for group in result.spectator_groups}
    assert "nitrile" in ids
    assert "aryl_halide" not in ids
    nitrile = next(group for group in result.spectator_groups if group.group_id == "nitrile")
    assert nitrile.graph_distance is not None
    assert nitrile.unchanged_evidence == "global_atom_correspondence_edit_exclusion"


def test_functional_group_registry_is_validated() -> None:
    assert validate_taxonomy() == []


def test_heteroaromatic_atoms_are_localized_functional_groups() -> None:
    examples = {
        "c1ccncc1": "heteroaryl_nitrogen",
        "c1cc[nH]c1": "pyrrolic_nitrogen",
        "c1ccoc1": "heteroaryl_oxygen",
        "c1ccsc1": "heteroaryl_sulfur",
        "[O-][n+]1ccccc1": "heteroaryl_n_oxide",
        "c1ccc[nH+]c1": "cationic_heteroaryl_nitrogen",
        "C[n+]1ccccc1": "cationic_heteroaryl_nitrogen",
    }
    for smiles, expected in examples.items():
        groups = analyze_molecule(smiles).motifs
        matching = [group for group in groups if group.motif_id == expected]
        assert len(matching) == 1
        assert len(matching[0].atom_indices) == (
            2 if expected == "heteroaryl_n_oxide" else 1
        )

    n_oxide_ids = {
        group.motif_id
        for group in analyze_molecule("[O-][n+]1ccccc1").motifs
    }
    assert "heteroaryl_nitrogen" not in n_oxide_ids


def test_pyridine_nitrogen_is_spectator_but_reactive_ring_is_context() -> None:
    reaction = (
        "Brc1ncccc1.OB(O)c1ccccc1"
        ">>c1ccc(-c2ncccc2)cc1"
    )
    result = featurize_reaction(reaction)
    assert result.interpretation is not None
    assert result.interpretation.family_environment is None
    assert result.interpretation.primary_pattern_id == (
        "organoboron_c_c_coupling_like"
    )
    spectator_ids = {
    group.group_id for group in result.spectator_groups
    }

    assert "heteroaryl_nitrogen" in spectator_ids
    assert "aryl_halide" not in spectator_ids
    pyridine_nitrogen = next(
        group
        for group in result.spectator_groups
        if group.group_id == "heteroaryl_nitrogen"
    )
    assert pyridine_nitrogen.graph_distance == 1
    assert pyridine_nitrogen.unchanged_evidence == (
        "global_atom_correspondence_edit_exclusion"
    )


def test_functional_group_ownership_suppresses_generic_sulfonamide_amine() -> None:
    result = analyze_molecule("CS(=O)(=O)NC")
    ids = [group.motif_id for group in result.motifs]
    assert "sulfonamide" in ids
    assert "secondary_amine" not in ids


def test_nitro_nitrogen_is_not_classified_as_a_tertiary_amine() -> None:
    result = analyze_molecule("O=[N+]([O-])c1ccccc1")
    ids = {group.motif_id for group in result.motifs}
    assert "nitro" in ids
    assert "tertiary_amine" not in ids


def test_nitro_spectator_is_not_reported_as_a_tertiary_amine() -> None:
    reaction = (
        "Sc1ccccc1.Sc1ccccc1.O=[N+]([O-])c1c(F)cccc1F"
        ">>O=[N+]([O-])c1c(Sc2ccccc2)cccc1Sc1ccccc1"
    )
    result = featurize_reaction(reaction)
    ids = [group.group_id for group in result.spectator_groups]
    assert ids.count("nitro") == 1
    assert "tertiary_amine" not in ids


def test_aldehyde_definition_includes_formaldehyde_but_excludes_formamides() -> None:
    assert [group.motif_id for group in analyze_molecule("C=O").motifs] == ["aldehyde"]
    assert [group.motif_id for group in analyze_molecule("CC=O").motifs] == ["aldehyde"]
    assert [group.motif_id for group in analyze_molecule("CN(C)C=O").motifs] == ["amide"]
    assert [group.motif_id for group in analyze_molecule("NC=O").motifs] == ["amide"]


def test_specific_carbonyl_groups_own_generic_overlaps() -> None:
    examples = {
        "CC(=O)OC(=O)C": {"acid_anhydride"},
        "COC(=O)NC": {"carbamate"},
        "O=C1C=CC(=O)N1c1ccccc1": {"imide", "alkene"},
    }
    for smiles, expected in examples.items():
        ids = {group.motif_id for group in analyze_molecule(smiles).motifs}
        assert ids == expected
    anhydride_sites = [
        site for site in analyze_molecule("CC(=O)OC(=O)C").reactive_site_hypotheses
        if site.site_type == "electrophilic_center"
    ]
    assert len(anhydride_sites) == 2
    assert {site.availability for site in anhydride_sites} == {"activated"}
    assert {site.details["acyl_subtype"] for site in anhydride_sites} == {"acid_anhydride"}
    protected_aniline = analyze_molecule("Brc1ccc(NC(=O)OC(C)(C)C)cc1")
    assert not [
        site for site in protected_aniline.reactive_site_hypotheses
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
        assert expected in {group.motif_id for group in analyze_molecule(smiles).motifs}


def test_site_ids_use_the_reactive_locus() -> None:
    aldehyde = next(
        site for site in analyze_molecule("CC=O").reactive_site_hypotheses
        if site.site_type == "electrophilic_center"
    )
    assert aldehyde.details["center_atom_index"] == 1
    assert aldehyde.hypothesis_id.startswith("mol0:atom1:")

    alkene = next(
        site for site in analyze_molecule("CC=C").reactive_site_hypotheses
        if site.site_type == "unsaturated_bond"
    )
    assert alkene.bond_indices == (1,)
    assert alkene.hypothesis_id.startswith("mol0:bond1:")


def test_suzuki_pattern_is_optional_and_does_not_assign_partner_roles() -> None:
    reaction = "Brc1ccccc1.OB(O)c1ccccn1>>c1ccc(-c2ccccn2)cc1"
    result = featurize_reaction(reaction)
    assert result.named_family == "suzuki_miyaura"
    assert result.family_environment is None
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == (
        "organoboron_c_c_coupling_like"
    )
    assert result.interpretation.partners == ()


def test_molecular_annotations_remain_available_beside_suzuki_pattern() -> None:
    reaction = "Brc1ccc(Br)cc1.OB(O)c1ccccc1>>Brc1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(reaction)
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == (
        "organoboron_c_c_coupling_like"
    )
    leaving_groups = [
        site
        for site in result.reactants[0].molecule_analysis.reactive_site_hypotheses
        if site.site_type == "leaving_group"
    ]
    assert len(leaving_groups) == 2


def test_suzuki_site_descriptors_distinguish_ortho_bulk_and_electronics() -> None:
    result = analyze_molecule("Brc1c(C)cccc1C#N")
    leaving_site = next(site for site in result.reactive_site_hypotheses if site.site_type == "leaving_group")
    environment = next(item for item in result.reactive_site_environments if item.hypothesis_id == leaving_site.hypothesis_id)
    profile = environment.reactivity_profile
    assert profile is not None
    assert profile.context.ortho_occupancy_count == 2
    assert profile.context.ortho_burden_class == "high"
    assert profile.electronic.activation_class == "slightly_poor"


def test_alkyl_attachment_sterics_are_distinct_from_amine_substitution() -> None:
    examples = {
        "CN": "methyl",
        "CCN": "primary",
        "CC(C)N": "secondary",
        "CC(C)(C)N": "tertiary",
    }
    for smiles, expected in examples.items():
        result = analyze_molecule(smiles)
        site = next(item for item in result.reactive_site_hypotheses if item.site_type == "pronucleophile_XH")
        environment = next(item for item in result.reactive_site_environments if item.hypothesis_id == site.hypothesis_id)
        profile = environment.reactivity_profile
        assert profile is not None
        assert profile.reactive_center.substitution_class == "primary"
        alkyl = next(
            group
            for group in profile.context.attached_groups
            if group.context == "Alkyl"
        )
        assert alkyl.attachment_carbon_class == expected

    tert_butyl_profile = (
        analyze_molecule("CC(C)(C)N")
        .reactive_site_environments[0]
        .reactivity_profile
    )
    assert tert_butyl_profile is not None
    assert tert_butyl_profile.context.attached_groups[0].alpha_branched is True


def test_nitrogen_substitution_class_does_not_count_sulfur_as_carbon() -> None:
    environment = analyze_molecule("CS(=O)(=O)NC").reactive_site_environments[0]
    assert environment.reactivity_profile is not None
    assert (
        environment.reactivity_profile.reactive_center.substitution_class
        == "secondary"
    )
