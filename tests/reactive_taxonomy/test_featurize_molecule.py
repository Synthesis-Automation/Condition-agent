from reactive_taxonomy import MatchIndex, ReactiveSiteCandidate, available_styles, analyze_molecule, load_handle_patterns, validate_taxonomy
from reactive_taxonomy.context import load_context_taxonomy
from reactive_taxonomy.sites import pronucleophiles


def signatures(smiles: str) -> set[str]:
    result = analyze_molecule(smiles)
    assert result.valid, result.error
    return {site.canonical_signature for site in result.reactive_site_hypotheses}


def test_taxonomy_bundle_validates() -> None:
    assert validate_taxonomy() == []


def test_handle_smarts_are_independent_and_mapped() -> None:
    patterns = load_handle_patterns()
    assert len(patterns) >= 12
    assert {pattern["site_type"] for pattern in patterns} == {
        "leaving_group", "pronucleophile_XH", "transfer_group",
        "electrophilic_center", "aromatic_CH", "unsaturated_bond",
        "dipolar_group", "heteroatom_bond", "nucleophile_anion",
        "addition_donor", "eliminable_pair",
    }
    assert all(pattern.get("atom_roles") for pattern in patterns)


def test_site_reports_pattern_provenance() -> None:
    site = analyze_molecule("Brc1ccccc1").reactive_site_hypotheses[0]
    assert site.details["matched_pattern"] == "terminal_carbon_halogen"
    assert site.details["alternative_patterns"] == []


def test_context_taxonomy_declares_classification_methods() -> None:
    payload = load_context_taxonomy()
    records = {record["id"]: record for record in payload["contexts"]}
    assert records["HetAr"]["classification_method"] == "aromatic_ring_system"
    assert records["Alkyl"]["classification_method"] == "atom_property"
    assert records["C(O)NR"]["classification_method"] == "mapped_smarts"
    assert records["C(O)NR"]["atom_roles"]["context_anchor"] == 1


def test_leaving_groups() -> None:
    assert "LG|Ar|Br" in signatures("Brc1ccccc1")
    assert "LG|HetAr|Cl" in signatures("Clc1ncccc1")
    assert "LG|Alkyl|OMs" in signatures("CCOS(=O)(=O)C")


def test_pronucleophiles() -> None:
    cases = {
        "CN": "XH|N|H2|Alkyl",
        "Nc1ccccc1": "XH|N|H2|Ar",
        "CCNCC": "XH|N|H1|Alkyl,Alkyl",
        "CNc1ccccc1": "XH|N|H1|Ar,Alkyl",
        "CC(=O)NC": "XH|N|H1|C(O)R,Alkyl",
        "CS(=O)(=O)NC": "XH|N|H1|SO2R,Alkyl",
        "Oc1ccccc1": "XH|O|H1|Ar",
        "CCS": "XH|S|H1|Alkyl",
    }
    for smiles, expected in cases.items():
        assert expected in signatures(smiles)


def test_explicit_sulfur_anions_are_separate_nucleophile_sites() -> None:
    thioacetate = analyze_molecule("CC(=O)[S-]")
    methylthiolate = analyze_molecule("C[S-]")

    assert {
        (site.site_type, site.canonical_signature, site.chemist_label)
        for site in thioacetate.reactive_site_hypotheses
    } == {
        ("nucleophile_anion", "NU-|S|-1|C(O)R", "R–C(O)–S⁻")
    }
    assert {
        (site.site_type, site.canonical_signature, site.chemist_label)
        for site in methylthiolate.reactive_site_hypotheses
    } == {
        ("nucleophile_anion", "NU-|S|-1|Alkyl", "Alk–S⁻")
    }


def test_hydrazine_has_two_sites() -> None:
    nh_sites = [s for s in analyze_molecule("CNN").reactive_site_hypotheses if s.site_type == "pronucleophile_XH" and s.details["center_element"] == "N"]
    assert len(nh_sites) == 2


def test_transfer_groups() -> None:
    assert "TM|Ar|B(OH)2" in signatures("OB(O)c1ccccc1")
    assert "TM|Ar|BPin" in signatures("CC1(C)OB(c2ccccc2)OC1(C)C")


def test_bpin_oxygens_are_not_reported_as_generic_ethers() -> None:
    result = analyze_molecule(
        "Cc1c(B2OC(C)(C)C(C)(C)O2)c2c(cnn2C3OCCCC3)cc1Cl"
    )

    assert any(
        site.canonical_signature == "TM|HetAr|BPin"
        for site in result.reactive_site_hypotheses
    )
    ether_matches = [
        group for group in result.motifs if group.motif_id == "ether"
    ]
    assert len(ether_matches) == 1


def test_generic_ether_still_requires_two_carbon_neighbors() -> None:
    assert [
        group.motif_id for group in analyze_molecule("COC").motifs
    ] == ["ether"]


def test_electrophilic_centers() -> None:
    assert "EC|Acyl|Alkyl|OH|latent" in signatures("CC(=O)O")
    assert "EC|Acyl|Alkyl|Cl|activated" in signatures("CC(=O)Cl")
    assert "EC|Sulfonyl|Ar|Cl|activated" in signatures("O=S(=O)(Cl)c1ccccc1")


def test_invalid_input_and_family_filter() -> None:
    assert analyze_molecule("not smiles").structure.error == "INVALID_SMILES"
    result = analyze_molecule("Brc1ccccc1N", site_types={"leaving_group"})
    assert result.valid
    assert {s.site_type for s in result.reactive_site_hypotheses} == {"leaving_group"}


def test_components_are_preserved() -> None:
    result = analyze_molecule("Nc1ccccc1.[K+]")
    assert result.valid
    assert len(result.components) == 2


def test_composite_handles_do_not_emit_internal_sites() -> None:
    result = analyze_molecule("c1ccc(B(O)O)cc1")
    non_aromatic_ch = [site for site in result.reactive_site_hypotheses if site.site_type != "aromatic_CH"]
    assert [site.canonical_signature for site in non_aromatic_ch] == ["TM|Ar|B(OH)2"]


def test_retained_fluorines_are_not_leaving_groups() -> None:
    sigs = {
        site.canonical_signature for site in analyze_molecule("Brc1ccc(C(F)(F)F)cc1").reactive_site_hypotheses
        if site.site_type == "leaving_group"
    }
    assert sigs == {"LG|Ar|Br"}


def test_sn_and_si_emit_only_transferable_carbon_site() -> None:
    tin_sites = [
        site.canonical_signature for site in analyze_molecule("c1ccc([Sn](C)(C)C)cc1").reactive_site_hypotheses
        if site.site_type == "transfer_group"
    ]
    assert tin_sites == ["TM|Ar|SnR3"]
    silicon_sites = {
        site.canonical_signature for site in analyze_molecule("C#C[Si](C)(C)C").reactive_site_hypotheses
        if site.site_type in {"pronucleophile_XH", "transfer_group"}
    }
    assert silicon_sites == {"XH|Csp|H1|Alkynyl", "TM|Alkynyl|SiR3"}


def test_silyl_ether_is_not_a_transfer_group() -> None:
    result = analyze_molecule("Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1")
    assert not [site for site in result.reactive_site_hypotheses if site.site_type == "transfer_group"]
    assert {site.canonical_signature for site in result.reactive_site_hypotheses if site.site_type == "leaving_group"} == {"LG|Ar|Br"}


def test_bridging_metal_halogen_is_not_a_leaving_group() -> None:
    assert "LG|Ar|Br" not in signatures("[Mg]Brc1ccccc1")


def test_acyl_halide_owns_halogen_site() -> None:
    result = analyze_molecule("CC(=O)Cl")
    assert [site.canonical_signature for site in result.reactive_site_hypotheses] == [
        "EC|Acyl|Alkyl|Cl|activated"
    ]


def test_ammonia_is_supported_nh_pronucleophile() -> None:
    result = analyze_molecule("N")
    assert [site.canonical_signature for site in result.reactive_site_hypotheses] == ["XH|N|H3|"]
    assert result.reactive_site_hypotheses[0].chemist_label == "NH₃"


def test_expanded_condensation_and_nitrogen_classes() -> None:
    expected = {
        "O=C(O)c1ccccc1": "EC|Acyl|Ar|OH|latent",
        "CS(=O)(=O)Cl": "EC|Sulfonyl|Alkyl|Cl|activated",
        "NC(=O)N": "XH|N|H2|C(O)NR",
        "CS(=O)(=O)NC": "XH|N|H1|SO2R,Alkyl",
        "NNC(=O)c1ccccc1": "XH|N|H1|C(O)R,N",
    }
    for smiles, signature in expected.items():
        assert signature in signatures(smiles)


def test_aromatic_nh_is_one_ring_context() -> None:
    for smiles in ("c1cc[nH]c1", "c1ccc2[nH]ccc2c1"):
        sites = [
            site for site in analyze_molecule(smiles).reactive_site_hypotheses
            if site.site_type == "pronucleophile_XH"
        ]
        assert len(sites) == 1
        assert sites[0].canonical_signature == "XH|N_aromatic|H1|HetAr"
        assert sites[0].chemist_label == "HetAr–[nH]"
        assert sites[0].details["center_element"] == "N"
        assert sites[0].details["derived_family"] == "aromatic_nh"


def test_bromopyrrole_keeps_both_distinct_sites() -> None:
    legacy_sites = {
        site.canonical_signature for site in analyze_molecule("Brc1cc[nH]c1").reactive_site_hypotheses
        if site.site_type in {"leaving_group", "pronucleophile_XH"}
    }
    assert legacy_sites == {
        "LG|HetAr|Br",
        "XH|N_aromatic|H1|HetAr",
    }


def test_rendering_styles_preserve_signature() -> None:
    unicode_site = analyze_molecule("Nc1ccccc1").reactive_site_hypotheses[0]
    ascii_site = analyze_molecule("Nc1ccccc1", label_style="ascii").reactive_site_hypotheses[0]
    hte_site = analyze_molecule("Nc1ccccc1", label_style="hte_legacy").reactive_site_hypotheses[0]
    assert available_styles() == ("unicode", "ascii", "hte_legacy")
    assert unicode_site.chemist_label == "Ar–NH₂"
    assert ascii_site.chemist_label == "Ar-NH2"
    assert hte_site.chemist_label == "Ar-NH2"
    assert len({unicode_site.canonical_signature, ascii_site.canonical_signature, hte_site.canonical_signature}) == 1
    amide_nh = next(
        site
        for site in analyze_molecule("CC(=O)NC").reactive_site_hypotheses
        if site.details.get("center_token") == "N"
    )
    sulfonamide_nh = next(
        site
        for site in analyze_molecule("CS(=O)(=O)NC").reactive_site_hypotheses
        if site.details.get("center_token") == "N"
    )
    assert amide_nh.chemist_label == "R–C(O)–NHR"
    assert sulfonamide_nh.chemist_label == "R–S(O)₂–NHR"


def test_availability_distinguishes_chemical_site_state() -> None:
    assert analyze_molecule("CCN").reactive_site_hypotheses[0].availability == "free"
    assert next(
        site
        for site in analyze_molecule("CC(=O)NC").reactive_site_hypotheses
        if site.details.get("center_token") == "N"
    ).availability == "deactivated"
    assert next(site for site in analyze_molecule("CC(=O)O").reactive_site_hypotheses if site.site_type == "electrophilic_center").availability == "latent"
    assert analyze_molecule("CC(=O)Cl").reactive_site_hypotheses[0].availability == "activated"
    assert analyze_molecule("OB(O)c1ccccc1").reactive_site_hypotheses[0].availability == "transferable"


def test_rich_context_record_is_exposed() -> None:
    site = analyze_molecule("Nc1ccccn1").reactive_site_hypotheses[0]
    context = site.context_features["contexts"][0]
    assert context["token"] == "HetAr"
    assert context["classification_method"] == "aromatic_ring_system"
    assert context["subtype"] == "pyridine_like"
    assert context["features"]["heteroatoms"] == ["N"]
    assert context["features"]["heteroatom_counts"] == {"N": 1}
    assert context["features"]["heteroatom_details"] == [{
        "atom_index": 6,
        "element": "N",
        "formal_charge": 0,
        "hydrogen_count": 0,
        "aromatic_type": "pyridine_like",
        "distance_from_attachment": 1,
    }]
    assert context["features"]["ring_system_key"] == (
        "HetAr|rings=6|hetero=N:pyridine_like:d1|fused=0"
    )
    assert len(context["fragment_atom_indices"]) == 6


def test_heteroaromatic_context_preserves_broad_token_and_nitrogen_position() -> None:
    examples = {
        "Brc1ncccc1": 1,
        "Brc1cccnc1": 2,
        "Brc1ccncc1": 3,
    }
    for smiles, expected_distance in examples.items():
        site = next(
            candidate
            for candidate in analyze_molecule(smiles).reactive_site_hypotheses
            if candidate.site_type == "leaving_group"
        )
        context = site.context_features["contexts"][0]
        assert site.canonical_signature == "LG|HetAr|Br"
        assert context["token"] == "HetAr"
        assert context["subtype"] == "pyridine_like"
        assert (
            context["features"]["heteroatom_details"][0][
                "distance_from_attachment"
            ]
            == expected_distance
        )


def test_graph_defined_heteroaromatic_subtypes_cover_common_ring_classes() -> None:
    examples = {
        "Brc1ccoc1": ("furan_like", "O"),
        "Brc1ccsc1": ("thiophene_like", "S"),
        "Brc1cc[nH]c1": ("pyrrole_like", "N"),
        "Brc1ccc2ncccc2c1": ("fused_pyridine_like", "N"),
    }
    for smiles, (expected_subtype, expected_element) in examples.items():
        site = next(
            candidate
            for candidate in analyze_molecule(smiles).reactive_site_hypotheses
            if candidate.site_type == "leaving_group"
        )
        context = site.context_features["contexts"][0]
        assert context["token"] == "HetAr"
        assert context["subtype"] == expected_subtype
        assert context["features"]["heteroatom_counts"] == {
            expected_element: 1
        }


def test_charged_heteroaromatic_nitrogen_states_are_distinct_contexts() -> None:
    examples = {
        "Brc1cccc[n+]1[O-]": "pyridine_n_oxide_like",
        "Brc1cccc[nH+]1": "pyridinium_like",
        "Brc1cccc[n+]1C": "pyridinium_like",
    }
    for smiles, expected_subtype in examples.items():
        site = next(
            candidate
            for candidate in analyze_molecule(smiles).reactive_site_hypotheses
            if candidate.site_type == "leaving_group"
        )
        context = site.context_features["contexts"][0]
        assert context["token"] == "HetAr"
        assert context["subtype"] == expected_subtype
        assert context["features"]["heteroatom_details"][0][
            "formal_charge"
        ] == 1


def test_alkyl_leaving_groups_preserve_benzylic_allylic_and_propargylic_subtypes() -> None:
    examples = {
        "ClCc1ccccc1": ("benzylic", "Bn–Cl"),
        "ClCC=C": ("allylic", "Allyl–Cl"),
        "ClCC#C": ("propargylic", "Propargyl–Cl"),
        "CCCCl": ("simple_alkyl", "Alk–Cl"),
    }
    for smiles, (subtype, label) in examples.items():
        site = next(
            item
            for item in analyze_molecule(smiles).reactive_site_hypotheses
            if item.site_type == "leaving_group"
        )
        context = site.context_features["contexts"][0]
        assert site.details["anchor_context"] == "Alkyl"
        assert site.details["anchor_subtype"] == subtype
        assert context["subtype"] == subtype
        assert site.chemist_label == label


def test_benzyl_chloride_and_aryl_bromide_are_separate_sites() -> None:
    result = analyze_molecule("ClCc1ccccc1Br")
    assert {
        site.chemist_label for site in result.reactive_site_hypotheses if site.site_type == "leaving_group"
    } == {"Bn–Cl", "Ar–Br"}


def test_aldehydes_and_ketones_are_carbonyl_addition_centers() -> None:
    aldehyde = next(site for site in analyze_molecule("CC=O").reactive_site_hypotheses if site.site_type == "electrophilic_center")
    assert aldehyde.canonical_signature == "EC|Carbonyl|aldehyde|Alkyl|addition"
    assert aldehyde.details["reaction_mode"] == "addition"
    assert aldehyde.details["atom_roles"] == {"center": [1], "heteroatom": [2], "substituents": [0]}
    assert aldehyde.chemist_label == "R–CH=O"

    ketone = next(site for site in analyze_molecule("CC(=O)C").reactive_site_hypotheses if site.site_type == "electrophilic_center")
    assert ketone.canonical_signature == "EC|Carbonyl|ketone|Alkyl,Alkyl|addition"
    assert ketone.chemist_label == "R²C=O"


def test_acyl_substitution_and_carbonyl_addition_are_distinct() -> None:
    acid_sites = [site for site in analyze_molecule("CC(=O)O").reactive_site_hypotheses if site.site_type == "electrophilic_center"]
    assert [site.details["center_family"] for site in acid_sites] == ["Acyl"]
    assert acid_sites[0].details["reaction_mode"] == "substitution"
    assert not [site for site in analyze_molecule("CC(=O)N").reactive_site_hypotheses if site.site_type == "electrophilic_center"]


def test_aromatic_ch_sites_are_atom_localized_and_ring_classified() -> None:
    benzene_sites = [site for site in analyze_molecule("c1ccccc1").reactive_site_hypotheses if site.site_type == "aromatic_CH"]
    assert len(benzene_sites) == 6
    assert {site.canonical_signature for site in benzene_sites} == {"CH|ArH"}
    assert {site.topology for site in benzene_sites} == {"atom"}
    assert {site.chemist_label for site in benzene_sites} == {"Ar–H"}

    pyridine_sites = [site for site in analyze_molecule("c1ccncc1").reactive_site_hypotheses if site.site_type == "aromatic_CH"]
    assert len(pyridine_sites) == 5
    assert {site.details["handle_token"] for site in pyridine_sites} == {"HetArH"}
    assert {site.chemist_label for site in pyridine_sites} == {"HetAr–H"}


def test_unsaturated_carbon_bonds_are_bond_localized() -> None:
    alkene = next(site for site in analyze_molecule("CC=C").reactive_site_hypotheses if site.site_type == "unsaturated_bond")
    assert alkene.topology == "bond"
    assert alkene.canonical_signature == "PI|Alkene"
    assert alkene.details["bond_order"] == 2
    assert set(alkene.details["atom_roles"]) == {"endpoint_a", "endpoint_b"}
    assert alkene.chemist_label == "H₂C=CHR¹"
    assert alkene.details["endpoint_h_counts"] == [2, 1]
    assert alkene.details["endpoint_substituent_counts"] == [0, 1]

    alkyne_sites = analyze_molecule("CC#C").reactive_site_hypotheses
    assert {site.site_type for site in alkyne_sites} == {"unsaturated_bond", "pronucleophile_XH"}
    assert next(
        site for site in alkyne_sites if site.site_type == "unsaturated_bond"
    ).chemist_label == "R¹–C≡C–H"


def test_alkene_labels_expose_all_hydrogen_and_substituent_positions() -> None:
    examples = {
        "C=C": "H₂C=CH₂",
        "CC=C": "H₂C=CHR¹",
        "CC=CC": "R¹HC=CHR²",
        "C=C(C)c1ccccc1": "H₂C=CR¹R²",
        "CC(C)=C(C)C": "R¹R²C=CR³R⁴",
    }

    for smiles, expected in examples.items():
        site = next(
            site
            for site in analyze_molecule(smiles).reactive_site_hypotheses
            if site.site_type == "unsaturated_bond"
        )
        assert site.chemist_label == expected
        assert sum(site.details["endpoint_h_counts"]) + sum(
            site.details["endpoint_substituent_counts"]
        ) == 4


def test_alkene_labels_retain_defined_e_z_stereochemistry() -> None:
    trans = next(
        site
        for site in analyze_molecule("c1ccc(/C=C/c2ccccc2)cc1").reactive_site_hypotheses
        if site.site_type == "unsaturated_bond"
    )
    cis = next(
        site
        for site in analyze_molecule("c1ccc(/C=C\\c2ccccc2)cc1").reactive_site_hypotheses
        if site.site_type == "unsaturated_bond"
    )

    assert trans.chemist_label == "R¹HC=CHR² (E)"
    assert cis.chemist_label == "R¹HC=CHR² (Z)"


def test_alkyne_labels_distinguish_acetylene_terminal_and_internal() -> None:
    examples = {
        "C#C": "H–C≡C–H",
        "CC#C": "R¹–C≡C–H",
        "CC#CC": "R¹–C≡C–R²",
        "c1ccccc1C#Cc1ccccc1": "R¹–C≡C–R²",
    }

    for smiles, expected in examples.items():
        site = next(
            site
            for site in analyze_molecule(smiles).reactive_site_hypotheses
            if site.site_type == "unsaturated_bond"
        )
        assert site.chemist_label == expected

    acetylene_xh = next(
        site
        for site in analyze_molecule("C#C").reactive_site_hypotheses
        if site.site_type == "pronucleophile_XH"
    )
    assert acetylene_xh.chemist_label == "H–C≡C–H"
    assert acetylene_xh.details["opposite_endpoint_substituent_count"] == 0


def test_unsaturated_labels_support_ascii_style_without_changing_signature() -> None:
    unicode_site = next(
        site
        for site in analyze_molecule("CC#C").reactive_site_hypotheses
        if site.site_type == "unsaturated_bond"
    )
    ascii_site = next(
        site
        for site in analyze_molecule("CC#C", label_style="ascii").reactive_site_hypotheses
        if site.site_type == "unsaturated_bond"
    )

    assert unicode_site.chemist_label == "R¹–C≡C–H"
    assert ascii_site.chemist_label == "R1-C#C-H"
    assert unicode_site.canonical_signature == ascii_site.canonical_signature


def test_organic_nitriles_are_typed_pi_handles() -> None:
    for smiles in ("CC#N", "N#Cc1ccccc1"):
        result = analyze_molecule(smiles)
        site = next(
            site
            for site in result.reactive_site_hypotheses
            if site.canonical_signature == "PI|Nitrile"
        )
        assert site.site_type == "unsaturated_bond"
        assert site.topology == "bond"
        assert site.chemist_label == "R–C≡N"
        assert set(site.details["atom_roles"]) >= {
            "attachment",
            "carbon_endpoint",
            "nitrogen_endpoint",
        }
        assert site.details["endpoint_elements"] == ["C", "N"]
        assert site.details["reaction_modes"] == ["addition", "reduction"]
        assert site.details["electrophilic_endpoint_atom_index"] == (
            site.details["carbon_endpoint_atom_index"]
        )
        group = next(group for group in result.motifs if group.motif_id == "nitrile")
        assert group.chemist_label == "R–C≡N"
        assert {"electrophilic", "reduction_sensitive"} <= set(group.tags)

    ascii_site = next(
        site
        for site in analyze_molecule("CC#N", label_style="ascii").reactive_site_hypotheses
        if site.canonical_signature == "PI|Nitrile"
    )
    assert ascii_site.chemist_label == "R-C#N"


def test_cyanide_and_isocyanide_are_not_organic_nitrile_handles() -> None:
    for smiles in ("[Na+].[C-]#N", "C[N+]#[C-]"):
        assert "PI|Nitrile" not in signatures(smiles)
        assert all(
            group.motif_id != "nitrile"
            for group in analyze_molecule(smiles).motifs
        )


def test_organic_azide_resonance_forms_share_one_dipolar_handle() -> None:
    examples = {
        "CN=[N+]=[N-]": "organic_azide_double_bond_resonance",
        "[N-]=[N+]=Nc1ccccc1": "organic_azide_double_bond_resonance",
        "C[N-][N+]#N": "organic_azide_triple_bond_resonance",
    }
    for smiles, matched_pattern in examples.items():
        result = analyze_molecule(smiles)
        sites = [site for site in result.reactive_site_hypotheses if site.site_type == "dipolar_group"]
        assert len(sites) == 1
        site = sites[0]
        assert site.canonical_signature == "DG|Azide|Organic"
        assert site.chemist_label == "R–N₃"
        assert site.details["matched_pattern"] == matched_pattern
        assert site.details["net_group_charge"] == 0
        assert site.details["reaction_modes"] == ["cycloaddition", "reduction"]
        assert set(site.details["atom_roles"]) >= {
            "attachment",
            "proximal_nitrogen",
            "central_nitrogen",
            "terminal_nitrogen",
        }
        assert any(group.motif_id == "azide" for group in result.motifs)

    ascii_site = next(
        site
        for site in analyze_molecule(
            "CN=[N+]=[N-]", label_style="ascii"
        ).reactive_site_hypotheses
        if site.site_type == "dipolar_group"
    )
    assert ascii_site.chemist_label == "R-N3"


def test_inorganic_azide_is_not_an_organic_dipolar_handle() -> None:
    result = analyze_molecule("[Na+].[N-]=[N+]=[N-]")

    assert result.valid
    assert all(site.site_type != "dipolar_group" for site in result.reactive_site_hypotheses)


def test_organic_heteroatom_pair_bonds_are_typed_handles() -> None:
    examples = {
        "CN=NC": (
            "HB|Azo",
            "R¹–N=N–R²",
            ["isomerization", "reduction"],
            "azo",
        ),
        "CSSC": (
            "HB|Disulfide",
            "R¹–S–S–R²",
            ["exchange", "reduction"],
            "disulfide",
        ),
        "COOC": (
            "HB|Peroxide",
            "R¹–O–O–R²",
            ["homolysis", "reduction"],
            "peroxide",
        ),
    }

    for smiles, (signature, label, modes, group_id) in examples.items():
        result = analyze_molecule(smiles)
        sites = [site for site in result.reactive_site_hypotheses if site.site_type == "heteroatom_bond"]
        assert len(sites) == 1
        site = sites[0]
        assert site.canonical_signature == signature
        assert site.chemist_label == label
        assert site.topology == "bond"
        assert site.details["reaction_modes"] == modes
        assert set(site.details["atom_roles"]) == {
            "attachment_a", "endpoint_a", "endpoint_b", "attachment_b"
        }
        group = next(
            group for group in result.motifs if group.motif_id == group_id
        )
        assert group.chemist_label == label
        if group_id == "peroxide":
            assert all(
                candidate.motif_id != "ether"
                for candidate in result.motifs
            )


def test_heteroatom_pair_handles_require_two_organic_attachments() -> None:
    for smiles in ("NN", "N=N", "SS", "CSC", "OO", "COO", "COC"):
        result = analyze_molecule(smiles)
        assert result.valid, smiles
        assert all(site.site_type != "heteroatom_bond" for site in result.reactive_site_hypotheses)


def test_heteroatom_pair_labels_support_ascii_style() -> None:
    expected = {
        "CN=NC": "R1-N=N-R2",
        "CSSC": "R1-S-S-R2",
        "COOC": "R1-O-O-R2",
    }
    for smiles, label in expected.items():
        site = next(
            site
            for site in analyze_molecule(smiles, label_style="ascii").reactive_site_hypotheses
            if site.site_type == "heteroatom_bond"
        )
        assert site.chemist_label == label


def test_detector_emits_typed_candidates_from_shared_match_index() -> None:
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCN")
    candidates = pronucleophiles.detect(mol, MatchIndex(mol))
    assert len(candidates) == 1
    assert isinstance(candidates[0], ReactiveSiteCandidate)
    assert candidates[0].atom_roles["center"] == (2,)
