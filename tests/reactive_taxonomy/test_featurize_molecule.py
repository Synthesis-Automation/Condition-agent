from reactive_taxonomy import MatchIndex, SiteCandidate, available_styles, featurize_molecule, load_handle_patterns, validate_taxonomy
from reactive_taxonomy.context import load_context_taxonomy
from reactive_taxonomy.sites import pronucleophiles


def signatures(smiles: str) -> set[str]:
    result = featurize_molecule(smiles)
    assert result.valid, result.error
    return {site.canonical_signature for site in result.sites}


def test_taxonomy_bundle_validates() -> None:
    assert validate_taxonomy() == []


def test_handle_smarts_are_independent_and_mapped() -> None:
    patterns = load_handle_patterns()
    assert len(patterns) >= 12
    assert {pattern["site_type"] for pattern in patterns} == {
        "leaving_group", "pronucleophile_XH", "transfer_group",
        "electrophilic_center", "aromatic_CH", "unsaturated_bond",
    }
    assert all(pattern.get("atom_roles") for pattern in patterns)


def test_site_reports_pattern_provenance() -> None:
    site = featurize_molecule("Brc1ccccc1").sites[0]
    assert site.details["matched_pattern"] == "terminal_carbon_halogen"
    assert site.details["alternative_patterns"] == []


def test_context_taxonomy_declares_classification_methods() -> None:
    payload = load_context_taxonomy()
    records = {record["id"]: record for record in payload["contexts"]}
    assert records["HeteroAr"]["classification_method"] == "aromatic_ring_system"
    assert records["Alkyl"]["classification_method"] == "atom_property"
    assert records["C(O)NR"]["classification_method"] == "mapped_smarts"
    assert records["C(O)NR"]["atom_roles"]["context_anchor"] == 1


def test_leaving_groups() -> None:
    assert "LG|Ar|Br" in signatures("Brc1ccccc1")
    assert "LG|HeteroAr|Cl" in signatures("Clc1ncccc1")
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


def test_hydrazine_has_two_sites() -> None:
    nh_sites = [s for s in featurize_molecule("CNN").sites if s.site_type == "pronucleophile_XH" and s.details["center_element"] == "N"]
    assert len(nh_sites) == 2


def test_transfer_groups() -> None:
    assert "TM|Ar|B(OH)2" in signatures("OB(O)c1ccccc1")
    assert "TM|Ar|BPin" in signatures("CC1(C)OB(c2ccccc2)OC1(C)C")


def test_electrophilic_centers() -> None:
    assert "EC|Acyl|Alkyl|OH|latent" in signatures("CC(=O)O")
    assert "EC|Acyl|Alkyl|Cl|activated" in signatures("CC(=O)Cl")
    assert "EC|Sulfonyl|Ar|Cl|activated" in signatures("O=S(=O)(Cl)c1ccccc1")


def test_invalid_input_and_family_filter() -> None:
    assert featurize_molecule("not smiles").error == "INVALID_SMILES"
    result = featurize_molecule("Brc1ccccc1N", site_types={"leaving_group"})
    assert result.valid
    assert {s.site_type for s in result.sites} == {"leaving_group"}


def test_components_are_preserved() -> None:
    result = featurize_molecule("Nc1ccccc1.[K+]")
    assert result.valid
    assert len(result.components) == 2


def test_composite_handles_do_not_emit_internal_sites() -> None:
    result = featurize_molecule("c1ccc(B(O)O)cc1")
    non_aromatic_ch = [site for site in result.sites if site.site_type != "aromatic_CH"]
    assert [site.canonical_signature for site in non_aromatic_ch] == ["TM|Ar|B(OH)2"]


def test_retained_fluorines_are_not_leaving_groups() -> None:
    sigs = {
        site.canonical_signature for site in featurize_molecule("Brc1ccc(C(F)(F)F)cc1").sites
        if site.site_type == "leaving_group"
    }
    assert sigs == {"LG|Ar|Br"}


def test_sn_and_si_emit_only_transferable_carbon_site() -> None:
    tin_sites = [
        site.canonical_signature for site in featurize_molecule("c1ccc([Sn](C)(C)C)cc1").sites
        if site.site_type == "transfer_group"
    ]
    assert tin_sites == ["TM|Ar|SnR3"]
    silicon_sites = {
        site.canonical_signature for site in featurize_molecule("C#C[Si](C)(C)C").sites
        if site.site_type in {"pronucleophile_XH", "transfer_group"}
    }
    assert silicon_sites == {"XH|Csp|H1|Alkynyl", "TM|Alkynyl|SiR3"}


def test_silyl_ether_is_not_a_transfer_group() -> None:
    result = featurize_molecule("Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1")
    assert not [site for site in result.sites if site.site_type == "transfer_group"]
    assert {site.canonical_signature for site in result.sites if site.site_type == "leaving_group"} == {"LG|Ar|Br"}


def test_bridging_metal_halogen_is_not_a_leaving_group() -> None:
    assert "LG|Ar|Br" not in signatures("[Mg]Brc1ccccc1")


def test_acyl_halide_owns_halogen_site() -> None:
    result = featurize_molecule("CC(=O)Cl")
    assert [site.canonical_signature for site in result.sites] == [
        "EC|Acyl|Alkyl|Cl|activated"
    ]


def test_ammonia_is_supported_nh_pronucleophile() -> None:
    result = featurize_molecule("N")
    assert [site.canonical_signature for site in result.sites] == ["XH|N|H3|"]
    assert result.sites[0].chemist_label == "NH3"


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
            site for site in featurize_molecule(smiles).sites
            if site.site_type == "pronucleophile_XH"
        ]
        assert len(sites) == 1
        assert sites[0].canonical_signature == "XH|N_aromatic|H1|HeteroAr"
        assert sites[0].chemist_label == "AromN–H"
        assert sites[0].details["center_element"] == "N"
        assert sites[0].details["derived_family"] == "aromatic_nh"


def test_bromopyrrole_keeps_both_distinct_sites() -> None:
    legacy_sites = {
        site.canonical_signature for site in featurize_molecule("Brc1cc[nH]c1").sites
        if site.site_type in {"leaving_group", "pronucleophile_XH"}
    }
    assert legacy_sites == {
        "LG|HeteroAr|Br",
        "XH|N_aromatic|H1|HeteroAr",
    }


def test_rendering_styles_preserve_signature() -> None:
    unicode_site = featurize_molecule("Nc1ccccc1").sites[0]
    ascii_site = featurize_molecule("Nc1ccccc1", label_style="ascii").sites[0]
    hte_site = featurize_molecule("Nc1ccccc1", label_style="hte_legacy").sites[0]
    assert available_styles() == ("unicode", "ascii", "hte_legacy")
    assert unicode_site.chemist_label == "Ar–NH2"
    assert ascii_site.chemist_label == "Ar-NH2"
    assert hte_site.chemist_label == "Ar-NH2"
    assert len({unicode_site.canonical_signature, ascii_site.canonical_signature, hte_site.canonical_signature}) == 1
    assert featurize_molecule("CC(=O)NC").sites[0].chemist_label == "R–C(O)–NHR"
    assert featurize_molecule("CS(=O)(=O)NC").sites[0].chemist_label == "R–S(O)2–NHR"


def test_availability_distinguishes_chemical_site_state() -> None:
    assert featurize_molecule("CCN").sites[0].availability == "free"
    assert next(site for site in featurize_molecule("CC(=O)NC").sites if site.site_type == "pronucleophile_XH").availability == "deactivated"
    assert next(site for site in featurize_molecule("CC(=O)O").sites if site.site_type == "electrophilic_center").availability == "latent"
    assert featurize_molecule("CC(=O)Cl").sites[0].availability == "activated"
    assert featurize_molecule("OB(O)c1ccccc1").sites[0].availability == "transferable"


def test_rich_context_record_is_exposed() -> None:
    site = featurize_molecule("Nc1ccccn1").sites[0]
    context = site.context_features["contexts"][0]
    assert context["token"] == "HeteroAr"
    assert context["classification_method"] == "aromatic_ring_system"
    assert context["subtype"] == "heteroaromatic_ring"
    assert context["features"]["heteroatoms"] == ["N"]
    assert len(context["fragment_atom_indices"]) == 6


def test_alkyl_leaving_groups_preserve_benzylic_allylic_and_propargylic_subtypes() -> None:
    examples = {
        "ClCc1ccccc1": ("benzylic", "Benzyl–Cl"),
        "ClCC=C": ("allylic", "Allyl–Cl"),
        "ClCC#C": ("propargylic", "Propargyl–Cl"),
        "CCCCl": ("simple_alkyl", "R–Cl"),
    }
    for smiles, (subtype, label) in examples.items():
        site = featurize_molecule(smiles).sites[0]
        context = site.context_features["contexts"][0]
        assert site.details["anchor_context"] == "Alkyl"
        assert site.details["anchor_subtype"] == subtype
        assert context["subtype"] == subtype
        assert site.chemist_label == label


def test_benzyl_chloride_and_aryl_bromide_are_separate_sites() -> None:
    result = featurize_molecule("ClCc1ccccc1Br")
    assert {
        site.chemist_label for site in result.sites if site.site_type == "leaving_group"
    } == {"Benzyl–Cl", "Ar–Br"}


def test_aldehydes_and_ketones_are_carbonyl_addition_centers() -> None:
    aldehyde = next(site for site in featurize_molecule("CC=O").sites if site.site_type == "electrophilic_center")
    assert aldehyde.canonical_signature == "EC|Carbonyl|aldehyde|Alkyl|addition"
    assert aldehyde.details["reaction_mode"] == "addition"
    assert aldehyde.details["atom_roles"] == {"center": [1], "heteroatom": [2], "substituents": [0]}
    assert aldehyde.chemist_label == "R–CH=O"

    ketone = next(site for site in featurize_molecule("CC(=O)C").sites if site.site_type == "electrophilic_center")
    assert ketone.canonical_signature == "EC|Carbonyl|ketone|Alkyl,Alkyl|addition"
    assert ketone.chemist_label == "R2C=O"


def test_acyl_substitution_and_carbonyl_addition_are_distinct() -> None:
    acid_sites = [site for site in featurize_molecule("CC(=O)O").sites if site.site_type == "electrophilic_center"]
    assert [site.details["center_family"] for site in acid_sites] == ["Acyl"]
    assert acid_sites[0].details["reaction_mode"] == "substitution"
    assert not [site for site in featurize_molecule("CC(=O)N").sites if site.site_type == "electrophilic_center"]


def test_aromatic_ch_sites_are_atom_localized_and_ring_classified() -> None:
    benzene_sites = [site for site in featurize_molecule("c1ccccc1").sites if site.site_type == "aromatic_CH"]
    assert len(benzene_sites) == 6
    assert {site.canonical_signature for site in benzene_sites} == {"CH|ArH"}
    assert {site.topology for site in benzene_sites} == {"atom"}
    assert {site.chemist_label for site in benzene_sites} == {"Ar–H"}

    pyridine_sites = [site for site in featurize_molecule("c1ccncc1").sites if site.site_type == "aromatic_CH"]
    assert len(pyridine_sites) == 5
    assert {site.details["handle_token"] for site in pyridine_sites} == {"HetArH"}
    assert {site.chemist_label for site in pyridine_sites} == {"HeteroAr–H"}


def test_unsaturated_carbon_bonds_are_bond_localized() -> None:
    alkene = next(site for site in featurize_molecule("CC=C").sites if site.site_type == "unsaturated_bond")
    assert alkene.topology == "bond"
    assert alkene.canonical_signature == "PI|Alkene"
    assert alkene.details["bond_order"] == 2
    assert set(alkene.details["atom_roles"]) == {"endpoint_a", "endpoint_b"}
    assert alkene.chemist_label == "H2C=CHR1"
    assert alkene.details["endpoint_h_counts"] == [2, 1]
    assert alkene.details["endpoint_substituent_counts"] == [0, 1]

    alkyne_sites = featurize_molecule("CC#C").sites
    assert {site.site_type for site in alkyne_sites} == {"unsaturated_bond", "pronucleophile_XH"}
    assert next(
        site for site in alkyne_sites if site.site_type == "unsaturated_bond"
    ).chemist_label == "R1–C≡C–H"


def test_alkene_labels_expose_all_hydrogen_and_substituent_positions() -> None:
    examples = {
        "C=C": "H2C=CH2",
        "CC=C": "H2C=CHR1",
        "CC=CC": "R1HC=CHR2",
        "C=C(C)c1ccccc1": "H2C=CR1R2",
        "CC(C)=C(C)C": "R1R2C=CR3R4",
    }

    for smiles, expected in examples.items():
        site = next(
            site
            for site in featurize_molecule(smiles).sites
            if site.site_type == "unsaturated_bond"
        )
        assert site.chemist_label == expected
        assert sum(site.details["endpoint_h_counts"]) + sum(
            site.details["endpoint_substituent_counts"]
        ) == 4


def test_alkene_labels_retain_defined_e_z_stereochemistry() -> None:
    trans = next(
        site
        for site in featurize_molecule("c1ccc(/C=C/c2ccccc2)cc1").sites
        if site.site_type == "unsaturated_bond"
    )
    cis = next(
        site
        for site in featurize_molecule("c1ccc(/C=C\\c2ccccc2)cc1").sites
        if site.site_type == "unsaturated_bond"
    )

    assert trans.chemist_label == "R1HC=CHR2 (E)"
    assert cis.chemist_label == "R1HC=CHR2 (Z)"


def test_alkyne_labels_distinguish_acetylene_terminal_and_internal() -> None:
    examples = {
        "C#C": "H–C≡C–H",
        "CC#C": "R1–C≡C–H",
        "CC#CC": "R1–C≡C–R2",
        "c1ccccc1C#Cc1ccccc1": "R1–C≡C–R2",
    }

    for smiles, expected in examples.items():
        site = next(
            site
            for site in featurize_molecule(smiles).sites
            if site.site_type == "unsaturated_bond"
        )
        assert site.chemist_label == expected

    acetylene_xh = next(
        site
        for site in featurize_molecule("C#C").sites
        if site.site_type == "pronucleophile_XH"
    )
    assert acetylene_xh.chemist_label == "H–C≡C–H"
    assert acetylene_xh.details["opposite_endpoint_substituent_count"] == 0


def test_unsaturated_labels_support_ascii_style_without_changing_signature() -> None:
    unicode_site = next(
        site
        for site in featurize_molecule("CC#C").sites
        if site.site_type == "unsaturated_bond"
    )
    ascii_site = next(
        site
        for site in featurize_molecule("CC#C", label_style="ascii").sites
        if site.site_type == "unsaturated_bond"
    )

    assert unicode_site.chemist_label == "R1–C≡C–H"
    assert ascii_site.chemist_label == "R1-C#C-H"
    assert unicode_site.canonical_signature == ascii_site.canonical_signature


def test_detector_emits_typed_candidates_from_shared_match_index() -> None:
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCN")
    candidates = pronucleophiles.detect(mol, MatchIndex(mol))
    assert len(candidates) == 1
    assert isinstance(candidates[0], SiteCandidate)
    assert candidates[0].atom_roles["center"] == (2,)
