from chemtools.reactive_taxonomy import featurize_molecule, load_handle_patterns, validate_taxonomy
from chemtools.reactive_taxonomy.context import load_context_taxonomy


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
        "leaving_group", "pronucleophile_XH", "transfer_group", "electrophilic_center"
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
    assert [site.canonical_signature for site in result.sites] == ["TM|Ar|B(OH)2"]


def test_retained_fluorines_are_not_leaving_groups() -> None:
    sigs = signatures("Brc1ccc(C(F)(F)F)cc1")
    assert sigs == {"LG|Ar|Br"}


def test_sn_and_si_emit_only_transferable_carbon_site() -> None:
    assert signatures("c1ccc([Sn](C)(C)C)cc1") == {"TM|Ar|SnR3"}
    assert signatures("C#C[Si](C)(C)C") == {"XH|Csp|H1|Alkynyl", "TM|Alkynyl|SiR3"}


def test_silyl_ether_is_not_a_transfer_group() -> None:
    sigs = signatures("Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1")
    assert sigs == {"LG|Ar|Br"}


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
