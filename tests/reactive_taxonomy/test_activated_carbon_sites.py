"""Regression coverage for explicitly activated sp3 C-H pronucleophiles."""

from reactive_taxonomy import featurize_molecule


def _activated_sites(smiles: str):
    result = featurize_molecule(smiles)
    assert result.valid, result.error
    return [
        site
        for site in result.sites
        if site.site_type == "pronucleophile_XH"
        and site.details.get("center_token") == "Csp3"
    ]


def test_common_electron_withdrawing_groups_activate_sp3_c_h() -> None:
    examples = {
        "CC(=O)c1ccccc1": "XH|Csp3|H3|alpha_to:ketone",
        "CC(=O)OCC": "XH|Csp3|H3|alpha_to:ester",
        "CC(=O)NC": "XH|Csp3|H3|alpha_to:amide",
        "CC#N": "XH|Csp3|H3|alpha_to:nitrile",
        "CS(=O)(=O)C": "XH|Csp3|H3|alpha_to:sulfone",
        "C[N+](=O)[O-]": "XH|Csp3|H3|alpha_to:nitro",
    }

    for smiles, signature in examples.items():
        sites = _activated_sites(smiles)
        assert signature in {site.canonical_signature for site in sites}
        site = next(
            item for item in sites if item.canonical_signature == signature
        )
        assert site.availability == "activated"
        assert site.chemist_label == "activated C–H"
        assert site.details["activation_relationship"] == "alpha_to"
        assert len(site.details["activation_records"]) == 1


def test_doubly_activated_carbon_is_one_site_with_two_evidence_records() -> None:
    sites = _activated_sites("COC(=O)CC(=O)OC")

    assert len(sites) == 1
    assert sites[0].canonical_signature == (
        "XH|Csp3|H2|alpha_to:ester*2"
    )
    assert sites[0].details["activation_count"] == 2
    assert len(sites[0].details["activation_records"]) == 2
    assert len(sites[0].details["activation_anchor_atom_indices"]) == 2


def test_symmetry_equivalent_activated_carbons_remain_atom_localized() -> None:
    sites = _activated_sites("CC(=O)C")

    assert len(sites) == 2
    assert {site.canonical_signature for site in sites} == {
        "XH|Csp3|H3|alpha_to:ketone"
    }
    assert len({tuple(site.atom_indices) for site in sites}) == 2


def test_unactivated_and_other_carbon_hydrogens_are_not_reclassified() -> None:
    assert _activated_sites("CCC") == []
    assert _activated_sites("Cc1ccccc1") == []

    terminal_alkyne = featurize_molecule("CC#C")
    carbon_sites = [
        site
        for site in terminal_alkyne.sites
        if site.site_type == "pronucleophile_XH"
    ]
    assert [site.canonical_signature for site in carbon_sites] == [
        "XH|Csp|H1|Alkynyl"
    ]


def test_activated_site_signature_is_serialization_invariant() -> None:
    forward = {
        site.canonical_signature for site in _activated_sites("CC(=O)OCC")
    }
    reverse = {
        site.canonical_signature for site in _activated_sites("CCOC(C)=O")
    }

    assert forward == reverse == {"XH|Csp3|H3|alpha_to:ester"}
