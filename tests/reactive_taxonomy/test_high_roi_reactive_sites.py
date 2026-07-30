"""Conservative high-ROI reactive-site coverage regressions."""

from __future__ import annotations

import pytest

from reactive_taxonomy import featurize_molecule
from reactive_taxonomy.reaction_candidates import enumerate_reaction_candidates
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def _sites(smiles: str, site_type: str):
    analysis = featurize_molecule(smiles)
    assert analysis.valid, analysis.error
    return tuple(
        (site, interface)
        for site, interface in zip(
            analysis.sites,
            analysis.connectivity_sites,
        )
        if site.site_type == site_type
    )


@pytest.mark.parametrize(
    ("smiles", "element"),
    (
        ("[Na+].[O-]CC", "O"),
        ("[NH2-]", "N"),
        ("[CH2-]C#N", "C"),
        ("C[S-]", "S"),
    ),
)
def test_explicit_anions_emit_charge_balanced_connection_endpoints(
    smiles: str,
    element: str,
) -> None:
    values = _sites(smiles, "nucleophile_anion")

    assert len(values) == 1
    site, interface = values[0]
    assert site.details["center_token"] == element
    assert site.details["formal_charge"] == -1
    assert site.availability == "ionic"
    assert not interface.reactive_links
    assert not interface.bond_capacities
    assert len(interface.connection_endpoints) == 1
    endpoint = interface.connection_endpoints[0]
    assert endpoint.endpoint.element == element
    assert endpoint.required_formal_charge_delta == 1


def test_anionic_sites_are_not_inferred_from_neutral_precursors() -> None:
    for smiles in ("CCO", "CCN", "CCC#N"):
        assert not _sites(smiles, "nucleophile_anion")


@pytest.mark.parametrize("smiles", ("CC=NCC", "CC=NO", "CC=NNC"))
def test_c_n_double_bond_is_a_polarized_capacity_without_alkene_identity(
    smiles: str,
) -> None:
    values = _sites(smiles, "unsaturated_bond")

    assert len(values) == 1
    site, interface = values[0]
    assert site.canonical_signature == "PI|PolarizedC=N"
    assert site.details["electrophilic_endpoint_atom_index"] in site.atom_indices
    assert site.details["endpoint_elements"] == ["C", "N"]
    assert not interface.reactive_links
    assert len(interface.bond_capacities) == 1
    capacity = interface.bond_capacities[0]
    assert capacity.current_order == "DOUBLE"
    assert capacity.maximum_decrement == 1
    assert capacity.bond_class == "localized_multiple_bond"
    assert {
        endpoint.endpoint.element
        for endpoint in interface.connection_endpoints
    } == {"C", "N"}


def test_resonance_anions_are_not_free_nucleophile_endpoints() -> None:
    for smiles in ("O=[N+]([O-])c1ccccc1", "CN=[N+]=[N-]"):
        assert not _sites(smiles, "nucleophile_anion")


def test_nitrile_and_carbonyl_do_not_collapse_to_polarized_c_n() -> None:
    for smiles in ("CC#N", "CC=O", "CC(=O)N"):
        assert all(
            site.canonical_signature != "PI|PolarizedC=N"
            for site in featurize_molecule(smiles).sites
        )


@pytest.mark.parametrize(
    ("smiles", "subtype", "heteroatom"),
    (
        ("C1CO1", "epoxide", "O"),
        ("C1CN1", "aziridine", "N"),
    ),
)
def test_strained_three_membered_rings_expose_both_opening_centers(
    smiles: str,
    subtype: str,
    heteroatom: str,
) -> None:
    values = _sites(smiles, "electrophilic_center")
    strained = tuple(
        value
        for value in values
        if value[0].details.get("center_family") == "StrainedRing"
    )

    assert len(strained) == 2
    assert len({site.details["center_atom_index"] for site, _ in strained}) == 2
    for site, interface in strained:
        assert site.details["strained_ring_type"] == subtype
        assert site.details["reaction_mode"] == "ring_opening"
        assert site.availability == "activated"
        assert len(interface.reactive_links) == 1
        assert not interface.bond_capacities
        link = interface.reactive_links[0]
        assert link.before_order == "SINGLE"
        assert (link.endpoint_a.element, link.endpoint_b.element) == (
            "C",
            heteroatom,
        )


def test_cyclopropane_is_not_a_strained_heterocycle_site() -> None:
    assert not _sites("C1CC1", "electrophilic_center")


def test_silyl_ether_exposes_only_the_oxygen_silicon_release_link() -> None:
    values = _sites("CCO[Si](C)(C)C", "leaving_group")

    assert len(values) == 1
    site, interface = values[0]
    assert site.canonical_signature == "LG|O|SiR3"
    assert site.details["release_mode"] == "deprotection"
    assert site.availability == "conditional"
    assert len(interface.reactive_links) == 1
    link = interface.reactive_links[0]
    assert (link.endpoint_a.element, link.endpoint_b.element) == ("O", "Si")
    assert link.before_order == "SINGLE"


def test_silanol_and_siloxane_are_not_silyl_ether_release_sites() -> None:
    for smiles in ("O[Si](C)(C)C", "C[Si](C)(C)O[Si](C)(C)C"):
        assert all(
            site.canonical_signature != "LG|O|SiR3"
            for site in featurize_molecule(smiles).sites
        )


@pytest.mark.parametrize(
    ("smiles", "token"),
    (
        ("C[Li]", "Li"),
        ("C[Cu]", "Cu"),
        ("C[Al](C)C", "Al"),
    ),
)
def test_additional_organometallic_links_are_transfer_sites(
    smiles: str,
    token: str,
) -> None:
    values = _sites(smiles, "transfer_group")

    assert values
    for site, interface in values:
        assert site.details["handle_token"] == token
        assert site.availability == "transferable"
        assert len(interface.reactive_links) == 1
        link = interface.reactive_links[0]
        assert link.endpoint_a.element == "C"
        assert link.endpoint_b.element == token


def test_dibal_carbon_groups_are_not_marked_as_transfer_partners() -> None:
    assert not _sites("CC(C)C[AlH]CC(C)C", "transfer_group")


def test_new_observations_do_not_create_unregistered_reaction_grammars() -> None:
    parsed = parse_reaction_smiles("C1CO1.[O-]C>>COCCO")

    assert parsed.valid
    assert enumerate_reaction_candidates(parsed.reactants) == []
