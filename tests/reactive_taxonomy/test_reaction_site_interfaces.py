"""Reactive-site interfaces are optional molecular annotations."""

from reactive_taxonomy import (
    BondCapacitySite,
    ConnectionEndpointSite,
    ReactiveLinkSite,
    analyze_molecule,
    normalize_compound_sites,
)


def _interfaces(smiles: str):
    analysis = analyze_molecule(smiles)
    assert analysis.valid
    normalized = normalize_compound_sites(analysis)
    assert len(normalized) == len(analysis.reactive_site_hypotheses)
    return analysis, normalized


def test_leaving_group_exposes_an_optional_reactive_link() -> None:
    analysis, normalized = _interfaces("Brc1ccccc1")
    source, view = next(
        (site, interface)
        for site, interface in zip(
            analysis.reactive_site_hypotheses, normalized
        )
        if site.site_type == "leaving_group"
    )
    link = view.reactive_links[0]

    assert isinstance(link, ReactiveLinkSite)
    assert link.source_kind == "explicit_bond"
    assert (link.endpoint_a.element, link.endpoint_b.element) == ("C", "Br")
    assert link.source_site_id == source.hypothesis_id


def test_boron_transfer_annotation_retains_carrier_identity() -> None:
    analysis, normalized = _interfaces("OB(O)c1ccccc1")
    source, view = next(
        (site, interface)
        for site, interface in zip(
            analysis.reactive_site_hypotheses, normalized
        )
        if site.site_type == "transfer_group"
    )
    link = view.reactive_links[0]

    assert (link.endpoint_a.element, link.endpoint_b.element) == ("C", "B")
    assert "B(OH)2" in link.annotation_tokens
    assert link.source_chemist_label == source.chemist_label


def test_xh_annotation_uses_a_virtual_hydrogen_endpoint() -> None:
    _, normalized = _interfaces("CN")
    view = next(interface for interface in normalized if interface.reactive_links)
    link = view.reactive_links[0]
    endpoint = view.connection_endpoints[0]

    assert link.endpoint_a.element == "N"
    assert link.endpoint_b.endpoint_kind == "virtual_hydrogen"
    assert link.endpoint_b.element == "H"
    assert isinstance(endpoint, ConnectionEndpointSite)
    assert endpoint.required_hydrogen_delta == -1


def test_unsaturated_bond_exposes_bounded_capacity() -> None:
    _, normalized = _interfaces("C#C")
    capacity = next(
        interface.bond_capacities[0]
        for interface in normalized
        if interface.bond_capacities
    )

    assert isinstance(capacity, BondCapacitySite)
    assert capacity.current_order == "TRIPLE"
    assert capacity.maximum_decrement == 2
    assert capacity.maximum_increment == 0
    assert capacity.bond_class == "carbon_carbon_pi"


def test_inert_atoms_are_not_promoted_to_reactive_sites() -> None:
    analysis, normalized = _interfaces("CCC")

    assert analysis.reactive_site_hypotheses == ()
    assert normalized == ()
