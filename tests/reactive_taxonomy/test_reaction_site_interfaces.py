"""Phase 3 regressions for normalized reactive-site interfaces."""

import pytest

from reactive_taxonomy import (
    BondCapacitySite,
    ConnectionEndpointSite,
    ReactiveLinkSite,
    featurize_molecule,
    normalize_compound_sites,
    normalize_reaction_assignment,
)
from reactive_taxonomy.reaction_candidates import enumerate_reaction_candidates
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def _interfaces(smiles: str):
    analysis = featurize_molecule(smiles)
    assert analysis.valid
    normalized = normalize_compound_sites(analysis)
    assert len(normalized) == len(analysis.sites)
    return analysis, normalized


def _assignment(reaction_smiles: str, grammar_id: str):
    parsed = parse_reaction_smiles(reaction_smiles)
    match = next(
        (
            (grammar, assignment)
            for grammar, assignment in enumerate_reaction_candidates(
                parsed.reactants
            )
            if grammar["id"] == grammar_id
        ),
        None,
    )
    assert match is not None
    return parsed, match[0], match[1]


@pytest.mark.parametrize(
    ("smiles", "expected_elements", "expected_handle"),
    [
        ("Brc1ccccc1", ("C", "Br"), "Br"),
        ("COS(=O)(=O)C", ("C", "O"), "OMs"),
    ],
)
def test_leaving_groups_adapt_to_explicit_reactive_links(
    smiles: str,
    expected_elements: tuple[str, str],
    expected_handle: str,
) -> None:
    analysis, normalized = _interfaces(smiles)
    pair = next(
        (site, view)
        for site, view in zip(analysis.sites, normalized)
        if site.site_type == "leaving_group"
    )
    source, view = pair
    link = view.reactive_links[0]

    assert isinstance(link, ReactiveLinkSite)
    assert link.source_kind == "explicit_bond"
    assert (link.endpoint_a.element, link.endpoint_b.element) == expected_elements
    assert link.before_order == "SINGLE"
    assert link.endpoint_b.source_atom_role == "handle"
    assert source.site_id == link.source_site_id
    assert source.site_type in link.annotation_tokens
    assert source.canonical_signature in link.annotation_tokens
    assert source.chemist_label == link.source_chemist_label
    assert source.chemist_label in link.annotation_tokens
    assert expected_handle in link.annotation_tokens


def test_transfer_and_addition_modes_keep_distinct_source_annotations() -> None:
    analysis, normalized = _interfaces("c1ccccc1[SiH](C)C")
    by_type = {
        site.site_type: view
        for site, view in zip(analysis.sites, normalized)
        if site.site_type in {"transfer_group", "addition_donor"}
    }

    transfer = by_type["transfer_group"].reactive_links[0]
    addition = by_type["addition_donor"].reactive_links[0]

    assert (
        transfer.endpoint_a.element,
        transfer.endpoint_b.element,
        transfer.source_kind,
    ) == ("C", "Si", "explicit_bond")
    assert (
        addition.endpoint_a.element,
        addition.endpoint_b.element,
        addition.source_kind,
    ) == ("Si", "H", "implicit_hydrogen")
    assert transfer.source_site_type == "transfer_group"
    assert addition.source_site_type == "addition_donor"
    assert transfer.source_site_id != addition.source_site_id


def test_boron_transfer_group_adapts_without_losing_carrier_annotation() -> None:
    analysis, normalized = _interfaces("OB(O)c1ccccc1")
    pair = next(
        (site, view)
        for site, view in zip(analysis.sites, normalized)
        if site.site_type == "transfer_group"
    )
    source, view = pair
    link = view.reactive_links[0]

    assert (link.endpoint_a.element, link.endpoint_b.element) == ("C", "B")
    assert link.source_kind == "explicit_bond"
    assert "B(OH)2" in link.annotation_tokens
    assert link.source_chemist_label == source.chemist_label


@pytest.mark.parametrize(
    ("smiles", "element", "available_units"),
    [
        ("CN", "N", 2),
        ("CO", "O", 1),
        ("CS", "S", 1),
        ("B", "B", 3),
        ("C[SiH](C)C", "Si", 1),
    ],
)
def test_xh_detectors_share_one_virtual_hydrogen_link_contract(
    smiles: str,
    element: str,
    available_units: int,
) -> None:
    analysis, normalized = _interfaces(smiles)
    view = next(item for item in normalized if item.reactive_links)
    link = view.reactive_links[0]
    connection = view.connection_endpoints[0]

    assert link.source_kind == "implicit_hydrogen"
    assert link.endpoint_a.element == element
    assert link.endpoint_b.endpoint_kind == "virtual_hydrogen"
    assert link.endpoint_b.element == "H"
    assert link.endpoint_b.atom_index is None
    assert (
        link.endpoint_b.carrier_atom_index
        == link.endpoint_a.atom_index
    )
    assert link.available_units == available_units
    assert link.symmetry_class == "virtual_hydrogen"
    assert isinstance(connection, ConnectionEndpointSite)
    assert connection.required_hydrogen_delta == -1
    assert connection.required_link_release_id is None


@pytest.mark.parametrize(
    ("smiles", "symmetry"),
    [
        ("BrBr", "equivalent_endpoints"),
        ("ClBr", "distinguishable_endpoints"),
        ("BB", "equivalent_endpoints"),
        ("[Si]B", "distinguishable_endpoints"),
    ],
)
def test_explicit_ab_donors_preserve_endpoint_symmetry(
    smiles: str,
    symmetry: str,
) -> None:
    _, normalized = _interfaces(smiles)
    view, link = next(
        (view, link)
        for view in normalized
        for link in view.reactive_links
        if link.source_kind == "explicit_bond"
    )

    assert link.source_kind == "explicit_bond"
    assert link.symmetry_class == symmetry
    assert all(
        endpoint.required_link_release_id == link.site_id
        for endpoint in view.connection_endpoints
    )


@pytest.mark.parametrize(
    ("smiles", "order", "decrement", "increment", "bond_class"),
    [
        ("C=C", "DOUBLE", 1, 1, "carbon_carbon_pi"),
        ("C#C", "TRIPLE", 2, 0, "carbon_carbon_pi"),
        ("CC#N", "TRIPLE", 2, 0, "carbon_nitrogen_pi"),
    ],
)
def test_unsaturated_bonds_adapt_to_bounded_capacity(
    smiles: str,
    order: str,
    decrement: int,
    increment: int,
    bond_class: str,
) -> None:
    _, normalized = _interfaces(smiles)
    capacity = next(
        view.bond_capacities[0]
        for view in normalized
        if view.bond_capacities
    )

    assert isinstance(capacity, BondCapacitySite)
    assert capacity.current_order == order
    assert capacity.maximum_decrement == decrement
    assert capacity.maximum_increment == increment
    assert capacity.bond_class == bond_class
    assert not capacity.aromatic
    assert capacity.endpoint_a.local_environment
    assert capacity.endpoint_b.local_environment


def test_capacity_endpoints_declare_required_bond_unit_consumption() -> None:
    _, normalized = _interfaces("CC=C")
    view = next(item for item in normalized if item.bond_capacities)
    capacity = view.bond_capacities[0]

    assert {
        endpoint.required_bond_capacity_id
        for endpoint in view.connection_endpoints
    } == {capacity.site_id}
    assert {
        endpoint.required_bond_capacity_decrement
        for endpoint in view.connection_endpoints
    } == {1}


def test_inert_atoms_are_not_promoted_to_generic_connection_candidates() -> None:
    analysis, normalized = _interfaces("CCC")

    assert analysis.sites == []
    assert normalized == ()


def test_normalized_ids_disambiguate_reaction_components() -> None:
    parsed, _, assignment = _assignment(
        "C=C.BrBr>>BrCCBr",
        "addend_pair_addition_to_alkene",
    )
    normalized = normalize_reaction_assignment(assignment, parsed.reactants)
    interface_ids = {
        item.site_id
        for view in normalized.values()
        for item in (
            *view.reactive_links,
            *view.bond_capacities,
            *view.connection_endpoints,
        )
    }

    assert len(interface_ids) == sum(
        len(view.reactive_links)
        + len(view.bond_capacities)
        + len(view.connection_endpoints)
        for view in normalized.values()
    )
    assert any("component0" in value for value in interface_ids)
    assert any("component1" in value for value in interface_ids)


def test_normalized_views_do_not_change_public_molecule_serialization() -> None:
    analysis, normalized = _interfaces("C=C.BrBr")
    before = analysis.to_dict()

    assert normalized
    assert analysis.to_dict() == before
    assert all(
        "reactive_links" not in site_payload
        and "bond_capacities" not in site_payload
        and "connection_endpoints" not in site_payload
        for site_payload in before["sites"]
    )


def test_interface_chemistry_is_component_order_invariant() -> None:
    forward = _assignment(
        "C=C.C[SiH](C)C>>CC[Si](C)(C)C",
        "addend_pair_addition_to_alkene",
    )
    reversed_partners = _assignment(
        "C[SiH](C)C.C=C>>CC[Si](C)(C)C",
        "addend_pair_addition_to_alkene",
    )

    def chemistry_key(case) -> tuple:
        parsed, _, assignment = case
        normalized = normalize_reaction_assignment(
            assignment, parsed.reactants
        )
        return tuple(
            sorted(
                (
                    role,
                    tuple(
                        (
                            link.source_kind,
                            link.endpoint_a.element,
                            link.endpoint_b.element,
                            link.before_order,
                            link.available_units,
                        )
                        for link in view.reactive_links
                    ),
                    tuple(
                        (
                            capacity.endpoint_a.element,
                            capacity.endpoint_b.element,
                            capacity.current_order,
                            capacity.maximum_decrement,
                        )
                        for capacity in view.bond_capacities
                    ),
                )
                for role, view in normalized.items()
            )
        )

    assert chemistry_key(forward) == chemistry_key(reversed_partners)


def test_activated_acyl_center_exposes_one_releasable_link() -> None:
    parsed, _, assignment = _assignment(
        "CC(=O)Cl.CN>>CC(=O)NC",
        "amide_formation",
    )
    normalized = normalize_reaction_assignment(
        assignment, parsed.reactants
    )
    link = normalized["acyl_partner"].reactive_links[0]

    assert (link.endpoint_a.element, link.endpoint_b.element) == ("C", "Cl")
    assert (
        link.endpoint_a.source_atom_role,
        link.endpoint_b.source_atom_role,
    ) == ("center", "leaving_or_activatable")


@pytest.mark.parametrize(
    ("reaction", "grammar_id", "role", "bond_class", "order"),
    [
        (
            "CC(=O)C>>CC(O)C",
            "carbonyl_reduction",
            "carbonyl",
            "polarized_multiple_bond",
            "DOUBLE",
        ),
        (
            "CC(O)C>>CC(=O)C",
            "alcohol_oxidation",
            "alcohol",
            "alcohol_carbon_heteroatom",
            "SINGLE",
        ),
    ],
)
def test_simple_polar_order_changes_expose_bond_capacity(
    reaction: str,
    grammar_id: str,
    role: str,
    bond_class: str,
    order: str,
) -> None:
    parsed, _, assignment = _assignment(reaction, grammar_id)
    normalized = normalize_reaction_assignment(
        assignment, parsed.reactants
    )
    capacity = normalized[role].bond_capacities[0]

    assert capacity.bond_class == bond_class
    assert capacity.current_order == order


def test_eliminable_pair_exposes_link_capacity_and_hydrogen_endpoint() -> None:
    parsed, _, assignment = _assignment(
        "CCC(C)Br>>CC=CC",
        "beta_halo_elimination",
    )
    normalized = normalize_reaction_assignment(
        assignment, parsed.reactants
    )["substrate"]

    assert len(normalized.reactive_links) == 1
    assert len(normalized.bond_capacities) == 1
    assert len(normalized.connection_endpoints) == 1
    assert normalized.reactive_links[0].endpoint_b.source_atom_role == (
        "departing_a"
    )
    assert normalized.bond_capacities[0].bond_class == (
        "eliminable_backbone"
    )
    assert normalized.connection_endpoints[0].required_hydrogen_delta == -1
    assert (
        normalized.connection_endpoints[0].endpoint.source_atom_role
        == "hydrogen_carrier_b"
    )


def test_anionic_partner_exposes_charge_normalizing_connection() -> None:
    parsed, _, assignment = _assignment(
        "CC(=O)[S-].[K+].Ic1ccccc1>>CC(=O)Sc1ccccc1",
        "sp2_c_s_anion_substitution",
    )
    normalized = normalize_reaction_assignment(
        assignment, parsed.reactants
    )["nucleophile"]

    assert not normalized.reactive_links
    assert len(normalized.connection_endpoints) == 1
    endpoint = normalized.connection_endpoints[0]
    assert endpoint.endpoint.formal_charge == -1
    assert endpoint.required_formal_charge_delta == 1
