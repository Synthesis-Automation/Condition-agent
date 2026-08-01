"""Regression coverage for reusable substitution/addition/elimination archetypes."""

from reactive_taxonomy import featurize_molecule, featurize_reaction
from reactive_taxonomy.connectivity_rewrite import load_connectivity_rewrites
from reactive_taxonomy.reaction_reconstruction_rules import (
    load_reaction_reconstruction_rules,
)


def test_every_reconstruction_rule_has_a_registered_operator() -> None:
    rules = load_reaction_reconstruction_rules()
    rewrite_ids = {
        rewrite.rewrite_id for rewrite in load_connectivity_rewrites()
    }

    assert rules
    assert {rule["operator_id"] for rule in rules} <= rewrite_ids
    assert all("transformation_class" not in rule for rule in rules)
    assert all("compatible_named_families" not in rule for rule in rules)


def test_explicit_and_implicit_addition_donors_are_distinct_sites() -> None:
    bromine = featurize_molecule("BrBr")
    hydrosilane = featurize_molecule("C[SiH](C)C")
    borane = featurize_molecule("B")

    assert [site.canonical_signature for site in bromine.sites] == [
        "AD|ATOM_ATOM|Br-Br|SINGLE"
    ]
    assert [site.canonical_signature for site in hydrosilane.sites] == [
        "AD|ATOM_HYDROGEN|Si-H|H1"
    ]
    assert [site.canonical_signature for site in borane.sites] == [
        "AD|ATOM_HYDROGEN|B-H|H3"
    ]


def test_hydrosilane_can_expose_transfer_and_addition_modes() -> None:
    result = featurize_molecule("c1ccccc1[SiH](C)C")

    assert {
        (site.site_type, site.canonical_signature)
        for site in result.sites
        if site.site_type in {"transfer_group", "addition_donor"}
    } == {
        ("transfer_group", "TM|Ar|SiR3"),
        ("addition_donor", "AD|ATOM_HYDROGEN|Si-H|H1"),
    }


def test_xh_addition_reconstructs_both_regioisomers_without_guessing() -> None:
    linear = featurize_reaction("CC=C.CN>>CCCNC")
    branched = featurize_reaction("CC=C.CN>>CC(C)NC")

    assert linear.selected_candidate is not None
    assert branched.selected_candidate is not None
    assert linear.selected_candidate.grammar_id == "xh_addition_to_alkene"
    assert branched.selected_candidate.grammar_id == "xh_addition_to_alkene"
    assert (
        linear.selected_candidate.rewrite_outcome_id
        != branched.selected_candidate.rewrite_outcome_id
    )
    assert linear.edit_archetype == branched.edit_archetype == "addition"
    assert linear.named_family is None
    assert branched.named_family is None
    assert {
        edit.edit_type for edit in linear.reaction_signature.edits
    } == {"formed", "order_changed", "hydrogen_change"}
    assert {
        candidate.rewrite_outcome_id
        for candidate in linear.candidates
        if candidate.grammar_id == "xh_addition_to_alkene"
    } == {
        "endpoint_a_addend_a__endpoint_b_addend_b",
        "endpoint_b_addend_a__endpoint_a_addend_b",
    }


def test_xh_addition_mapping_agreement_and_conflict_are_preserved() -> None:
    agreed = featurize_reaction(
        "[CH3:1][CH:2]=[CH2:3].[NH2:4][CH3:5]>>"
        "[CH3:1][CH2:2][CH2:3][NH:4][CH3:5]"
    )
    conflicted = featurize_reaction(
        "[CH3:1][CH:2]=[CH2:3].[NH2:4][CH3:5]>>"
        "[CH3:1][CH2:3][CH2:2][NH:4][CH3:5]"
    )

    assert agreed.evidence_quality == (
        "validated_mapping_and_exact_reconstruction"
    )
    assert conflicted.evidence_quality == "conflicting_edit_evidence"
    assert "MAPPING_RECONSTRUCTION_CONFLICT" in conflicted.warnings
    assert conflicted.reaction_signature is not None
    assert conflicted.reaction_signature.edit_archetype == "addition"


def test_dihalogen_and_hydrosilane_use_the_same_pair_addition_rewrite() -> None:
    bromination = featurize_reaction("C=C.BrBr>>BrCCBr")
    hydrosilylation = featurize_reaction(
        "C=C.C[SiH](C)C>>CC[Si](C)(C)C"
    )

    assert bromination.selected_candidate is not None
    assert hydrosilylation.selected_candidate is not None
    assert bromination.selected_candidate.grammar_id == (
        "addend_pair_addition_to_alkene"
    )
    assert hydrosilylation.selected_candidate.grammar_id == (
        "addend_pair_addition_to_alkene"
    )
    assert bromination.edit_archetype == hydrosilylation.edit_archetype == (
        "addition"
    )
    assert {
        edit.edit_type for edit in bromination.reaction_signature.edits
    } == {"formed", "broken", "order_changed"}
    assert {
        edit.edit_type for edit in hydrosilylation.reaction_signature.edits
    } == {"formed", "order_changed", "hydrogen_change"}


def test_one_equivalent_alkyne_pair_addition_stops_at_alkene() -> None:
    result = featurize_reaction("C#C.BrBr>>BrC=CBr")

    assert result.selected_candidate is not None
    assert result.selected_candidate.grammar_id == (
        "addend_pair_addition_to_alkyne"
    )
    assert result.reaction_signature.order_changes == (
        "C-C:TRIPLE>DOUBLE",
    )


def test_pair_addition_signature_is_partner_order_invariant() -> None:
    forward = featurize_reaction("C=C.BrBr>>BrCCBr")
    reversed_partners = featurize_reaction("BrBr.C=C>>BrCCBr")

    assert forward.reaction_signature is not None
    assert reversed_partners.reaction_signature is not None
    assert (
        forward.reaction_signature.signature_id
        == reversed_partners.reaction_signature.signature_id
    )


def test_beta_elimination_resolves_each_available_beta_site() -> None:
    internal = featurize_reaction("CCC(C)Br>>CC=CC")
    terminal = featurize_reaction("CCC(C)Br>>C=CCC")

    assert internal.selected_candidate is not None
    assert terminal.selected_candidate is not None
    assert internal.selected_candidate.grammar_id == "beta_halo_elimination"
    assert terminal.selected_candidate.grammar_id == "beta_halo_elimination"
    assert internal.edit_archetype == terminal.edit_archetype == "elimination"
    assert {
        edit.edit_type for edit in internal.reaction_signature.edits
    } == {"broken", "order_changed", "hydrogen_change"}


def test_beta_elimination_requires_an_adjacent_hydrogen() -> None:
    molecule = featurize_molecule("CC(C)(C)CBr")
    result = featurize_reaction("CC(C)(C)CBr>>CC(C)(C)C")

    assert all(site.site_type != "eliminable_pair" for site in molecule.sites)
    assert all(
        candidate.grammar_id != "beta_halo_elimination"
        for candidate in result.candidates
    )
