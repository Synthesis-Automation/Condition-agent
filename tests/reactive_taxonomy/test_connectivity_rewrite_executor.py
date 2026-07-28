"""Phase 2 regression tests for the bounded connectivity rewrite executor."""

from dataclasses import replace

import pytest

from reactive_taxonomy import (
    apply_connectivity_rewrite,
    compare_connectivity_rewrite,
    compile_connectivity_rewrite_definitions,
    load_connectivity_rewrites,
    reaction_signature_definition_versions,
)
from reactive_taxonomy.reaction_candidates import enumerate_reaction_candidates
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def _assignment(
    reaction_smiles: str,
    grammar_id: str,
):
    parsed = parse_reaction_smiles(reaction_smiles)
    matches = [
        (grammar, assignment)
        for grammar, assignment in enumerate_reaction_candidates(parsed.reactants)
        if grammar["id"] == grammar_id
    ]
    assert matches
    return parsed, matches[0][0], matches[0][1]


def test_rewrite_definitions_are_separately_versioned_and_compiled() -> None:
    rewrites = load_connectivity_rewrites()

    assert {rewrite.rewrite_id for rewrite in rewrites} == {
        "alkene_split_and_distribute",
        "beta_depart_and_unsaturate",
        "c_n_release_and_connect",
        "suzuki_release_and_connect",
    }
    assert {
        grammar_id
        for rewrite in rewrites
        for grammar_id in rewrite.grammar_ids
    } == {
        "addend_pair_addition_to_alkene",
        "beta_halo_elimination",
        "boron_transfer_coupling",
        "sp2_c_n_substitution",
        "sp3_c_n_substitution",
    }
    assert (
        "connectivity_rewrites.v1.json"
        not in reaction_signature_definition_versions()
    )


@pytest.mark.parametrize(
    "instruction",
    [
        {"op": "run_python", "callable": "module.function"},
        {
            "op": "change_localized_bond_state",
            "endpoints": ["left.anchor", "right.anchor"],
            "before": "AROMATIC",
            "after": "SINGLE",
        },
    ],
)
def test_compiler_rejects_unbounded_or_unsupported_instructions(
    instruction,
) -> None:
    payload = {
        "schema_version": "1.0",
        "instruction_set_version": "1.0",
        "rewrites": [
            {
                "id": "invalid",
                "template": "release_and_connect",
                "grammar_ids": ["example"],
                "variants": [
                    {
                        "id": "default",
                        "instructions": [
                            instruction,
                            {
                                "op": "declare_product_seed",
                                "selector": "left.anchor",
                            },
                        ],
                    }
                ],
            }
        ],
    }

    with pytest.raises(ValueError):
        compile_connectivity_rewrite_definitions(payload)


def test_compiler_rejects_projection_without_a_declared_bond_break() -> None:
    payload = {
        "schema_version": "1.0",
        "instruction_set_version": "1.0",
        "rewrites": [
            {
                "id": "unsafe_projection",
                "template": "release_and_connect",
                "grammar_ids": ["example"],
                "variants": [
                    {
                        "id": "default",
                        "instructions": [
                            {
                                "op": (
                                    "declare_projection_discardable_attachment"
                                ),
                                "retained": "left.anchor",
                                "discarded": "right.handle",
                            },
                            {
                                "op": "declare_product_seed",
                                "selector": "left.anchor",
                            },
                        ],
                    }
                ],
            }
        ],
    }

    with pytest.raises(ValueError, match="declared broken attachment"):
        compile_connectivity_rewrite_definitions(payload)


@pytest.mark.parametrize(
    ("reaction_smiles", "grammar_id"),
    [
        (
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            "boron_transfer_coupling",
        ),
        ("Brc1ccccc1.CN>>CNc1ccccc1", "sp2_c_n_substitution"),
        ("CCBr.CN>>CCNC", "sp3_c_n_substitution"),
        ("C=C.BrBr>>BrCCBr", "addend_pair_addition_to_alkene"),
        (
            "C=C.C[SiH](C)C>>CC[Si](C)(C)C",
            "addend_pair_addition_to_alkene",
        ),
        ("CC=C.ClBr>>CC(Cl)CBr", "addend_pair_addition_to_alkene"),
        ("C=C.B>>CCB", "addend_pair_addition_to_alkene"),
        ("CCC(C)Br>>CC=CC", "beta_halo_elimination"),
    ],
)
def test_pilot_rewrites_have_exact_legacy_parity(
    reaction_smiles: str,
    grammar_id: str,
) -> None:
    parsed, grammar, assignment = _assignment(reaction_smiles, grammar_id)

    report = compare_connectivity_rewrite(
        grammar, assignment, parsed.reactants
    )

    assert report.migrated
    assert report.parity
    assert report.product_parity
    assert report.edit_parity
    assert report.warning_parity
    assert report.ambiguity_parity
    assert report.ordering_parity
    assert report.rewrite_outcomes == report.legacy_outcomes


def test_pair_addition_enumerates_and_orders_both_regioisomers() -> None:
    parsed, grammar, assignment = _assignment(
        "CC=C.C[SiH](C)C>>CCC[Si](C)(C)C",
        "addend_pair_addition_to_alkene",
    )

    outcomes = apply_connectivity_rewrite(
        grammar, assignment, parsed.reactants
    )

    assert [
        (outcome.outcome_id, outcome.predicted_product_smiles)
        for outcome in outcomes
    ] == [
        (
            "endpoint_a_addend_a__endpoint_b_addend_b",
            "CCC[Si](C)(C)C",
        ),
        (
            "endpoint_b_addend_a__endpoint_a_addend_b",
            "CC(C)[Si](C)(C)C",
        ),
    ]


def test_rewrite_product_is_invariant_to_reactant_component_order() -> None:
    forward = _assignment(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "boron_transfer_coupling",
    )
    reversed_partners = _assignment(
        "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "boron_transfer_coupling",
    )

    forward_outcomes = apply_connectivity_rewrite(
        forward[1], forward[2], forward[0].reactants
    )
    reversed_outcomes = apply_connectivity_rewrite(
        reversed_partners[1],
        reversed_partners[2],
        reversed_partners[0].reactants,
    )

    assert forward_outcomes == reversed_outcomes


def test_each_ambiguous_beta_site_retains_operator_parity() -> None:
    parsed = parse_reaction_smiles("CCC(C)Br>>CC=CC")
    candidates = [
        (grammar, assignment)
        for grammar, assignment in enumerate_reaction_candidates(parsed.reactants)
        if grammar["id"] == "beta_halo_elimination"
    ]

    reports = [
        compare_connectivity_rewrite(grammar, assignment, parsed.reactants)
        for grammar, assignment in candidates
    ]

    assert len(reports) == 2
    assert all(report.parity for report in reports)
    assert {
        report.rewrite_outcomes[0].predicted_product_smiles
        for report in reports
    } == {"CC=CC", "C=CCC"}


def test_release_and_connect_supports_intramolecular_ring_closure() -> None:
    parsed, grammar, assignment = _assignment(
        "NCCCCBr>>C1CCCN1",
        "sp3_c_n_substitution",
    )

    report = compare_connectivity_rewrite(
        grammar, assignment, parsed.reactants
    )

    assert report.parity
    assert report.rewrite_outcomes[0].predicted_product_smiles == "C1CCNC1"


def test_release_and_connect_preserves_defined_alkene_stereo() -> None:
    parsed, grammar, assignment = _assignment(
        "Br/C=C/c1ccccc1.OB(O)c1ccccc1>>"
        "C(=C\\c1ccccc1)/c1ccccc1",
        "boron_transfer_coupling",
    )

    report = compare_connectivity_rewrite(
        grammar, assignment, parsed.reactants
    )

    assert report.parity
    assert report.rewrite_outcomes[0].predicted_product_smiles == (
        "C(=C\\c1ccccc1)/c1ccccc1"
    )


def test_executor_rejects_an_invalid_declared_before_state() -> None:
    parsed, grammar, assignment = _assignment(
        "CCBr.CN>>CCNC",
        "sp3_c_n_substitution",
    )
    electrophile = assignment["electrophile"]
    invalid_assignment = {
        **assignment,
        "electrophile": replace(
            electrophile,
            atom_roles={
                **electrophile.atom_roles,
                "handle": electrophile.atom_roles["anchor"],
            },
        ),
    }

    assert (
        apply_connectivity_rewrite(
            grammar, invalid_assignment, parsed.reactants
        )
        == ()
    )


def test_executor_does_not_guess_for_an_unregistered_grammar() -> None:
    parsed, _, assignment = _assignment(
        "CCBr.CN>>CCNC",
        "sp3_c_n_substitution",
    )

    assert (
        apply_connectivity_rewrite(
            {"id": "not_migrated"}, assignment, parsed.reactants
        )
        == ()
    )
