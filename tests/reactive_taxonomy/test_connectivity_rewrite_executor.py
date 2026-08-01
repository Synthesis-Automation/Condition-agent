"""V3 contract tests for interpretation-independent connectivity rewrites."""

from dataclasses import replace

import pytest

from reactive_taxonomy import (
    SITE_INTERFACE_SCHEMA_VERSION,
    apply_reaction_operator,
    compile_connectivity_rewrite_definitions,
    featurize_molecule,
    featurize_reaction,
    load_connectivity_rewrites,
    reaction_signature_definition_versions,
)
from reactive_taxonomy.connectivity_rewrite import (
    CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION,
    CONNECTIVITY_REWRITE_SCHEMA_VERSION,
)
from reactive_taxonomy.reaction_parser import parse_reaction_smiles
from reactive_taxonomy.reaction_reconstruction import (
    enumerate_reconstruction_candidates,
)
from reactive_taxonomy.reaction_reconstruction_rules import (
    load_reaction_reconstruction_rules,
)


def _assignment(reaction_smiles: str, rule_id: str):
    parsed = parse_reaction_smiles(reaction_smiles)
    matches = [
        (rule, assignment)
        for rule, assignment in enumerate_reconstruction_candidates(parsed.reactants)
        if rule["id"] == rule_id
    ]
    assert matches
    return parsed, matches[0][0], matches[0][1]


def _apply(rule, assignment, components):
    bindings = rule["operator_slot_bindings"]
    operator_assignment = {
        operator_slot: assignment[rule_slot]
        for operator_slot, rule_slot in bindings.items()
    }
    return apply_reaction_operator(
        rule["operator_id"],
        operator_assignment,
        components,
        output_role_labels={
            operator_slot: rule_slot
            for operator_slot, rule_slot in bindings.items()
        },
    )


MIGRATION_CORPUS = (
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "Brc1ccccc1.Cl[Zn]c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    "Brc1ccccc1.CN>>CNc1ccccc1",
    "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
    "Brc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1",
    "CCBr.CN>>CCNC",
    "CCCCBr.Oc1ccccc1>>CCCCOc1ccccc1",
    "CI.SC>>CSC",
    "Brc1ccccc1.C#Cc1ccccc1>>c1ccc(C#Cc2ccccc2)cc1",
    "Brc1ccccc1.CC#N>>N#CCc1ccccc1",
    "CC(=O)Cl.CN>>CC(=O)NC",
    "CC(=O)Cl.CCO>>CCOC(C)=O",
    "CC(=O)Cl.CCS>>CCSC(C)=O",
    "CC(=O)Cl.c1ccccc1>>CC(=O)c1ccccc1",
    "CS(=O)(=O)Cl.CN>>CS(=O)(=O)NC",
    "Oc1ccccc1.O=S(=O)(F)CC>>O=S(=O)(CC)Oc1ccccc1",
    "CC(=O)[S-].[K+].Ic1ccccc1>>CC(=O)Sc1ccccc1",
    "C=C.BrBr>>BrCCBr",
    "CC=C.CN>>CCCNC",
    "C#C.BrBr>>BrC=CBr",
    "CC#C.CN>>CC=CNC",
    "CCC(C)Br>>CC=CC",
    "CC(=O)C>>CC(O)C",
    "CC(O)C>>CC(=O)C",
    "CC=CC>>CCCC",
    "CC#CC>>CC=CC",
    "CC#CC>>CCCC",
)


def test_v3_definitions_separate_rules_from_graph_operators() -> None:
    rules = load_reaction_reconstruction_rules()
    rewrites = load_connectivity_rewrites()
    rewrite_ids = {rewrite.rewrite_id for rewrite in rewrites}

    assert {str(rule["operator_id"]) for rule in rules} <= rewrite_ids
    assert all(not hasattr(rewrite, "annotation_ids") for rewrite in rewrites)
    assert {
        rewrite.site_interface_schema_version for rewrite in rewrites
    } == {SITE_INTERFACE_SCHEMA_VERSION}
    assert SITE_INTERFACE_SCHEMA_VERSION == "2.0"

    versions = reaction_signature_definition_versions()
    assert {
        "site_patterns.v2.json",
        "context_facets.v2.json",
        "site_interfaces.v2.json",
        "reaction_reconstruction_rules.v1.json",
        "connectivity_rewrites.v3.json",
        "reactivity_descriptor_rules.v1.json",
        "aromatic_systems.v1.json",
        "signature_features.v3.json",
        "taxonomy_manifest.v3.json",
    } <= set(versions)
    assert all("@sha256:" in value for key, value in versions.items() if key.endswith(".json"))


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
        {
            "op": "change_localized_bond_state",
            "endpoints": [
                "left.reactive_link.unknown",
                "right.connection_endpoint.atom",
            ],
            "before": "NONE",
            "after": "SINGLE",
        },
    ],
)
def test_compiler_rejects_unbounded_or_invalid_instructions(instruction) -> None:
    payload = {
        "schema_version": CONNECTIVITY_REWRITE_SCHEMA_VERSION,
        "instruction_set_version": CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION,
        "site_interface_schema_version": SITE_INTERFACE_SCHEMA_VERSION,
        "rewrites": [
            {
                "id": "invalid",
                "template": "release_and_connect",
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
        "schema_version": CONNECTIVITY_REWRITE_SCHEMA_VERSION,
        "instruction_set_version": CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION,
        "site_interface_schema_version": SITE_INTERFACE_SCHEMA_VERSION,
        "rewrites": [
            {
                "id": "unsafe_projection",
                "template": "release_and_connect",
                "variants": [
                    {
                        "id": "default",
                        "instructions": [
                            {
                                "op": "declare_projection_discardable_attachment",
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


def test_compiler_rejects_unknown_tetrahedral_outcome() -> None:
    payload = {
        "schema_version": CONNECTIVITY_REWRITE_SCHEMA_VERSION,
        "instruction_set_version": CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION,
        "site_interface_schema_version": SITE_INTERFACE_SCHEMA_VERSION,
        "rewrites": [
            {
                "id": "invalid_stereo",
                "template": "release_and_connect",
                "variants": [
                    {
                        "id": "default",
                        "instructions": [
                            {
                                "op": "set_tetrahedral_outcome",
                                "selector": "left.anchor",
                                "outcome": "invert_unconditionally",
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

    with pytest.raises(ValueError, match="Invalid tetrahedral outcome"):
        compile_connectivity_rewrite_definitions(payload)


def test_compiler_rejects_interpretation_metadata_in_graph_operators() -> None:
    payload = {
        "schema_version": CONNECTIVITY_REWRITE_SCHEMA_VERSION,
        "instruction_set_version": CONNECTIVITY_REWRITE_INSTRUCTION_SET_VERSION,
        "site_interface_schema_version": SITE_INTERFACE_SCHEMA_VERSION,
        "rewrites": [
            {
                "id": "invalid_interpretation_coupling",
                "template": "release_and_connect",
                "annotation_ids": ["named_reaction"],
                "variants": [],
            }
        ],
    }

    with pytest.raises(ValueError, match="cannot contain interpretation metadata"):
        compile_connectivity_rewrite_definitions(payload)


@pytest.mark.parametrize("reaction_smiles", MIGRATION_CORPUS)
def test_v2_rewrite_corpus_reconstructs_the_observed_product(
    reaction_smiles: str,
) -> None:
    analysis = featurize_reaction(reaction_smiles)

    assert analysis.valid
    assert analysis.evidence_quality in {
        "exact_product_reconstruction",
        "exact_multi_event_reconstruction",
        "validated_mapping",
        "validated_mapping_and_reconstruction",
    }
    assert analysis.reaction_signature is not None
    assert analysis.reaction_signature.schema_version == "3.2"
    assert analysis.reaction_signature.signature_id.startswith("RS3:")


@pytest.mark.parametrize(
    ("reaction_smiles", "expected_elements"),
    [
        ("C=C.BrBr>>BrCCBr", {"Br"}),
        ("C=C.ICl>>ICCCl", {"I", "Cl"}),
        ("C=C.C[SiH](C)C>>CC[Si](C)(C)C", {"Si", "H"}),
        ("C=C.B[Si](C)(C)C>>BCC[Si](C)(C)C", {"Si", "B"}),
    ],
)
def test_general_ab_additions_share_one_edit_derived_archetype(
    reaction_smiles: str,
    expected_elements: set[str],
) -> None:
    analysis = featurize_reaction(reaction_smiles)

    assert analysis.edit_archetype == "addition"
    assert analysis.selected_candidate is not None
    assert analysis.selected_candidate.annotation_id == "addend_pair_addition_to_alkene"
    edits = analysis.reaction_signature.edits
    observed_elements = {
        atom.element
        for edit in edits
        for atom in (edit.atom_1, edit.atom_2)
        if atom is not None
    }
    assert expected_elements <= observed_elements | {"H"}


def test_si_b_is_exposed_as_an_explicit_reactive_link() -> None:
    analysis = featurize_molecule("B[Si](C)(C)C")
    links = [
        link
        for site in analysis.connectivity_sites
        for link in site.reactive_links
        if {link.endpoint_a.element, link.endpoint_b.element} == {"Si", "B"}
    ]

    assert links
    assert all(link.source_kind == "explicit_bond" for link in links)
    assert all(link.before_order == "SINGLE" for link in links)


def test_pair_addition_enumerates_and_orders_both_regioisomers() -> None:
    parsed, interpretation, assignment = _assignment(
        "CC=C.C[SiH](C)C>>CCC[Si](C)(C)C",
        "addend_pair_addition_to_alkene",
    )

    outcomes = _apply(interpretation, assignment, parsed.reactants)

    assert [
        (outcome.outcome_id, outcome.predicted_product_smiles)
        for outcome in outcomes
    ] == [
        ("endpoint_a_addend_a__endpoint_b_addend_b", "CCC[Si](C)(C)C"),
        ("endpoint_b_addend_a__endpoint_a_addend_b", "CC(C)[Si](C)(C)C"),
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

    assert _apply(
        forward[1], forward[2], forward[0].reactants
    ) == _apply(
        reversed_partners[1],
        reversed_partners[2],
        reversed_partners[0].reactants,
    )


def test_release_and_connect_supports_intramolecular_ring_closure() -> None:
    parsed, interpretation, assignment = _assignment(
        "NCCCCBr>>C1CCCN1",
        "sp3_c_n_substitution",
    )

    outcomes = _apply(interpretation, assignment, parsed.reactants)

    assert outcomes[0].predicted_product_smiles == "C1CCNC1"


def test_executor_rejects_an_invalid_declared_before_state() -> None:
    parsed, interpretation, assignment = _assignment(
        "CCBr.CN>>CCNC",
        "sp3_c_n_substitution",
    )
    leaving_slot = next(
        slot
        for slot, constraint in interpretation["slots"].items()
        if constraint["site_type"] == "leaving_group"
    )
    electrophile = assignment[leaving_slot]
    invalid_assignment = {
        **assignment,
        leaving_slot: replace(
            electrophile,
            atom_roles={
                **electrophile.atom_roles,
                "handle": electrophile.atom_roles["anchor"],
            },
        ),
    }

    assert _apply(interpretation, invalid_assignment, parsed.reactants) == ()


def test_executor_does_not_guess_for_an_unregistered_operator() -> None:
    parsed, _, assignment = _assignment(
        "CCBr.CN>>CCNC",
        "sp3_c_n_substitution",
    )

    assert apply_reaction_operator(
        "not_registered", assignment, parsed.reactants
    ) == ()
