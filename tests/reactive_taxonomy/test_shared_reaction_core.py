"""Shared graph projection positives, close negatives and evidence conflicts."""

from copy import deepcopy
from dataclasses import asdict, replace

import pytest
from rdkit import Chem

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.shared_reaction_core import (
    LEVELS,
    SharedReactionCore,
    build_shared_reaction_core,
    compare_reaction_cores,
    load_shared_core_rules,
)

QUERY = "Ic1ccccc1.N#C[Cu]>>N#Cc1ccccc1"


def project(smiles):
    analysis = featurize_reaction(smiles)
    return build_shared_reaction_core(
        smiles,
        asdict(analysis.reaction_signature) if analysis.reaction_signature else {},
        asdict(analysis.reaction_core) if analysis.reaction_core else {},
    )


@pytest.mark.parametrize("leaving", ["Br", "Cl", "CS(=O)(=O)O", "OS(=O)(=O)O"])
def test_observed_leaving_fragments_generalize_without_virtual_reactions(leaving):
    query = project(QUERY)
    precedent = project(f"{leaving}c1ccccc1.N#C[Cu]>>N#Cc1ccccc1")
    comparison = compare_reaction_cores(query, precedent)
    assert comparison.eligible and comparison.level == "retained_local"
    assert comparison.relation == "analogue_evidence"
    assert query.levels[0].key != precedent.levels[0].key
    assert query.realization_key != precedent.realization_key


def test_hcn_and_copper_source_are_related_but_remain_literal_different_inputs():
    query, precedent = project(QUERY), project("Ic1ccccc1.C#N>>N#Cc1ccccc1")
    assert compare_reaction_cores(query, precedent).level == "retained_local"
    assert query.input_identity != precedent.input_identity
    assert "[H]" in precedent.realization_details
    assert "[Cu]" in query.realization_details


def test_different_aromatic_environment_uses_typed_level():
    query = project(QUERY)
    precedent = project("Brc1ncccc1.N#C[Cu]>>N#Cc1ncccc1")
    comparison = compare_reaction_cores(query, precedent)
    assert comparison.eligible and comparison.level == "retained_typed"


@pytest.mark.parametrize(
    "partner,product",
    [
        ("CN", "CNc1ccccc1"),
        ("CO", "COc1ccccc1"),
        ("CS", "CSc1ccccc1"),
        ("OB(O)c1ccccc1", "c1ccc(-c2ccccc2)cc1"),
    ],
)
def test_same_contract_supports_cn_co_cs_and_suzuki(partner, product):
    query = project(f"Ic1ccccc1.{partner}>>{product}")
    precedent = project(f"Brc1ccccc1.{partner}>>{product}")
    assert compare_reaction_cores(query, precedent).level == "retained_local"


@pytest.mark.parametrize(
    "reaction",
    [
        "Brc1ccccc1.CN>>CNc1ccccc1",
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "N#Cc1ccccc1>>NCc1ccccc1",
        "c1ccccc1.C#N>>N#Cc1ccccc1",
        "NC=O>>C#N",
    ],
)
def test_unrelated_joins_reductions_and_carbon_hydrogen_edits_do_not_match(reaction):
    comparison = compare_reaction_cores(project(QUERY), project(reaction))
    assert not comparison.eligible


def test_partner_atom_order_and_map_numbers_do_not_change_keys():
    baseline = project(QUERY)
    reordered = project("[Cu]C#N.c1ccc(I)cc1>>c1ccc(C#N)cc1")
    assert [level.key for level in baseline.levels] == [
        level.key for level in reordered.levels
    ]
    mapped = "[I:7][c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[N:8]#[C:9][Cu:10]>>[N:8]#[C:9][c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1"
    remapped = mapped
    for n in range(10, 0, -1):
        remapped = remapped.replace(f":{n}]", f":{n + 20}]")
    assert project(mapped).levels == project(remapped).levels


def test_aromatic_kekule_serialization_is_invariant():
    parts = QUERY.split(">")
    kekule = ">".join(
        Chem.MolToSmiles(Chem.MolFromSmiles(part), kekuleSmiles=True) if part else ""
        for part in parts
    )
    assert [p.key for p in project(QUERY).levels] == [
        p.key for p in project(kekule).levels
    ]


@pytest.mark.parametrize(
    "mutation", ["blocked", "hypothesis", "conflict", "wrong_atom"]
)
def test_ambiguous_conflicting_and_contradictory_observations_abstain(mutation):
    analysis = featurize_reaction(QUERY)
    signature, core = (
        asdict(analysis.reaction_signature),
        deepcopy(asdict(analysis.reaction_core)),
    )
    if mutation == "blocked":
        core["quality"]["status"] = "blocked"
    elif mutation == "hypothesis":
        core["evidence_status"] = "hypothesis"
    elif mutation == "conflict":
        core["warnings"] = ("MAPPED_OPERATOR_CONFLICT",)
    else:
        core["atom_transitions"][0]["before_state"]["element"] = "Si"
    projection = build_shared_reaction_core(QUERY, signature, core)
    assert not projection.levels and projection.unavailable_reasons


def test_order_changes_receive_broad_views_but_conflicting_event_counts_abstain():
    for reaction in ("CC=O>>CCO", "CCBr.CCO>>CCOCC"):
        projection = project(reaction)
        assert tuple(p.level for p in projection.levels) == LEVELS
    analysis = featurize_reaction(QUERY)
    signature, core = (
        asdict(analysis.reaction_signature),
        asdict(analysis.reaction_core),
    )
    signature["event_count"] = 2
    projection = build_shared_reaction_core(QUERY, signature, core)
    assert not projection.levels
    assert projection.unavailable_reasons == ("EVENT_COUNT_CONTRADICTS_OBSERVATION",)


def test_schema_roundtrip_and_version_mismatch():
    projection = project(QUERY)
    assert SharedReactionCore.from_dict(projection.to_dict()) == projection
    assert not compare_reaction_cores(
        projection, replace(projection, definition_hash="other")
    ).eligible
    with pytest.raises(ValueError, match="incompatible"):
        SharedReactionCore.from_dict({**projection.to_dict(), "schema_version": "0"})
    assert tuple(load_shared_core_rules()["levels"]) == LEVELS


def test_noop_and_invalid_maps_never_provide_a_core():
    assert not project("CC>>CC").levels
    assert not project("[CH3:1][Br:1]>>[CH4:1]").levels
