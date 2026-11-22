"""Graph-derived partner assignment and conserved precursor context regressions."""

from dataclasses import replace

import pytest

from reactive_taxonomy.step_precedent_context import (
    build_step_context,
    compare_step_contexts,
)


SALICYLATE = "O=C(O)c1ccc(Br)cc1O.CN>>O=C(O)c1ccc(NC)cc1O"
REORDERED = "NC.Oc1cc(Br)ccc1C(=O)O>>CNc1ccc(C(=O)O)c(O)c1"
META = "O=C(O)c1cc(O)ccc1Br.CN>>O=C(O)c1cc(O)ccc1NC"
SPLIT = "Brc1cc(O)ccc1.NCC(=O)O>>O=C(O)CNc1cc(O)ccc1"
PROTECTED = "CC(=O)Oc1cc(Br)ccc1C(=O)OC.CN>>CC(=O)Oc1cc(NC)ccc1C(=O)OC"


def features(requirements: tuple[str, ...]) -> set[str]:
    return {r.split(" / ")[1] for r in requirements if " to " not in r}


def test_unchanged_query_order_independent_of_atom_indices_or_serialization():
    query = build_step_context(SALICYLATE)
    reordered = build_step_context(REORDERED)
    assert query == reordered and query.status == "resolved"
    assert features(query.requirements) == {
        "phenol",
        "carboxylic_acid",
        "ortho_free_hydroxybenzoic_acid_or_ion",
    }
    result = compare_step_contexts(query, reordered)
    assert result.status == "covered" and result.covered == query.requirements
    assert result.minimum_partner_similarity == 1


@pytest.mark.parametrize(
    ("record", "missing"),
    [
        (META, {"ortho_free_hydroxybenzoic_acid_or_ion"}),
        (SPLIT, {"carboxylic_acid", "ortho_free_hydroxybenzoic_acid_or_ion"}),
        (
            PROTECTED,
            {"phenol", "carboxylic_acid", "ortho_free_hydroxybenzoic_acid_or_ion"},
        ),
    ],
)
def test_anchor_placement_and_protection_mismatches_are_explained(record, missing):
    result = compare_step_contexts(
        build_step_context(SALICYLATE), build_step_context(record)
    )
    assert result.status == "missing"
    assert features(result.missing) == missing
    assert 0 < result.minimum_partner_similarity < 1


def test_ions_are_distinct_and_do_not_become_esters():
    ion = build_step_context("O=C([O-])c1ccc(Br)cc1[O-].CN>>O=C([O-])c1ccc(NC)cc1[O-]")
    assert features(ion.requirements) == {
        "phenoxide",
        "carboxylate",
        "ortho_free_hydroxybenzoic_acid_or_ion",
    }
    compared = compare_step_contexts(build_step_context(SALICYLATE), ion)
    assert features(compared.missing) == {"carboxylic_acid", "phenol"}


def test_distance_to_reacting_atom_distinguishes_three_nitrile_regioisomers():
    contexts = [
        build_step_context(f"Br{ring}.CN>>CN{ring}")
        for ring in ("c1ccccc1C#N", "c1cccc(C#N)c1", "c1ccc(C#N)cc1")
    ]
    for ctx, distance in zip(contexts, ("= 2 bonds", "= 3 bonds", "> 3 bonds")):
        assert any(f"nitrile to edit {distance}" in r for r in ctx.requirements)
        assert compare_step_contexts(ctx, ctx).status == "covered"
        assert all(
            compare_step_contexts(ctx, s).status == "missing"
            for s in contexts
            if s != ctx
        )


def test_missing_pairwise_proximity_despite_both_groups_far_from_the_reacting_center():
    near, far = [
        build_step_context(f"Br{ring}.CN>>CN{ring}")
        for ring in ("c1ccc(CC(O)C#N)cc1", "c1ccc(CC(O)CCCC#N)cc1")
    ]
    compared = compare_step_contexts(near, far)
    assert compared.status == "missing"
    assert [r.split(" / ")[1] for r in compared.missing] == [
        "alcohol to nitrile = 1 bond"
    ]
    assert features(compared.covered) == {"alcohol", "nitrile"}
    reversed_result = compare_step_contexts(far, near)
    assert (
        reversed_result.status == "covered"
        and reversed_result.reverse_status == "missing"
    )
    assert reversed_result.source_only == compared.missing


def test_source_only_nitrile_on_aryl_chloride_cannot_supply_bidirectional_support():
    rings = ("c1ccccc1", "c1ccc(C#N)cc1")
    plain, cyano = [build_step_context(f"Cl{ring}.CN>>CN{ring}") for ring in rings]
    result = compare_step_contexts(plain, cyano)
    assert result.status == "covered" and result.reverse_status == "missing"
    assert features(result.source_only) == {"nitrile"}
    assert any("nitrile to edit > 3 bonds" in r for r in result.source_only)
    assert "SOURCE_ONLY_RETAINED_CONTEXT" in result.warnings


@pytest.mark.parametrize(
    ("precursor", "product", "feature"),
    [
        ("Brc1ccc(Br)cc1", "CNc1ccc(Br)cc1", "aryl_halide(Br)"),
        ("Brc1ccc(CBr)cc1", "CNc1ccc(CBr)cc1", "LG|Alkyl|Br"),
        ("Brc1ccc(B(O)O)cc1", "CNc1ccc(B(O)O)cc1", "TM|Ar|B(OH)2"),
    ],
)
def test_only_the_unconsumed_handle_is_required(precursor, product, feature):
    result = build_step_context(f"{precursor}.CN>>{product}")
    assert features(result.requirements) == {feature}


def test_empty_registry_requirements_and_unresolved_graphs_are_distinct():
    empty = build_step_context("BrCC.CN>>CCNC")
    result = compare_step_contexts(empty, empty)
    assert result.status == "covered" and not result.covered
    assert "NO_TRACKED_SPECTATOR_REQUIREMENTS" in result.warnings
    invalid = build_step_context("[CH3:1][Br:2].[NH3:1]>>[CH3:1][NH2:1]")
    assert invalid.status == "unresolved"
    assert compare_step_contexts(empty, invalid).status == "unresolved"


def test_duplicate_anchor_assignments_do_not_pool_evidence():
    a = "[Br:1][c:2]1[cH:3][cH:4][cH:5][cH:6][c:7]1[OH:8]"
    b = "[Br:9][c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1"
    p = "[OH:8][c:7]1[cH:6][cH:5][cH:4][cH:3][c:2]1-[c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1"
    query = build_step_context(f"{a}.{b}>>{p}")
    changed = build_step_context(f"{a}.{b}>>{p}".replace("[OH:8]", "[F:8]"))
    assert query.status == changed.status == "resolved"
    assert compare_step_contexts(query, query).status == "covered"
    result = compare_step_contexts(query, changed)
    assert result.status == "unresolved" and result.minimum_partner_similarity is None
    assert result.reverse_status == "unresolved"
    assert "AMBIGUOUS_PARTNER_ALIGNMENT" in result.warnings


def test_unknown_retention_remains_unresolved_even_with_identical_precursor():
    query = build_step_context(SALICYLATE)
    source_partners = tuple(
        replace(p, uncertain_features=p.retained_features, retained_features=())
        for p in query.partners
    )
    source = replace(query, partners=source_partners)
    result = compare_step_contexts(query, source)
    assert result.status == "unresolved"
    assert result.unresolved == query.requirements and not result.missing
    assert (
        result.reverse_status == "unresolved"
        and result.source_unresolved == source.requirements
    )


@pytest.mark.parametrize(
    "rxn",
    [
        "BrCC.N[C@H](C)CC>>CCN[C@H](C)CC",
        "BrCC.NCC/C=C/C>>CCNCC/C=C/C",
    ],
)
def test_stereo_sensitive_queries_are_flagged_for_separate_transfer_review(rxn):
    assert build_step_context(rxn).has_specified_stereo
