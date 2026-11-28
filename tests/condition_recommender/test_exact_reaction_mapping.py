"""Regression tests for source-backed recovery on the benzoxazole acid example."""

import re

import pytest

from reactive_taxonomy.chemistry.source_mapping import source_material_graphs
from condition_recommender.exact_reaction_mapping import recover_exact_reaction_mapping
from condition_recommender.step_plausibility import assess_step_plausibility
from tests.condition_recommender.test_shared_core_retrieval import engine, record


QUERY = "Nc1cccc(C(=O)O)c1O.CCOC(C)(OCC)OCC>>Cc1nc2cccc(C(=O)O)c2o1"
LEFT = ("[C:5]([CH3:6])([O:16][CH2:15][CH3:14])([O:17][CH2:18][CH3:19])[O:20][CH2:21][CH3:22]."
        "[NH2:1][c:2]1[c:3]([OH:4])[c:7]([C:8]([OH:9])=[O:10])[cH:11][cH:12][cH:13]1")
RIGHT = "[n:1]1[c:2]2[c:3]([o:4][c:5]1[CH3:6])[c:7]([C:8]([OH:9])=[O:10])[cH:11][cH:12][cH:13]2"
MAPPED = f"{LEFT}>>{RIGHT}"


@pytest.fixture(scope="module")
def indexes(tmp_path_factory):
    renumbered = re.sub(r":(\d+)\]", lambda m: f":{int(m[1]) + 100}]", MAPPED)
    rows = [record(1, MAPPED), record(2, renumbered)]
    for row in rows:
        row.update(yield_pct=None, outcome_status="missing")
    recommender = engine(tmp_path_factory.mktemp("source-maps"), rows, sqlite=True)
    return recommender.index, recommender.shared_core_index


def test_unresolved_query_recovers_documented_mapping_but_no_numeric_yield(indexes):
    assessment = assess_step_plausibility(QUERY, *indexes, match_limit=1)
    assert assessment.query_reaction_smiles == QUERY
    assert "QUERY_GRAPH_EVIDENCE_UNAVAILABLE" in assessment.original_query_warnings
    assert assessment.grade == "transferable_precedent"
    assert assessment.support_counts.supporting_distinct_reaction_count == 1
    assert assessment.support_counts.unknown_yield_count == 2
    assert assessment.measured_evidence_status == "no_comparable_numeric_outcome"
    recovered = assessment.mapping_recovery
    assert recovered.status == "recovered"
    assert {s.observation_id for s in recovered.sources} == {"obs-1", "obs-2"}
    assert len({s.mapping_equivalence_key for s in recovered.sources}) == 1
    assert recover_exact_reaction_mapping(QUERY.replace(">>", ">O>"), *indexes).mapped_reaction_smiles == f"{LEFT}>O>{RIGHT}"


@pytest.mark.parametrize("query", [QUERY.replace(".CCOC", ".COC"), QUERY.replace(".CCOC", ".O.CCOC"),
                                   QUERY.replace(">>C", ">>[13CH3]"), QUERY + ".O"])
def test_identical_product_or_similar_reactants_do_not_supply_a_mapping(indexes, query):
    recovered = recover_exact_reaction_mapping(query, *indexes)
    assert recovered.status == "not_found" and recovered.mapped_reaction_smiles is None


def test_opposite_or_unspecified_stereochemistry_is_a_different_full_graph(tmp_path):
    stereo = "[C@H:6]([F:23])[Cl:24]"
    left, right = (side.replace("[CH3:6]", stereo) for side in (LEFT, RIGHT))
    recommender = engine(tmp_path, [record(1, f"{left}>>{right}")])
    inputs = ".".join(source_material_graphs(left))
    product, = source_material_graphs(right)
    alternatives = (product.replace("@@", "@") if "@@" in product else product.replace("@", "@@"),
                    product.replace("@", ""))
    for different in alternatives:
        result = recover_exact_reaction_mapping(f"{inputs}>>{different}", recommender.index,
                                                recommender.shared_core_index)
        assert result.status == "not_found" and not result.mapped_reaction_smiles


def test_conflicting_recorded_edit_graphs_leave_step_unresolved(tmp_path):
    other = RIGHT.replace(":2]", ":900]").replace(":7]", ":2]").replace(":900]", ":7]")
    recommender = engine(tmp_path, [record(1, MAPPED), record(2, f"{LEFT}>>{other}")])
    assessment = assess_step_plausibility(QUERY, recommender.index, recommender.shared_core_index)
    assert assessment.grade == "unresolved"
    assert assessment.mapping_recovery.status == "conflicting_mappings"
    assert assessment.mapping_recovery.mapped_reaction_smiles is None
    assert len({s.mapping_equivalence_key for s in assessment.mapping_recovery.sources}) == 2


def test_capped_and_supplied_mappings_are_never_replaced(indexes):
    limited = assess_step_plausibility(QUERY, *indexes, candidate_limit=1)
    assert limited.grade == "unresolved"
    assert limited.mapping_recovery.status == "truncated"
    assert limited.mapping_recovery.candidates_examined == 1
    first, = limited.mapping_recovery.sources
    assert first.observation_id == "obs-1" and first.mapping_equivalence_key is None
    damaged = MAPPED.replace("[n:1]1", "[n:999]1")
    recovered = recover_exact_reaction_mapping(damaged, *indexes)
    assert recovered.status == "supplied_mapping_preserved" and not recovered.mapped_reaction_smiles
