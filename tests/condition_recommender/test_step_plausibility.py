"""Regression cases for bounded, source-backed structural plausibility grades."""

from dataclasses import dataclass, replace
import json

import pytest

from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.shared_core_index import build_shared_core_index, load_shared_core_index
from condition_recommender.step_plausibility import assess_step_plausibility
import condition_recommender.step_plausibility as step_module
from tests.condition_recommender.test_shared_core_retrieval import QUERY, record


@dataclass(frozen=True)
class Source:
    reference: str
    yield_pct: float | None
    recipe: str = "recipe-A"
    status: str = "usable"
    reaction: str = QUERY
    observation: str | None = None


def evaluate(tmp_path, sources, *, query=QUERY, **limits):
    rows = []
    for n, source in enumerate(sources, 1):
        row = record(n, source.reaction)
        row.update(reference_id=source.reference, yield_pct=source.yield_pct,
                   outcome_status=source.status, resolved_recipe_id=source.recipe)
        row["resolved_recipe"]["recipe_id"] = source.recipe
        if source.observation:
            row["observation_id"] = source.observation
        rows.append(row)
    index = build_generic_index(rows)
    path = tmp_path / "core.sqlite"
    build_shared_core_index(index, path)
    shared = load_shared_core_index(path, index)
    return assess_step_plausibility(query, index, shared, **limits)


def test_display_limit_does_not_change_repeated_support(tmp_path):
    result = evaluate(tmp_path, [Source("paper-A", 61), Source("paper-B", 63, "recipe-B"),
                                 Source("paper-C", 84, "recipe-C")], match_limit=1)
    assert result.grade == "repeated_transferable_precedent"
    assert result.evidence_rank == 4
    assert result.qualified_count == 3 and result.matches_truncated
    assert result.support_counts.supporting_source_group_count == 3
    assert result.support_counts.supporting_distinct_reaction_count == 1
    assert len(result.matches) == 1


def test_same_reference_and_copied_experiment_do_not_inflate_depth(tmp_path):
    result = evaluate(tmp_path, [Source("paper-A", 61), Source("paper-A", 73),
                                 Source("translation-of-A", 73)])
    assert result.grade == "transferable_precedent"
    assert result.support_counts.supporting_source_group_count == 1
    assert result.support_counts.supporting_reference_count == 2


def test_zero_threshold_missing_and_invalid_outcomes_stay_separate(tmp_path):
    result = evaluate(tmp_path, [Source("A", 0), Source("B", 10), Source("E", 9.99),
                                 Source("C", None, status="missing"),
                                 Source("D", 98, status="invalid")])
    assert result.evidence_rank == 3
    assert result.support_counts.supporting_source_group_count == 1
    assert result.support_counts.supporting_observation_count == 3
    assert result.support_counts.at_or_above_threshold_yield_count == 1
    assert result.support_counts.zero_yield_count == 1
    assert result.support_counts.low_yield_count == 1
    assert result.support_counts.unknown_yield_count == 2
    assert result.measured_evidence_status == "mixed_at_or_above_threshold_and_low_or_zero"


@pytest.mark.parametrize(("number", "status", "expected", "outcome"), [
    (None, "missing", 3, "unreported"),
    (10, "usable", 3, "at_or_above_threshold"),
    (100, "usable", 3, "at_or_above_threshold"),
    (0, "usable", 2, "zero_yield"),
    (9.99, "usable", 2, "low_yield"),
    (2, "invalid", 3, "unusable"),
    (None, "invalid", 3, "unusable"),
])
def test_only_reliable_report_below_ten_excludes_comparable_source(tmp_path, number, status, expected, outcome):
    result = evaluate(tmp_path, [Source("A", number, status=status)])
    assert result.evidence_rank == expected
    match, = result.matches
    assert match.structurally_comparable
    assert match.outcome == outcome
    assert match.contributes_to_support == (expected == 3)
    assert match.excluded_by_low_yield == (expected == 2)
    assert result.support_counts.supporting_source_group_count == int(expected == 3)


@pytest.mark.parametrize("values", [(None, None, None, None), (10, 11.1, 75, 100), (60, None, 10, None)])
def test_missing_and_nonlow_yield_change_neither_copy_groups_nor_rank(tmp_path, values):
    sources = [
        Source(ref, n, recipe=recipe, status="usable" if n is not None else "missing")
        for n, (ref, recipe) in zip(values, (("A", "1"), ("copy-A", "1"), ("B", "2"), ("C", "3")))
    ]
    result = evaluate(tmp_path, sources)
    assert result.grade == "repeated_transferable_precedent"
    assert result.supporting_source_groups == (("A", "copy-A"), ("B",), ("C",))
    assert result.support_counts.supporting_observation_count == 4
    assert result.support_counts.at_or_above_threshold_yield_count == sum(n is not None for n in values)


def test_missing_recipe_cannot_bridge_distinct_documented_recipes(tmp_path):
    result = evaluate(tmp_path, [Source("A", None, "recipe-A", "missing"),
                                 Source("B", None, "recipe-B", "missing")])
    first = result.matches[0]
    incomplete = replace(first, precedent=replace(first.precedent, reference_id="C",
                                                  observation_id=None, resolved_recipe=None))
    groups = step_module._source_groups([*result.matches, incomplete], step_module._normalized_reaction)
    assert groups == {"A": "A", "B": "B", "C": "C"}


def test_copied_observation_with_measured_low_yield_cannot_regain_support(tmp_path):
    result = evaluate(tmp_path, [Source("A", 4, observation="identical-experiment"),
                                 Source("B", None, status="missing", observation="identical-experiment")])
    assert result.support_counts.supporting_source_group_count == 0 and result.evidence_rank == 2
    assert result.support_counts.unreported_yield_count == result.support_counts.low_yield_count == 1
    assert all(m.excluded_by_low_yield and not m.contributes_to_support for m in result.matches)
    assert "LOW_YIELD_REPORTED_FOR_SAME_OBSERVATION" in result.warnings


def test_product_without_matching_transformation_is_not_positive_support(tmp_path):
    result = evaluate(tmp_path, [Source("A", 92, reaction="NC(=O)c1ccccc1>>N#Cc1ccccc1")])
    assert result.grade == "no_qualified_precedent"
    assert result.evidence_rank == 0 and result.qualified_count == 0


def test_distant_substrates_do_not_earn_measured_support_by_repetition(tmp_path):
    distant = "Ic1ccc(CCCCCCCCCCCCCCCC)cc1.N#C[Cu]>>N#Cc1ccc(CCCCCCCCCCCCCCCC)cc1"
    result = evaluate(tmp_path, [Source(k, v, reaction=distant)
                                 for k, v in (("A", 61), ("B", 65), ("C", 78))])
    assert result.strongest_level == "observed_local"
    assert result.support_counts.comparable_observation_count == 0
    assert result.evidence_rank < 3


def test_additional_para_nitrile_on_every_source_limits_transfer(tmp_path):
    query_aryl = "Clc1cc(-c2ccc(OCC)cc2)ccc1"
    source_aryl = "Clc1cc(-c2ccc(OCC)cc2)c(C#N)cc1"

    def amination(aryl_chloride):
        return f"{aryl_chloride}.N1CCOCC1>>N1(CCOCC1){aryl_chloride[2:]}"

    query, source = amination(query_aryl), amination(source_aryl)
    result = evaluate(
        tmp_path, [Source(ref, y, reaction=source) for ref, y in (("A", 60), ("B", 75), ("C", 85))],
        query=query,
    )
    context = result.matches[0].context
    assert result.strongest_level == "observed_local"
    assert context.minimum_partner_similarity >= result.policy.minimum_partner_similarity
    assert context.status == "covered" and context.reverse_status == "missing"
    assert any("nitrile" in requirement for requirement in context.source_only)
    assert result.grade == "context_supported"
    assert result.support_counts.supporting_source_group_count == 0


def test_multiple_close_but_nonexact_examples_can_earn_repeated_support(tmp_path):
    query = "CCCc1ccc(I)cc1.N#C[Cu]>>CCCc1ccc(C#N)cc1"
    analogue = "CCCCc1ccc(I)cc1.N#C[Cu]>>CCCCc1ccc(C#N)cc1"
    homologue = "CCCCCc1ccc(I)cc1.N#C[Cu]>>CCCCCc1ccc(C#N)cc1"
    result = evaluate(tmp_path, [Source("A", None, status="missing", reaction=analogue),
                                 Source("B", 64, recipe="recipe-B", reaction=analogue),
                                 Source("C", 80, reaction=homologue)],
                      query=query)
    assert result.strongest_level == "observed_local"
    assert result.evidence_rank == 4 and result.support_counts.supporting_source_group_count == 3
    assert result.support_counts.supporting_distinct_reaction_count == 2


def test_highly_similar_opposite_stereoisomer_remains_analogue_evidence(tmp_path):
    query = "C[C@H](O)CCc1ccc(I)cc1.N#C[Cu]>>C[C@H](O)CCc1ccc(C#N)cc1"
    other = query.replace("[C@H]", "[C@@H]")
    result = evaluate(tmp_path, [Source("A", 90, reaction=other)], query=query)
    assert result.strongest_level == "observed_local"
    assert result.matches[0].context.minimum_partner_similarity > 0.5
    assert result.support_counts.supporting_source_group_count == 0
    assert "ANALOGUE_STEREOCHEMICAL_TRANSFER_UNVERIFIED" in result.warnings


def test_anonymous_measured_observation_cannot_create_attributed_support(tmp_path):
    result = evaluate(tmp_path, [Source("", 88)])
    assert result.support_counts.at_or_above_threshold_yield_count == 1
    assert result.support_counts.missing_reference_supporting_count == 1
    assert result.support_counts.supporting_source_group_count == 0
    assert result.evidence_rank < 3


def test_displays_a_measured_counterexample_beside_a_positive(tmp_path):
    result = evaluate(tmp_path, [Source("A", 70), Source("B", 71), Source("C", 0)],
                      match_limit=2)
    assert result.support_counts.zero_yield_count == 1
    assert result.measured_evidence_status == "mixed_at_or_above_threshold_and_low_or_zero"
    assert {m.outcome for m in result.matches} == {"at_or_above_threshold", "zero_yield"}
    assert "COMPARABLE_LOW_OR_ZERO_OUTCOMES_INSPECT_SOURCE_CONDITIONS" in result.warnings


def test_unknown_graph_stays_unresolved_even_with_a_successful_source(tmp_path):
    result = evaluate(tmp_path, [Source("A", 95)], query="CC>>CC")
    assert result.grade == "unresolved"
    assert result.evidence_rank is None


def test_reordered_query_gives_same_grade_and_support_counts(tmp_path):
    sources = [Source("A", 61), Source("B", 65)]
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    first = evaluate(a, sources)
    second = evaluate(b, sources, query="[Cu]C#N.c1ccc(I)cc1>>c1ccc(C#N)cc1")
    assert first.grade == second.grade
    assert first.support_counts == second.support_counts


@pytest.mark.parametrize("limit", [{"candidate_limit": False}, {"match_limit": 0}])
def test_limits_reject_invalid_inputs(tmp_path, limit):
    with pytest.raises(ValueError):
        evaluate(tmp_path, [Source("A", 70)], **limit)


@pytest.mark.parametrize(("field", "invalid_value"), [
    ("repeated_support_minimum_source_groups", 1),
    ("minimum_partner_similarity", 1.1),
    ("low_yield_below_pct", 150),
    ("transferable_levels", ["retained_typed"]),
])
def test_invalid_definition_rejected(tmp_path, monkeypatch, field, invalid_value):
    definition = json.loads(step_module._POLICY_PATH.read_text(encoding="utf-8"))
    definition[field] = invalid_value
    invalid = tmp_path / "invalid.json"
    invalid.write_text(json.dumps(definition), encoding="utf-8")
    monkeypatch.setattr(step_module, "_POLICY_PATH", invalid)
    with pytest.raises(ValueError, match="policy"):
        step_module.load_step_plausibility_policy()
