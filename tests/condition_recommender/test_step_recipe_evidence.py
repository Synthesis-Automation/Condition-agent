"""Regression tests for mixed source outcomes under a queried recipe."""

from copy import deepcopy
import pytest

from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.processed_catalog import ExperimentLookup
from condition_recommender.shared_core_index import build_shared_core_index, load_shared_core_index
from condition_recommender.step_plausibility import assess_step_plausibility
from condition_recommender.step_recipe_evidence import assess_step_recipe_evidence
from tests.condition_recommender.test_shared_core_retrieval import QUERY, record


def assay(number, *, product=None, kind="reported_product_yield_pct"):
    return {"outcome_type": kind, "unit": "percent",
            "value": number if isinstance(number, (float, int)) else None,
            "raw_value": str(number), "source_field": "assay",
            "metadata": {"measurement_basis": "NMR", **({"product_smiles": product} if product else {})}}


class Catalog:
    supports_experiment_analogues = True

    def __init__(self, rows):
        self.rows = rows

    def observation(self, oid, fields):
        for row in self.rows:
            if row["observation_id"] == oid:
                return row
        raise KeyError(oid)

    def experiment_candidates(self, smiles, *, limit, include_analogues, fragment_index):
        assert include_analogues
        return ExperimentLookup(tuple(self.rows[:limit]), len(self.rows) > limit,
                                min(limit, len(self.rows)), len(self.rows) <= limit, "exact_and_local")


def run_source(row, *outcomes, well="A1"):
    item = {**row, "source": {"experimental_procedure": "source procedure",
                             "raw_fields": {"source_observation": {
                                 "source": {"source_file_sha256": "123", "source_record_id": well},
                                 "raw_fields": {"source_location": f"Plate 1 well {well}"},
                                 "outcomes": list(outcomes)}}}}
    return item


def pair(tmp_path, rows):
    index = build_generic_index(rows)
    path = tmp_path / "shared.sqlite"
    build_shared_core_index(index, path)
    return index, load_shared_core_index(path, index)


def test_recorded_zero_and_trace_under_query_without_changing_grade(tmp_path):
    admitted = record(1, QUERY)
    indexes = pair(tmp_path, [admitted])
    zero = run_source({**record(2, QUERY), "index_eligibility": "review_only"}, assay(0), well="B1")
    zero["resolved_recipe"]["declared_absences"] = ["base"]
    trace = run_source({**zero, "observation_id": "obs-3"}, assay("Trace"), well="C1")
    reference = run_source(admitted, assay(70))
    query = {"declared_absences": ["base"], "temperature_c": 25}
    result = assess_step_recipe_evidence(QUERY, query, *indexes, experiment_catalog=Catalog([reference, trace, zero]),
                                         match_limit=1)
    assert result.step_assessment == assess_step_plausibility(QUERY, *indexes, match_limit=1)
    assert result.experiment_count == 3 and result.display_truncated
    assert result.experiments[0].outcome in ("zero_yield", "trace")
    assert result.experiments[0].source_locator["source_location"].startswith("Plate")
    groups = {c.status: c for c in result.outcome_counts}
    assert groups["matches_stated_fields"].zero_yield == groups["matches_stated_fields"].trace == 1
    assert groups["partially_reported"].at_or_above_threshold == 1
    assert result.experiments[0].index_eligibility == "review_only"


def test_target_assays_survive_zero_conversion_byproduct_and_conflicting_trace(tmp_path):
    src = record(1, QUERY)
    indexes = pair(tmp_path, [src])
    wrong = assay(0, product="CCO")
    conversion = assay(0, kind="conversion_pct")
    rows = [run_source(src, wrong, assay(41, product="N#Cc1ccccc1"), conversion),
            run_source({**src, "observation_id": "b"}, wrong, conversion, well="B1"),
            run_source({**src, "observation_id": "c"}, assay(0), assay("Trace"), well="C1")]
    result = assess_step_recipe_evidence(QUERY, {"temperature_c": 25}, *indexes, catalog=Catalog(rows))
    assert result.experiment_count == 3  # Indexed projection and exact catalog hit are one experiment.
    group = result.outcome_counts[0]
    assert group.at_or_above_threshold == group.unreported == group.unusable == 1 and group.zero_yield == 0
    by_id = {e.precedent.observation_id: e for e in result.experiments}
    assert by_id["obs-1"].precedent.yield_pct == 41 and len(by_id["obs-1"].reported_outcomes) == 3
    assert by_id["b"].precedent.yield_pct is None and by_id["c"].precedent.yield_pct is None


def test_reused_observation_ids_preserve_record_hashes_and_reordered_counts(tmp_path):
    row = run_source(record(1, QUERY), assay(70))
    distinct = deepcopy(row)
    distinct["source"]["raw_fields"]["source_observation"]["source"]["source_file_sha256"] = "789"
    indexes = pair(tmp_path, [row])
    query = {"temperature_c": 25}
    first = assess_step_recipe_evidence(QUERY, query, *indexes, catalog=Catalog([row]),
                                        experiment_catalog=Catalog([distinct, row]))
    reordered = assess_step_recipe_evidence(QUERY, query, *indexes, catalog=Catalog([row]),
                                            experiment_catalog=Catalog([row, distinct]))
    assert first == reordered and first.experiment_count == 2
    assert "OBSERVATION_ID_REUSED_FOR_DISTINCT_EXPERIMENT_RECORDS" in first.warnings
    assert len({e.record_sha256 for e in first.experiments}) == 2


def test_more_explicit_recipe_agreement_precedes_adverse_vaguely_reported_conditions(tmp_path):
    one = run_source(record(1, QUERY), assay(0))
    two = run_source(record(2, QUERY), assay(63), well="B1")
    two["resolved_recipe"]["time_h"] = 1
    indexes = pair(tmp_path, [one])
    query = {"temperature_c": 25, "time_h": 1, "illumination": "white LEDs"}
    result = assess_step_recipe_evidence(QUERY, query, *indexes, catalog=Catalog([one, two]), match_limit=1)
    assert result.outcome_counts[1].transferable_count == 2
    assert result.experiments[0].outcome == "at_or_above_threshold"


def test_unresolved_target_explains_that_catalogs_were_not_searched(tmp_path):
    indexes = pair(tmp_path, [record(1, QUERY)])
    result = assess_step_recipe_evidence("CC>>CCN", {"temperature_c": 25}, *indexes, catalog=Catalog([]))
    assert result.step_assessment.grade == "unresolved" and result.experiment_count == 0
    assert all(s.scope == "query_graph_unresolved" and not s.lookup_complete for s in result.search_scopes)
    assert "QUERY_GRAPH_UNRESOLVED_EXPERIMENT_SEARCH_NOT_PERFORMED" in result.warnings


def test_inconsistent_index_and_catalog_fields_cannot_create_hybrid_experiment(tmp_path):
    good = record(1, QUERY)
    indexes = pair(tmp_path, [good])
    stale = run_source(deepcopy(good), assay(0))
    stale["resolved_recipe"]["temperature_c"] = 100
    stale["yield_pct"] = 0
    with pytest.raises(ValueError, match="disagrees with indexed observation"):
        assess_step_recipe_evidence(QUERY, {"temperature_c": 25}, *indexes, catalog=Catalog([stale]))


def test_unattributed_scalar_cannot_override_product_specific_offtarget_yield(tmp_path):
    row = record(1, QUERY)
    indexes = pair(tmp_path, [row])
    observed = run_source(row, assay(98, kind="reported_yield_pct"), assay(98, product="CCO"))
    result = assess_step_recipe_evidence(QUERY, {"temperature_c": 25}, *indexes, catalog=Catalog([observed]))
    assert result.experiments[0].outcome == "unreported"
    assert len(result.experiments[0].reported_outcomes) == 2 and result.outcome_counts[0].at_or_above_threshold == 0


def test_recipe_yield_boundary_and_missing_measurements(tmp_path):
    source = record(1, QUERY)
    values = (9.99, None, 10, "Trace")
    rows = []
    for i, value in enumerate(values):
        numeric = isinstance(value, (int, float))
        observation = {**source, "observation_id": f"well-{i}", "yield_pct": value if numeric else None,
                       "outcome_status": "usable" if numeric else "unreported"}
        outcomes = () if value is None else (assay(value),)
        rows.append(run_source(observation, *outcomes, well=f"A{i + 1}"))
    indexes = pair(tmp_path, [rows[0]])
    result = assess_step_recipe_evidence(QUERY, {"temperature_c": 25}, *indexes, catalog=Catalog(rows))
    assert result.schema_version == "step_recipe_evidence.v2"
    counts = result.outcome_counts[0]
    assert counts.transferable_count == 4
    assert counts.low_yield == counts.at_or_above_threshold == counts.unreported == counts.trace == 1
    assert counts.zero_yield == counts.unusable == 0
    assert [e.outcome for e in result.experiments] == ["low_yield", "unreported", "at_or_above_threshold", "trace"]
    assert result.experiments[2].precedent.yield_pct == 10
