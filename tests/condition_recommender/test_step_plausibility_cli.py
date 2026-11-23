"""Exercise standalone CLI validation and both output forms on bound SQLite files."""

import json
from pathlib import Path

import pytest

from condition_recommender.experiment_catalog import build_experiment_catalog
from condition_recommender.generic_indexing import build_generic_index, load_generic_index
from condition_recommender.shared_core_index import build_shared_core_index
from condition_recommender.sqlite_indexing import save_sqlite_generic_index
from condition_recommender.step_plausibility_cli import main
from tests.condition_recommender.test_shared_core_retrieval import QUERY, record


@pytest.mark.parametrize(("percent", "status", "printed"), [
    (70, "usable", "source yield 70%"),
    (98, "invalid", "yield unusable (raw 98"),
    (None, "missing", "yield unreported"),
])
def test_explicit_bound_files_supply_summary_and_json(tmp_path, capsys, percent, status, printed):
    canonical = tmp_path / "index.sqlite"
    shared = tmp_path / "core.sqlite"
    source_row = {**record(1, QUERY), "yield_pct": percent, "outcome_status": status}
    save_sqlite_generic_index(build_generic_index([source_row]), canonical)
    build_shared_core_index(load_generic_index(canonical), shared)
    args = [QUERY, "--index", str(canonical), "--shared-index", str(shared)]
    assert main(args + ["--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["qualified_count"] == 1
    assert payload["matches"][0]["precedent"]["observation_id"] == "obs-1"
    assert main(args) == 0
    summary = capsys.readouterr().out
    assert payload["grade"] in summary and "REF1:1" in summary
    assert printed in summary
    if status == "invalid":
        assert "source yield 98" not in summary


@pytest.mark.parametrize("argument", ["--index", "--shared-index"])
def test_incomplete_custom_index_pair_exits_cleanly(argument, capsys):
    with pytest.raises(SystemExit) as error:
        main([QUERY, argument, "some.sqlite"])
    assert error.value.code == 2
    assert "must be supplied together" in capsys.readouterr().err


@pytest.mark.parametrize(("query", "excess", "warning", "message"), [
    ("Ic1n[nH]c2ccccc12>>N#Cc1n[nH]c2ccccc12", [["C", 1], ["N", 1]],
     "MISSING_REACTANT_SUSPECTED", "submitted reactants: C (1), N (1)"),
    ("CC(C)(C)OC(=O)N(CC=C)CC=C>>CC(C)(C)OC(=O)N1CC=CC1", [],
     "REACTION_COMPLETENESS_UNRESOLVED", "Product atom correspondence could not be verified"),
    ("CC", [], "INVALID_REACTION_FORMAT", "INVALID_REACTION_FORMAT"),
])
def test_ungraded_queries_explain_graph_evidence(tmp_path, capsys, query, excess, warning, message):
    canonical = tmp_path / "index.sqlite"
    shared = tmp_path / "core.sqlite"
    save_sqlite_generic_index(build_generic_index([record(1, QUERY)]), canonical)
    build_shared_core_index(load_generic_index(canonical), shared)
    args = [query, "--index", str(canonical), "--shared-index", str(shared)]
    assert main(args + ["--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["evidence_rank"] is None
    assert payload["query_context"]["product_element_excess"] == excess
    assert warning in payload["warnings"]
    assert main(args) == 0
    assert message in capsys.readouterr().out


def test_experiment_catalog_requires_an_explicit_recipe(capsys):
    with pytest.raises(SystemExit) as exc:
        main([QUERY, "--experiment-catalog", "controls.sqlite"])
    assert exc.value.code == 2
    assert "--experiment-catalog requires --recipe" in capsys.readouterr().err


def test_optional_recipe_view_preserves_grade_and_returns_source_evidence(tmp_path, capsys):
    canonical = tmp_path / "index.sqlite"
    shared = tmp_path / "core.sqlite"
    source = {**record(1, QUERY), "yield_pct": 0.0, "outcome_status": "usable"}
    save_sqlite_generic_index(build_generic_index([source]), canonical)
    build_shared_core_index(load_generic_index(canonical), shared)
    recipe = tmp_path / "proposed.json"
    recipe.write_text(json.dumps({"temperature_c": 25}), encoding="utf-8")
    args = [QUERY, "--index", str(canonical), "--shared-index", str(shared), "--match-limit", "1"]
    assert main(args + ["--json"]) == 0
    original = json.loads(capsys.readouterr().out)
    assert main(args + ["--recipe", str(recipe), "--json"]) == 0
    view = json.loads(capsys.readouterr().out)
    assert view["step_assessment"] == original
    assert view["experiment_count"] == 1
    assert view["experiments"][0]["precedent"]["observation_id"] == "obs-1"
    assert view["experiments"][0]["precedent"]["yield_pct"] == 0
    assert main(args + ["--recipe", str(recipe)]) == 0
    rendered = capsys.readouterr().out
    assert "0%" in rendered and "obs-1" in rendered


def test_published_no_base_experiment_in_text_and_json(tmp_path, capsys):
    examples = Path(__file__).resolve().parents[2] / "examples" / "common"
    source = examples / "liu_thioether_controls.jsonl"
    reaction = json.loads(source.read_bytes().splitlines()[0])["reaction"]["reaction_smiles"]
    canonical = tmp_path / "index.sqlite"
    shared = tmp_path / "core.sqlite"
    supplement = tmp_path / "experiments.sqlite"
    save_sqlite_generic_index(build_generic_index([record(1, QUERY)]), canonical)
    build_shared_core_index(load_generic_index(canonical), shared)
    build_experiment_catalog(source, supplement)
    args = [
        reaction, "--index", str(canonical), "--shared-index", str(shared),
        "--experiment-catalog", str(supplement),
        "--recipe", str(examples / "liu_thioether_no_base_recipe.json"), "--match-limit", "1",
    ]
    assert main(args + ["--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    first = result["experiments"][0]
    assert first["outcome"] == "zero_yield"
    assert first["recipe_comparison"]["status"] == "matches_stated_fields"
    assert first["index_eligibility"] == "ineligible"
    assert result["display_truncated"]
    assert main(args) == 0
    summary = capsys.readouterr().out
    for detail in ("0%", "NMR yield, internal standard", "S5 Table S1 entry 5", "10.1021/jacs.7b07390.s001"):
        assert detail in summary
