"""Exercise standalone CLI validation and both output forms on bound SQLite files."""

import json

import pytest

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
