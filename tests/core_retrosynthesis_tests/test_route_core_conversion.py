"""Route-core corpus conversion regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from core_retrosynthesis.route_conversion import build_observed_route_tree
from core_retrosynthesis.route_core_conversion import (
    convert_route_core_corpus,
    iter_route_core_projections,
)

from .test_route_core import _two_step_record


def _write_tree(path: Path) -> None:
    tree = build_observed_route_tree(_two_step_record())
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(tree.to_dict()) + "\n")


def test_route_core_corpus_is_deterministic_and_reloadable(
    tmp_path: Path,
) -> None:
    source = tmp_path / "trees.jsonl.gz"
    first = tmp_path / "first.jsonl.gz"
    second = tmp_path / "second.jsonl.gz"
    _write_tree(source)

    first_report = convert_route_core_corpus(source, first)
    second_report = convert_route_core_corpus(source, second)

    assert first.read_bytes() == second.read_bytes()
    assert first_report["conversion"]["route_count"] == 1
    assert first_report["conversion"]["reaction_count"] == 2
    assert first_report["conversion"]["reaction_signature_count"] == 2
    assert first_report["conversion"]["reaction_core_count"] == 2
    assert first_report["conversion"]["rejected_count"] == 0
    assert first_report["output"]["sha256"] == second_report["output"]["sha256"]
    projections = tuple(iter_route_core_projections(first))
    assert len(projections) == 1
    assert projections[0].fully_lineage_connected
