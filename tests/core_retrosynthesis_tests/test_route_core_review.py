"""Minimized route-core HTML review regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from core_retrosynthesis.route_conversion import build_observed_route_tree
from core_retrosynthesis.route_core import build_route_core_projection
from core_retrosynthesis.route_core_review import (
    render_route_core_review_html,
    write_route_core_review_html,
)

from .test_route_core import _two_step_record


def _projection():
    return build_route_core_projection(
        build_observed_route_tree(_two_step_record())
    )


def test_route_core_review_contains_minimized_chemistry_and_lineage() -> None:
    document = render_route_core_review_html(
        [_projection()],
        sample_size=1,
        seed=4,
    )

    assert "Minimized route core" in document
    assert "Core event 1" in document
    assert "Target-forming core event" in document
    assert "unique atom lineage" in document
    assert "Show retrosynthesis" in document
    assert "not executable multistep templates" in document
    assert document.count("<svg") >= 3


def test_route_core_review_writer_loads_gzip_jsonl(tmp_path: Path) -> None:
    source = tmp_path / "route-cores.jsonl.gz"
    output = tmp_path / "review.html"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps(_projection().to_dict()) + "\n")

    report = write_route_core_review_html(
        source,
        output,
        sample_size=1,
        seed=8,
    )

    assert output.is_file()
    assert report["sample_size"] == 1
    assert report["html_bytes"] == output.stat().st_size
