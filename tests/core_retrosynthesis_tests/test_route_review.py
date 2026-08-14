"""Visual route-review report regressions."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from core_retrosynthesis import (
    render_route_review_html,
    sample_route_records,
    write_route_review_html,
)


def _record(route_id: str, patent_id: str, target: str = "CCO") -> dict:
    return {
        "schema_version": "1.0",
        "route_id": route_id,
        "patent_id": patent_id,
        "split": "train",
        "target_smiles": target,
        "abstraction_reduction": 1,
        "steps": [
            {
                "retrosynthetic_position": 0,
                "source_reaction_id": f"{route_id}-step",
                "reaction_smiles": f"CC.O>O>{target}",
                "abstracted_reaction_smiles": f"CC.O>>{target}",
                "reagents_smiles": "O",
                "internal_precursor_smiles": [],
                "terminal_precursor_smiles": ["CC", "O"],
            }
        ],
    }


def test_random_sample_is_reproducible_and_file_order_invariant() -> None:
    records = [
        _record("route-c", "patent-c"),
        _record("route-a", "patent-a"),
        _record("route-b", "patent-b"),
    ]

    first = sample_route_records(records, sample_size=2, seed=17)
    second = sample_route_records(reversed(records), sample_size=2, seed=17)

    assert [row["route_id"] for row in first] == [
        row["route_id"] for row in second
    ]


def test_route_review_contains_reaction_drawings_and_review_controls() -> None:
    document = render_route_review_html(
        [_record("route-a", "patent-a"), _record("route-b", "patent-b")],
        sample_size=2,
        seed=7,
    )

    assert document.count("class=\"route-card\"") == 2
    assert document.count("<svg") >= 6
    assert "Original recorded reaction" in document
    assert "Higher-level abstraction" in document
    assert "Synthetic sequence" in document
    assert "Starting material(s)" in document
    assert "Final product" in document
    assert "Show retrosynthesis" in document
    assert "class=\"sequence-track\"" in document
    assert "Export review JSON" in document
    assert "localStorage" in document


def test_route_review_writer_loads_gzip_jsonl(tmp_path: Path) -> None:
    source = tmp_path / "routes.jsonl.gz"
    output = tmp_path / "review.html"
    records = [_record("route-a", "patent-a"), _record("route-b", "patent-b")]
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")

    report = write_route_review_html(
        source,
        output,
        sample_size=1,
        seed=9,
    )

    assert output.is_file()
    assert report["sample_size"] == 1
    assert report["html_bytes"] == output.stat().st_size
    assert "<!doctype html>" in output.read_text(encoding="utf-8")
