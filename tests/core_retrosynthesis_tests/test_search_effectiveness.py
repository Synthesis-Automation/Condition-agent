"""End-to-end report checks on deterministic planner fixtures."""

from contextlib import nullcontext
from dataclasses import replace
import json

import pytest

from core_retrosynthesis import search_effectiveness as effectiveness
from core_retrosynthesis.multistep_dataset_evaluation import (
    MultistepDatasetEvaluationConfig,
    evaluate_partition_review_routes,
)
from core_retrosynthesis.multistep import plan_multistep_routes
from .test_multistep import _LiteratureIndex, _candidate, _expander
from .test_operator_coverage_comparison import _tree


def test_effectiveness_report_pairs_routes_and_exposes_compute(
    tmp_path, monkeypatch
) -> None:
    import cas_tools
    import core_retrosynthesis.generic_library as libraries

    source, library, stock = (
        tmp_path / name for name in ("review.json", "library", "stock")
    )
    source.write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "case_id": "fixture",
                        "selection_rank": 1,
                        "route_tree": _tree().to_dict(),
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    library.write_bytes(b"library fixture")
    stock.write_bytes(b"stock fixture")
    monkeypatch.setattr(
        libraries,
        "load_generic_library",
        lambda path: type("Library", (), {"definition": {}})(),
    )
    monkeypatch.setattr(
        cas_tools, "open_stock_lookup", lambda path: nullcontext(_LiteratureIndex())
    )
    captured = []

    def evaluate(review, loaded_library, index, **kwargs):
        def planner(target, library, index, **controls):
            return plan_multistep_routes(
                target,
                library,
                index,
                **controls,
                expander=_expander({target: (_candidate(target, "CC(=O)O.CCN"),)}),
            )

        result = evaluate_partition_review_routes(
            review,
            loaded_library,
            index,
            **kwargs,
            planner=planner,
        )
        captured.append(result)
        return result

    monkeypatch.setattr(effectiveness, "evaluate_partition_review_routes", evaluate)
    output = tmp_path / "output"
    report = effectiveness.run_search_effectiveness(
        source,
        library,
        stock,
        output,
        config=MultistepDatasetEvaluationConfig(max_depth=1, per_step_top_k=1),
    )
    assert set(report["variants"]) == {"baseline", "fixed_wide", "deferred_wide"}
    for variant in report["variants"].values():
        assert variant["metrics"]["planner_solved_count"] == 1
        assert variant["metrics"]["supplier_complete_targets"] == 0
        assert variant["paired_changes"][0]["known_action_delta"] == 0
    assert captured[0].config.per_step_top_k == 1
    assert captured[1].config.per_step_top_k == 3
    assert captured[2].config.widening_factor == 3
    document = (output / "deferred_wide_paired.html").read_text(encoding="utf-8")
    assert "<svg" in document
    assert "Export review JSON" in document
    assert "Known actions: baseline / variant" in document
    assert "Validations: baseline / variant" in document
    assert (output / "comparison.html").exists()
    assert json.loads((output / "comparison.json").read_text()) == report
    with pytest.raises(ValueError, match="identical source"):
        effectiveness.paired_changes(
            captured[0], replace(captured[1], library_sha256="different")
        )
    with pytest.raises(ValueError, match="identical ordered cases"):
        effectiveness.paired_changes(captured[0], replace(captured[1], cases=()))


def test_process_workers_preserve_serial_metrics(tmp_path) -> None:
    from cas_tools import build_canonical_molecule_index
    from core_retrosynthesis.generic_library import save_generic_library
    from core_retrosynthesis.generic_models import GenericTemplateLibrary

    source = tmp_path / "review.json"
    source.write_text(
        json.dumps(
            {
                "cases": [
                    {
                        "case_id": "fixture",
                        "selection_rank": 1,
                        "route_tree": _tree().to_dict(),
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    library = tmp_path / "library.json.gz"
    save_generic_library(
        GenericTemplateLibrary(
            templates=(),
            source_row_count=0,
            accepted_observation_count=0,
            rejection_counts={},
            definition={},
        ),
        library,
    )
    csv_path = tmp_path / "stock.csv"
    csv_path.write_text("compound_smiles\nCC\n", encoding="utf-8")
    stock = tmp_path / "stock.sqlite"
    build_canonical_molecule_index(csv_path, stock)
    config = MultistepDatasetEvaluationConfig(max_depth=1, max_expansions=1)
    serial = effectiveness.run_search_effectiveness(
        source,
        library,
        stock,
        tmp_path / "serial",
        config=config,
    )
    concurrent = effectiveness.run_search_effectiveness(
        source,
        library,
        stock,
        tmp_path / "concurrent",
        config=config,
        workers=2,
    )
    assert concurrent["workers"] == 2
    for name in serial["variants"]:
        assert (
            serial["variants"][name]["metrics"]
            == concurrent["variants"][name]["metrics"]
        )
        assert (
            serial["variants"][name]["paired_changes"]
            == concurrent["variants"][name]["paired_changes"]
        )
