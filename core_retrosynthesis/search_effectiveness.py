"""Reproducible search ablations using the canonical multistep evaluator.

Run with ``python -m core_retrosynthesis.search_effectiveness --help``.
The development panel is reused deliberately; this is not a blind benchmark.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from dataclasses import replace
import hashlib
import html
import json
from pathlib import Path
import time
from typing import Any

from .multistep_dataset_evaluation import (
    MultistepDatasetEvaluation,
    MultistepDatasetEvaluationConfig,
    evaluate_partition_review_routes,
    write_multistep_dataset_evaluation,
)
from .multistep_panel_review import render_multistep_panel_html
from .search_exploration import load_search_exploration_policy


def evaluation_metrics(evaluation: MultistepDatasetEvaluation) -> dict[str, Any]:
    """Measure work, termination evidence, and reference recovery separately."""

    results = [case.panel_case.baseline for case in evaluation.cases]
    diagnostics = [result.diagnostics for result in results]
    return {
        **evaluation.summary,
        "supplier_complete_targets": sum(
            any(
                all(
                    leaf.terminal_evidence == "supplier_stock_portfolio"
                    for leaf in route.leaves
                )
                for route in result.routes
            )
            for result in results
        ),
        "reference_depth_exceeds_limit": sum(
            (case.panel_case.reference_maximum_depth or 0) > evaluation.config.max_depth
            for case in evaluation.cases
        ),
        **{
            field: sum(getattr(item, field) for item in diagnostics)
            for field in (
                "expanded_states",
                "one_step_calls",
                "validation_attempts",
                "beam_pruned_states",
                "dead_end_states",
                "widening_revisits",
            )
        },
        "expansion_limit_targets": sum(
            d.stopped_by_expansion_limit for d in diagnostics
        ),
    }


def paired_changes(
    baseline: MultistepDatasetEvaluation,
    variant: MultistepDatasetEvaluation,
) -> list[dict[str, Any]]:
    """Pair identical targets and expose improvements and regressions equally."""

    if (
        baseline.source_review_sha256,
        baseline.library_sha256,
        baseline.stock_index_path,
    ) != (
        variant.source_review_sha256,
        variant.library_sha256,
        variant.stock_index_path,
    ):
        raise ValueError(
            "paired evaluations require identical source, library and stock"
        )
    if [case.case_id for case in baseline.cases] != [
        case.case_id for case in variant.cases
    ]:
        raise ValueError("paired evaluations require identical ordered cases")
    rows = []
    for before, after in zip(baseline.cases, variant.cases):
        left, right = before.panel_case.baseline, after.panel_case.baseline
        rows.append(
            {
                "case_id": before.case_id,
                "selection_rank": before.selection_rank,
                "split": before.split,
                "known_actions_before": before.maximum_observed_action_matches,
                "known_actions_after": after.maximum_observed_action_matches,
                "known_action_delta": (
                    after.maximum_observed_action_matches
                    - before.maximum_observed_action_matches
                ),
                "heuristic_solved_before": bool(left.routes),
                "heuristic_solved_after": bool(right.routes),
                "validation_delta": (
                    right.diagnostics.validation_attempts
                    - left.diagnostics.validation_attempts
                ),
            }
        )
    return rows


def _evaluate_worker(
    name: str,
    config: MultistepDatasetEvaluationConfig,
    guidance: Any,
    source_review: Path,
    library_path: Path,
    stock_path: Path,
    output_dir: Path,
    hashes: dict[str, str],
) -> tuple[MultistepDatasetEvaluation, float]:
    """Evaluate one independent mode with its own index connection for lookups."""

    from cas_tools import open_stock_lookup
    from .generic_library import load_generic_library

    review = json.loads(source_review.read_text(encoding="utf-8"))
    library = load_generic_library(library_path)
    with open_stock_lookup(stock_path) as stock:
        started = time.perf_counter()
        evaluation = evaluate_partition_review_routes(
            review,
            library,
            stock,
            config=config,
            source_review_path=str(source_review.resolve()),
            source_review_sha256=hashes["review"],
            library_path=str(library_path.resolve()),
            library_sha256=hashes["library"],
            stock_index_path=str(stock_path.resolve()),
            search_guidance=guidance,
            progress=lambda i, n, target: print(
                f"[{name}] {i}/{n}: {target}", flush=True
            ),
        )
        elapsed = time.perf_counter() - started
    write_multistep_dataset_evaluation(
        evaluation,
        output_dir / f"{name}.json",
        output_dir / f"{name}.html",
    )
    return evaluation, elapsed


def run_search_effectiveness(
    source_review: Path,
    library_path: Path,
    stock_path: Path,
    output_dir: Path,
    *,
    config: MultistepDatasetEvaluationConfig,
    widening_factor: int = 3,
    catalog_path: Path | None = None,
    workers: int = 1,
) -> dict[str, Any]:
    """Run narrow, fixed-wide, deferred-wide and optional ordering ablations.

    All modes share depth/expansion/beam/per-level validation limits. Actual
    validation work is reported, not claimed equal. No chemistry or ranker is
    fitted to this panel. Each completed mode is checkpointed before continuing.
    """

    from cas_tools import open_stock_lookup
    from .generic_library import load_generic_library
    from .route_state_learning import (
        LiteratureRouteOrderingGuidance,
        load_route_state_learning_catalog,
    )

    if type(widening_factor) is not int or widening_factor < 2:
        raise ValueError("comparison widening factor must be at least two")
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError("comparison workers must be between one and four")
    if config.widening_factor != 1 or config.route_state_ordering_enabled:
        raise ValueError("comparison config must describe the ordinary baseline")
    if config.per_step_top_k * widening_factor > config.max_candidates_to_validate:
        raise ValueError("wide pool must fit within the per-level validation limit")
    output_dir.mkdir(parents=True, exist_ok=True)
    review = json.loads(source_review.read_text(encoding="utf-8"))
    library = load_generic_library(library_path)
    hashes = {}
    for name, path in (
        ("review", source_review),
        ("library", library_path),
        ("stock", stock_path),
    ):
        checksum = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                checksum.update(block)
        hashes[name] = checksum.hexdigest()
    variants = [
        ("baseline", config, None),
        (
            "fixed_wide",
            replace(config, per_step_top_k=config.per_step_top_k * widening_factor),
            None,
        ),
        ("deferred_wide", replace(config, widening_factor=widening_factor), None),
    ]
    if catalog_path is not None:
        catalog = load_route_state_learning_catalog(catalog_path)
        hashes["catalog"] = hashlib.sha256(catalog_path.read_bytes()).hexdigest()
        ordered_config = replace(
            config,
            widening_factor=widening_factor,
            route_state_definition_id=catalog.definition_id,
            route_state_catalog_sha256=hashes["catalog"],
            route_state_ordering_enabled=True,
        )
        # Isolate ordering: no state-reservation selector in this ablation.
        variants.append(
            (
                "deferred_ordered",
                ordered_config,
                LiteratureRouteOrderingGuidance(catalog),
            )
        )
    report: dict[str, Any] = {
        "schema_version": "1.0",
        "definition_id": "search_effectiveness.v1",
        "input_sha256": hashes,
        "workers": workers,
        "exploration_policy": vars(load_search_exploration_policy()),
        "warnings": [
            "Development panel; includes previously inspected targets and train patents.",
            "Solved means completion under configured terminal heuristics, not feasibility.",
            "Equal expansion limits do not imply equal validation work; inspect both.",
            (
                "Wall time is descriptive: modes run concurrently in independent processes."
                if workers > 1
                else "Wall time is descriptive: modes run sequentially with shared process caches."
            ),
            "Reference action recovery is not topology-verified route recovery.",
        ],
        "variants": {},
    }
    baseline = None
    with ExitStack() as stack:
        stock = stack.enter_context(open_stock_lookup(stock_path))
        pool = (
            stack.enter_context(ProcessPoolExecutor(max_workers=workers))
            if workers > 1
            else None
        )
        futures = (
            {
                name: pool.submit(
                    _evaluate_worker,
                    name,
                    variant_config,
                    guidance,
                    source_review,
                    library_path,
                    stock_path,
                    output_dir,
                    hashes,
                )
                for name, variant_config, guidance in variants
            }
            if pool
            else {}
        )
        for name, variant_config, guidance in variants:
            started = time.perf_counter()
            if pool:
                evaluation, elapsed = futures[name].result()
            else:
                evaluation = evaluate_partition_review_routes(
                    review,
                    library,
                    stock,
                    config=variant_config,
                    source_review_path=str(source_review.resolve()),
                    source_review_sha256=hashes["review"],
                    library_path=str(library_path.resolve()),
                    library_sha256=hashes["library"],
                    stock_index_path=str(stock_path.resolve()),
                    search_guidance=guidance,
                    progress=lambda i, n, target: print(
                        f"[{name}] {i}/{n}: {target}", flush=True
                    ),
                )
                elapsed = time.perf_counter() - started
            write_multistep_dataset_evaluation(
                evaluation,
                output_dir / f"{name}.json",
                output_dir / f"{name}.html",
            )
            if baseline is None:
                baseline = evaluation
            changes = paired_changes(baseline, evaluation)
            report["variants"][name] = {
                "config": variant_config.to_dict(),
                "metrics": evaluation_metrics(evaluation),
                "elapsed_seconds": round(elapsed, 3),
                "paired_changes": changes,
            }
            if name != "baseline":
                paired = tuple(
                    replace(
                        left.panel_case,
                        policy=right.panel_case.baseline,
                        evaluation_metrics=(
                            (
                                "Known actions: baseline / variant",
                                f"{left.maximum_observed_action_matches} / {right.maximum_observed_action_matches}",
                            ),
                            (
                                "Validations: baseline / variant",
                                f"{left.panel_case.baseline.diagnostics.validation_attempts} / "
                                f"{right.panel_case.baseline.diagnostics.validation_attempts}",
                            ),
                        ),
                    )
                    for left, right in zip(baseline.cases, evaluation.cases)
                )
                document = render_multistep_panel_html(
                    paired,
                    title=f"Search effectiveness: baseline vs {name}",
                    top_k=config.top_k_routes,
                    metadata={
                        "panel_id": f"{evaluation.evaluation_id}:{name}",
                        "warnings": report["warnings"],
                    },
                )
                (output_dir / f"{name}_paired.html").write_text(
                    document, encoding="utf-8"
                )
            (output_dir / "comparison.json").write_text(
                json.dumps(report, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    columns = (
        ("planner_solved_count", "Heuristic solved"),
        ("supplier_complete_targets", "Supplier complete"),
        ("matched_observed_action_count", "Known actions"),
        ("exact_root_action_recovery_count", "Root recovery"),
        ("validation_attempts", "Validations"),
        ("expanded_states", "Expansions"),
        ("widening_revisits", "Widening revisits"),
    )
    rows = []
    for name, result in report["variants"].items():
        link = f"{name}.html" if name == "baseline" else f"{name}_paired.html"
        rows.append(
            f'<tr><td><a href="{link}">{name}</a></td>'
            + "".join(f"<td>{result['metrics'][key]}</td>" for key, _ in columns)
            + f"<td>{result['elapsed_seconds']:.1f}</td></tr>"
        )
    document = (
        '<!doctype html><html lang="en"><meta charset="utf-8"><title>Search effectiveness</title>'
        "<style>body{font:16px system-ui;max-width:1300px;margin:40px auto;padding:20px;"
        "background:#f6f8f7;color:#233b32}table{border-collapse:collapse;background:white;"
        "width:100%}td,th{padding:12px;border:1px solid #ccd8d1;text-align:left}"
        "a{color:#165e8b}li{margin:10px 0}</style><h1>Search effectiveness</h1>"
        "<p>Open each variant for paired molecular routes, review notes, and review JSON export.</p>"
        "<ul>"
        + "".join(f"<li>{html.escape(w)}</li>" for w in report["warnings"])
        + "</ul>"
        f"<p>Depth {config.max_depth}; expansion limit {config.max_expansions}; beam {config.beam_width}; "
        f"initial top-k {config.per_step_top_k}; wide pool {config.per_step_top_k * widening_factor}.</p>"
        "<table><tr><th>Mode</th>"
        + "".join(f"<th>{label}</th>" for _, label in columns)
        + "<th>Seconds</th></tr>"
        + "".join(rows)
        + "</table>"
        '<p><a href="comparison.json">Machine-readable metrics, paired deltas and input hashes</a></p></html>'
    )
    (output_dir / "comparison.html").write_text(document, encoding="utf-8")
    return report


def main() -> None:
    """Run the fixed development ablation from local artifacts."""

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source_review", "library", "stock", "output_dir"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--catalog", type=Path)
    parser.add_argument("--max-depth", type=int, default=6)
    parser.add_argument("--max-expansions", type=int, default=30)
    parser.add_argument("--widening-factor", type=int, default=3)
    parser.add_argument("--workers", type=int, choices=(1, 2, 3, 4), default=1)
    parser.add_argument("--allow-untyped-literature-terminals", action="store_true")
    args = parser.parse_args()
    report = run_search_effectiveness(
        args.source_review,
        args.library,
        args.stock,
        args.output_dir,
        config=MultistepDatasetEvaluationConfig(
            max_depth=args.max_depth,
            max_expansions=args.max_expansions,
            allow_untyped_literature_terminals=args.allow_untyped_literature_terminals,
        ),
        widening_factor=args.widening_factor,
        catalog_path=args.catalog,
        workers=args.workers,
    )
    print(
        json.dumps(
            {key: value["metrics"] for key, value in report["variants"].items()},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
