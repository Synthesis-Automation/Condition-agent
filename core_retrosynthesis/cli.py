"""CLI for core-derived retrosynthesis and paired evaluation."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from itertools import islice
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Optional, Sequence

from cas_tools import open_stock_lookup
from .baselines.cx_rdchiral import load_library as load_baseline_library

from .comparison import run_comparison
from .condition_ranking import rank_retrosynthesis_candidates_with_conditions
from .diverse_benchmark import (
    load_diverse_rows,
    load_stress_rows,
    run_diverse_benchmark,
)
from .ensemble import disconnect_ensemble
from .coverage_audit import audit_operator_library_coverage
from .full_scale import FullScaleBuildConfig, build_full_scale_operator_library
from .generic_library import load_generic_library
from .generic_search import (
    disconnect_generic_target,
    disconnect_generic_target_detailed,
    disconnect_operator_ladder,
)
from .html_report import DEFAULT_METHODS, write_comparison_html
from .library import build_library, load_library, save_library
from .multistep import plan_multistep_routes
from .operator_benchmark import (
    load_operator_rows,
    run_operator_coverage_benchmark,
)
from .precedent_route_expansion import run_precedent_route_expansion_poc
from .route_curation import RouteSubsetPolicy, curate_route_subset
from .route_conversion import (
    DEFAULT_OBSERVED_ROUTE_DATASET_ID,
    DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED,
    convert_observed_route_corpus,
)
from .route_core_conversion import (
    DEFAULT_ROUTE_CORE_SAMPLE_SEED,
    convert_route_core_corpus,
)
from .route_action_conversion import (
    DEFAULT_ROUTE_ACTION_SAMPLE_SEED,
    convert_route_action_corpus,
    iter_route_action_evaluations,
    merge_route_action_shards,
)
from .route_action_evaluation import RouteActionEvaluationConfig
from .route_action_policy import (
    load_route_action_policy,
    train_route_action_policy_from_replay,
)
from .route_policy_calibration import (
    RoutePolicyCalibrationConfig,
    build_route_policy_calibration_targets,
    calibrate_route_action_policy,
    save_route_policy_calibration,
)
from .multistep_panel_review import load_multistep_panel_targets
from .route_review import DEFAULT_ROUTE_REVIEW_SEED, write_route_review_html
from .route_core_review import (
    DEFAULT_ROUTE_CORE_REVIEW_SEED,
    write_route_core_review_html,
)
from .route_action_review import (
    DEFAULT_ROUTE_ACTION_REVIEW_SEED,
    write_route_action_review_html,
)
from .search import disconnect_target
from .sources import LIBRARY_MODES, iter_library_rows, resolve_library_mode


def _format_elapsed(seconds: object) -> str:
    value = max(0, int(float(seconds or 0)))
    hours, remainder = divmod(value, 3600)
    minutes, remaining_seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{remaining_seconds:02d}"


def _print_build_progress(event: Dict[str, Any]) -> None:
    """Write compact full-build progress without contaminating JSON stdout."""

    phase = str(event.get("phase") or "unknown")
    elapsed = _format_elapsed(event.get("elapsed_seconds"))
    completed = int(event.get("completed_shards") or 0)
    total = int(event.get("total_shards") or event.get("source_shards") or 0)
    percent = 100.0 * completed / total if total else 0.0
    rows = int(event.get("source_rows") or 0)
    if phase == "compile":
        accepted = int(event.get("accepted_observations") or 0)
        reused = int(event.get("reused_shards") or 0)
        active = int(event.get("active_shards") or 0)
        queued = int(event.get("queued_shards") or 0)
        newly_completed = int(event.get("newly_completed_shards") or 0)
        newly_reused = int(event.get("newly_reused_shards") or 0)
        prefix = (
            f"[compile] {completed}/{total} shards ({percent:.1f}%), "
            f"{rows:,} completed rows, {accepted:,} accepted, {reused} reused"
        )
        if newly_completed > newly_reused or completed == total:
            rate = rows / max(1.0, float(event.get("elapsed_seconds") or 0))
            message = (
                f"{prefix}, completed-output average {rate:.1f} rows/s, "
                f"{active} active, {queued} queued, elapsed {elapsed}"
            )
        else:
            message = (
                f"{prefix}, {active} active, {queued} queued, elapsed {elapsed}"
            )
    elif phase == "merge":
        templates = int(event.get("template_count") or 0)
        message = (
            f"[merge] {completed}/{total} shards ({percent:.1f}%), "
            f"{rows:,} rows, {templates:,} templates, elapsed {elapsed}"
        )
    elif phase == "finalize":
        templates = int(event.get("template_count") or 0)
        message = (
            f"[finalize] calculating support and index for "
            f"{templates:,} templates, elapsed {elapsed}"
        )
    else:
        accepted = int(event.get("accepted_observations") or 0)
        templates = int(event.get("template_count") or 0)
        message = (
            f"[complete] {rows:,} rows, {accepted:,} accepted, "
            f"{templates:,} templates, elapsed {elapsed}"
        )
    print(message, file=sys.stderr, flush=True)


def _add_library_mode_argument(command: argparse.ArgumentParser) -> None:
    command.add_argument(
        "--library-mode",
        choices=LIBRARY_MODES,
        default="full",
        help=(
            "recommendation-library mode when SOURCE contains full/ and "
            "compact/ directories (default: full)"
        ),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Chemistry-first, type-agnostic retrosynthesis"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build-library", help="compile a core library")
    build.add_argument("source")
    build.add_argument("output")
    build.add_argument("--include", action="append", default=[])
    build.add_argument("--level", action="append", choices=("L1", "L2"))
    build.add_argument("--max-rows", type=int)
    build.add_argument("--max-rows-per-include", type=int)
    build.add_argument("--max-precedents-per-template", type=int, default=8)
    _add_library_mode_argument(build)

    disconnect = commands.add_parser(
        "disconnect",
        help="generate core-derived precursors for one target",
    )
    disconnect.add_argument("library")
    disconnect.add_argument("target")
    disconnect.add_argument("--bond", action="append", choices=("C-N", "C-O", "C-S"))
    disconnect.add_argument("--level", action="append", choices=("L1", "L2"))
    disconnect.add_argument("--max-templates", type=int, default=250)
    disconnect.add_argument("--max-candidates-to-validate", type=int, default=50)
    disconnect.add_argument("--top-k", type=int, default=20)
    disconnect.add_argument("--concise", action="store_true")
    disconnect.add_argument("--no-context", action="store_true")
    disconnect.add_argument("--skip-forward-validation", action="store_true")

    ensemble = commands.add_parser(
        "disconnect-ensemble",
        help="use the RDChiral baseline with contextual L1 fallback",
    )
    ensemble.add_argument("baseline_library")
    ensemble.add_argument("core_library")
    ensemble.add_argument("target")
    ensemble.add_argument("--bond", action="append", choices=("C-N", "C-O", "C-S"))
    ensemble.add_argument("--top-k", type=int, default=20)
    ensemble.add_argument("--max-candidates-to-validate", type=int, default=50)
    ensemble.add_argument("--concise", action="store_true")
    ensemble.add_argument("--skip-forward-validation", action="store_true")

    compare = commands.add_parser(
        "compare",
        help="build both methods and evaluate a held-out reference split",
    )
    compare.add_argument("source")
    compare.add_argument("output_directory")
    compare.add_argument("--include", action="append", default=[])
    compare.add_argument("--max-rows", type=int)
    compare.add_argument("--max-rows-per-include", type=int)
    compare.add_argument("--test-fraction", type=float, default=0.2)
    compare.add_argument("--top-k", type=int, default=10)
    compare.add_argument("--max-test-targets", type=int)
    compare.add_argument("--max-candidates-to-validate", type=int, default=20)
    _add_library_mode_argument(compare)

    diverse = commands.add_parser(
        "compare-diverse",
        help="evaluate generic operators across diverse edit archetypes",
    )
    diverse.add_argument("source")
    diverse.add_argument("output_directory")
    diverse.add_argument("--max-rows-per-cohort", type=int, default=200)
    diverse.add_argument("--test-fraction", type=float, default=0.2)
    diverse.add_argument("--max-targets-per-transformation", type=int, default=10)
    diverse.add_argument("--top-k", type=int, default=10)
    diverse.add_argument("--max-candidates-to-validate", type=int, default=30)
    _add_library_mode_argument(diverse)

    stress = commands.add_parser(
        "compare-stress",
        help="stress-test C-C coupling and conjugate-addition site selection",
    )
    stress.add_argument("source")
    stress.add_argument("output_directory")
    stress.add_argument("--max-rows-per-cohort", type=int, default=1_000)
    stress.add_argument("--test-fraction", type=float, default=0.25)
    stress.add_argument("--max-targets-per-transformation", type=int, default=50)
    stress.add_argument("--top-k", type=int, default=10)
    stress.add_argument("--max-candidates-to-validate", type=int, default=75)
    _add_library_mode_argument(stress)

    operators = commands.add_parser(
        "compare-operators",
        help="mine unrestricted graph operators and benchmark coverage",
    )
    operators.add_argument("source")
    operators.add_argument("output_directory")
    operators.add_argument("--max-rows", type=int, default=1_000)
    operators.add_argument("--test-fraction", type=float, default=0.25)
    operators.add_argument("--max-targets", type=int, default=80)
    operators.add_argument("--max-targets-per-operator", type=int, default=3)
    operators.add_argument("--top-k", type=int, default=25)
    operators.add_argument("--max-templates", type=int, default=500)
    operators.add_argument("--max-candidates-to-validate", type=int, default=100)
    _add_library_mode_argument(operators)

    full_build = commands.add_parser(
        "build-operators-full",
        help="compile and merge a resumable data-derived operator library",
    )
    full_build.add_argument("source")
    full_build.add_argument("output_directory")
    full_build.add_argument("--max-shards", type=int)
    full_build.add_argument("--max-rows-per-shard", type=int)
    full_build.add_argument("--max-precedents-per-template", type=int, default=8)
    full_build.add_argument("--workers", type=int, default=1)
    full_build.add_argument("--progress-interval", type=float, default=30.0)
    full_build.add_argument("--quiet-progress", action="store_true")
    full_build.add_argument("--skip-l0", action="store_true")
    full_build.add_argument("--force", action="store_true")
    _add_library_mode_argument(full_build)

    coverage_audit = commands.add_parser(
        "audit-operator-coverage",
        help="attribute held-out coverage misses to exact pipeline stages",
    )
    coverage_audit.add_argument("library")
    coverage_audit.add_argument("source")
    coverage_audit.add_argument("output_directory")
    coverage_audit.add_argument("--max-rows", type=int, default=1_000)
    coverage_audit.add_argument("--top-k", type=int, default=25)
    coverage_audit.add_argument("--max-templates", type=int, default=500)
    coverage_audit.add_argument(
        "--max-candidates-to-validate",
        type=int,
        default=100,
    )
    _add_library_mode_argument(coverage_audit)

    generic = commands.add_parser(
        "disconnect-generic",
        help="apply a structurally diverse generic template library",
    )
    generic.add_argument("library")
    generic.add_argument("target")
    generic.add_argument("--transformation", action="append")
    generic.add_argument(
        "--level",
        action="append",
        choices=("RDCHIRAL", "L0", "L1", "L2"),
    )
    generic.add_argument("--top-k", type=int, default=20)
    generic.add_argument("--max-candidates-to-validate", type=int, default=50)
    generic.add_argument("--concise", action="store_true")
    generic.add_argument("--no-context", action="store_true")
    generic.add_argument("--diversify-sites", action="store_true")

    operator_search = commands.add_parser(
        "disconnect-operators",
        help="apply data-derived operators with L2-to-L1-to-L0 fallback",
    )
    operator_search.add_argument("library")
    operator_search.add_argument("target")
    operator_search.add_argument("--top-k", type=int, default=20)
    operator_search.add_argument("--max-templates", type=int, default=500)
    operator_search.add_argument(
        "--max-candidates-to-validate",
        type=int,
        default=100,
    )
    operator_search.add_argument("--concise", action="store_true")
    operator_search.add_argument("--no-context", action="store_true")
    operator_search.add_argument("--skip-l0", action="store_true")
    operator_search.add_argument(
        "--no-diversity",
        action="store_true",
        help=(
            "disable hierarchical SITE1/SYN1/REAL1 ranking and "
            "operator/site/synthon diversity for an ablation"
        ),
    )
    operator_search.add_argument(
        "--no-hierarchical-ranking",
        action="store_true",
        help="retain legacy diversity ranking without SITE1/SYN1/REAL1 evidence",
    )
    operator_search.add_argument(
        "--diagnostics",
        action="store_true",
        help="include retrieval and validation stage counters",
    )
    operator_search.add_argument(
        "--condition-index",
        help=(
            "generic condition index used to enrich and rerank verified "
            "retrosynthesis candidates"
        ),
    )
    operator_search.add_argument("--condition-top-k", type=int, default=3)
    operator_search.add_argument("--condition-minimum-pool-size", type=int)
    operator_search.add_argument(
        "--condition-unrestricted-fallback",
        action="store_true",
        help="allow review-core condition precedents when the index permits it",
    )
    operator_search.add_argument(
        "--condition-use-rxnmapper",
        action="store_true",
        help="use the configured external mapper for condition queries",
    )
    operator_search.add_argument(
        "--keep-retrosynthesis-order",
        action="store_true",
        help="attach condition evidence without condition-informed reranking",
    )

    route_search = commands.add_parser(
        "plan-routes",
        help="search deterministic two- or three-depth retrosynthesis routes",
    )
    route_search.add_argument("library")
    route_search.add_argument("stock_index")
    route_search.add_argument("target")
    route_search.add_argument("--max-depth", type=int, choices=(2, 3), default=3)
    route_search.add_argument("--molecular-weight-threshold", type=float, default=150.0)
    route_search.add_argument("--top-k-routes", type=int, default=5)
    route_search.add_argument("--per-step-top-k", type=int, default=5)
    route_search.add_argument("--beam-width", type=int, default=20)
    route_search.add_argument("--max-expansions", type=int, default=100)
    route_search.add_argument("--max-templates", type=int, default=300)
    route_search.add_argument("--max-candidates-to-validate", type=int, default=50)
    route_search.add_argument("--no-context", action="store_true")
    route_search.add_argument("--skip-l0", action="store_true")
    route_search.add_argument("--no-diversity", action="store_true")
    route_search.add_argument("--no-hierarchical-ranking", action="store_true")
    route_search.add_argument(
        "--route-action-policy",
        help="optional learned policy over already validated one-step actions",
    )
    route_search.add_argument(
        "--allow-untyped-literature-terminals",
        action="store_true",
        help=(
            "temporarily accept legacy catalog matches without source-role "
            "provenance"
        ),
    )

    route_curation = commands.add_parser(
        "curate-route-corpus",
        help="validate and extract a patent-disjoint multistep route corpus",
    )
    route_curation.add_argument("source_routes")
    route_curation.add_argument("output_jsonl")
    route_curation.add_argument("--manifest")
    route_curation.add_argument("--testset")
    route_curation.add_argument("--minimum-steps", type=int, default=3)
    route_curation.add_argument("--maximum-steps", type=int, default=6)
    route_curation.add_argument("--maximum-routes", type=int, default=5_000)
    route_curation.add_argument(
        "--allow-no-abstraction-reduction",
        action="store_true",
    )
    route_curation.add_argument("--overwrite", action="store_true")

    route_review = commands.add_parser(
        "render-route-review",
        help="render a random curated-route sample as self-contained HTML",
    )
    route_review.add_argument("source_routes")
    route_review.add_argument("output_html")
    route_review.add_argument("--sample-size", type=int, default=50)
    route_review.add_argument("--seed", type=int, default=DEFAULT_ROUTE_REVIEW_SEED)
    route_review.add_argument("--title", default="Random route review")
    route_review.add_argument(
        "--sequence-only",
        action="store_true",
        help="omit detailed retrosynthetic step panels",
    )

    route_conversion = commands.add_parser(
        "convert-route-trees",
        help="convert curated observed routes to the canonical tree schema",
    )
    route_conversion.add_argument("source_routes")
    route_conversion.add_argument("output_jsonl")
    route_conversion.add_argument("--manifest")
    route_conversion.add_argument(
        "--source-dataset-id",
        default=DEFAULT_OBSERVED_ROUTE_DATASET_ID,
    )
    route_conversion.add_argument("--sample-size", type=int)
    route_conversion.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_OBSERVED_ROUTE_SAMPLE_SEED,
    )
    route_conversion.add_argument("--allow-rejections", action="store_true")
    route_conversion.add_argument("--overwrite", action="store_true")

    route_core = commands.add_parser(
        "build-route-cores",
        help="build minimized chemistry projections from canonical route trees",
    )
    route_core.add_argument("source_trees")
    route_core.add_argument("output_jsonl")
    route_core.add_argument("--manifest")
    route_core.add_argument("--sample-size", type=int)
    route_core.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_ROUTE_CORE_SAMPLE_SEED,
    )
    route_core.add_argument("--allow-rejections", action="store_true")
    route_core.add_argument("--overwrite", action="store_true")
    route_core.add_argument("--workers", type=int, default=1)

    route_core_review = commands.add_parser(
        "render-route-core-review",
        help="render minimized route-core projections as self-contained HTML",
    )
    route_core_review.add_argument("source_route_cores")
    route_core_review.add_argument("output_html")
    route_core_review.add_argument("--sample-size", type=int, default=50)
    route_core_review.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_ROUTE_CORE_REVIEW_SEED,
    )
    route_core_review.add_argument(
        "--title",
        default="Minimized route-core review",
    )

    route_action = commands.add_parser(
        "evaluate-route-actions",
        help="replay observed route steps through validated single-step search",
    )
    route_action.add_argument("source_trees")
    route_action.add_argument("operator_library")
    route_action.add_argument("output_jsonl")
    route_action.add_argument("--manifest")
    route_action.add_argument("--sample-size", type=int)
    route_action.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_ROUTE_ACTION_SAMPLE_SEED,
    )
    route_action.add_argument("--top-k", type=int, default=25)
    route_action.add_argument("--max-templates-to-apply", type=int, default=500)
    route_action.add_argument(
        "--max-candidates-to-validate", type=int, default=100
    )
    route_action.add_argument("--minimum-candidates-per-level", type=int, default=0)
    route_action.add_argument("--no-context", action="store_true")
    route_action.add_argument("--skip-l0", action="store_true")
    route_action.add_argument("--no-diversity", action="store_true")
    route_action.add_argument("--no-hierarchical-ranking", action="store_true")
    route_action.add_argument("--lazy-validation", action="store_true")
    route_action.add_argument(
        "--labels-only",
        action="store_true",
        help="extract supervision facets without running candidate search",
    )
    route_action.add_argument("--workers", type=int, default=1)
    route_action.add_argument("--shard-count", type=int, default=1)
    route_action.add_argument("--shard-index", type=int, default=0)
    route_action.add_argument(
        "--resume",
        action="store_true",
        help="reuse an existing shard only after validating its manifest",
    )
    route_action.add_argument("--allow-rejections", action="store_true")
    route_action.add_argument("--overwrite", action="store_true")

    route_action_merge = commands.add_parser(
        "merge-route-action-shards",
        help="validate and merge a complete deterministic route-action shard set",
    )
    route_action_merge.add_argument("output_jsonl")
    route_action_merge.add_argument("shard_outputs", nargs="+")
    route_action_merge.add_argument("--manifest")
    route_action_merge.add_argument("--overwrite", action="store_true")

    route_action_review = commands.add_parser(
        "render-route-action-review",
        help="render promoted observed-action labels for chemist review",
    )
    route_action_review.add_argument("source_evaluations")
    route_action_review.add_argument("output_html")
    route_action_review.add_argument("--sample-size", type=int, default=100)
    route_action_review.add_argument(
        "--seed", type=int, default=DEFAULT_ROUTE_ACTION_REVIEW_SEED
    )
    route_action_review.add_argument(
        "--title", default="Promoted route-action label review"
    )

    route_policy = commands.add_parser(
        "train-route-action-policy",
        help="train a listwise multistep action policy from replay choices",
    )
    route_policy.add_argument("source_replay")
    route_policy.add_argument("output_model")
    route_policy.add_argument("--report")
    route_policy.add_argument("--overwrite", action="store_true")

    route_policy_calibration = commands.add_parser(
        "calibrate-route-action-policy",
        help="freeze learned policy influence using whole-route validation",
    )
    route_policy_calibration.add_argument("source_replay")
    route_policy_calibration.add_argument("library")
    route_policy_calibration.add_argument("stock_index")
    route_policy_calibration.add_argument("input_model")
    route_policy_calibration.add_argument("validation_panel")
    route_policy_calibration.add_argument("output_model")
    route_policy_calibration.add_argument("--report")
    route_policy_calibration.add_argument("--max-depth", type=int, default=3)
    route_policy_calibration.add_argument(
        "--molecular-weight-threshold", type=float, default=80.0
    )
    route_policy_calibration.add_argument("--top-k-routes", type=int, default=3)
    route_policy_calibration.add_argument("--per-step-top-k", type=int, default=4)
    route_policy_calibration.add_argument("--beam-width", type=int, default=15)
    route_policy_calibration.add_argument("--max-expansions", type=int, default=15)
    route_policy_calibration.add_argument("--max-templates", type=int, default=100)
    route_policy_calibration.add_argument(
        "--max-candidates-to-validate", type=int, default=20
    )
    route_policy_calibration.add_argument(
        "--minimum-validation-targets", type=int, default=5
    )
    route_policy_calibration.add_argument("--overwrite", action="store_true")

    precedent_expansion = commands.add_parser(
        "expand-precedent-routes",
        help="enumerate nested chemical space from two-step precedent routes",
    )
    precedent_expansion.add_argument("definition")
    precedent_expansion.add_argument("output_json")

    report = commands.add_parser(
        "render-report",
        help="render comparison JSON as a self-contained chemistry HTML review",
    )
    report.add_argument("comparison_json")
    report.add_argument("output_html")
    report.add_argument("--method", action="append")
    report.add_argument("--top-k", type=int, default=5)
    report.add_argument("--max-targets", type=int)
    report.add_argument("--title", default="Retrosynthesis comparison review")
    return parser


def _selected_rows(
    source: str | Path,
    includes: Sequence[str],
    *,
    max_rows: int | None,
    max_rows_per_include: int | None,
) -> Iterable[Dict[str, Any]]:
    if max_rows_per_include is not None and max_rows_per_include < 1:
        raise ValueError("max rows per include must be positive")
    if max_rows_per_include is not None and includes:
        def balanced() -> Iterator[Dict[str, Any]]:
            for pattern in includes:
                yield from islice(
                    iter_library_rows(source, include=(pattern,)),
                    max_rows_per_include,
                )

        values: Iterable[Dict[str, Any]] = balanced()
    else:
        values = iter_library_rows(source, include=includes)
    return islice(values, max_rows) if max_rows is not None else values


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run library construction, disconnection, or paired evaluation."""

    arguments = _parser().parse_args(argv)
    if arguments.command == "curate-route-corpus":
        report = curate_route_subset(
            arguments.source_routes,
            arguments.output_jsonl,
            manifest_path=arguments.manifest,
            testset_path=arguments.testset,
            policy=RouteSubsetPolicy(
                minimum_steps=arguments.minimum_steps,
                maximum_steps=arguments.maximum_steps,
                maximum_routes=arguments.maximum_routes,
                require_abstraction_reduction=(
                    not arguments.allow_no_abstraction_reduction
                ),
            ),
            overwrite=arguments.overwrite,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "render-route-review":
        report = write_route_review_html(
            arguments.source_routes,
            arguments.output_html,
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            title=arguments.title,
            include_step_details=not arguments.sequence_only,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "convert-route-trees":
        report = convert_observed_route_corpus(
            arguments.source_routes,
            arguments.output_jsonl,
            manifest_path=arguments.manifest,
            source_dataset_id=arguments.source_dataset_id,
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            overwrite=arguments.overwrite,
            strict=not arguments.allow_rejections,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "build-route-cores":
        report = convert_route_core_corpus(
            arguments.source_trees,
            arguments.output_jsonl,
            manifest_path=arguments.manifest,
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            overwrite=arguments.overwrite,
            strict=not arguments.allow_rejections,
            workers=arguments.workers,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "render-route-core-review":
        report = write_route_core_review_html(
            arguments.source_route_cores,
            arguments.output_html,
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            title=arguments.title,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "evaluate-route-actions":
        report = convert_route_action_corpus(
            arguments.source_trees,
            arguments.operator_library,
            arguments.output_jsonl,
            manifest_path=arguments.manifest,
            config=RouteActionEvaluationConfig(
                top_k=arguments.top_k,
                max_templates_to_apply=arguments.max_templates_to_apply,
                max_candidates_to_validate=arguments.max_candidates_to_validate,
                use_context=not arguments.no_context,
                include_l0=not arguments.skip_l0,
                diversify=not arguments.no_diversity,
                use_hierarchical_ranking=(
                    not arguments.no_hierarchical_ranking
                ),
                minimum_candidates_per_level=(
                    arguments.minimum_candidates_per_level
                ),
                lazy_validation=arguments.lazy_validation,
                run_search=not arguments.labels_only,
            ),
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            overwrite=arguments.overwrite,
            strict=not arguments.allow_rejections,
            workers=arguments.workers,
            shard_count=arguments.shard_count,
            shard_index=arguments.shard_index,
            resume=arguments.resume,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "merge-route-action-shards":
        report = merge_route_action_shards(
            arguments.shard_outputs,
            arguments.output_jsonl,
            manifest_path=arguments.manifest,
            overwrite=arguments.overwrite,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "render-route-action-review":
        report = write_route_action_review_html(
            arguments.source_evaluations,
            arguments.output_html,
            sample_size=arguments.sample_size,
            seed=arguments.seed,
            title=arguments.title,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "train-route-action-policy":
        report = train_route_action_policy_from_replay(
            arguments.source_replay,
            arguments.output_model,
            report_path=arguments.report,
            overwrite=arguments.overwrite,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "calibrate-route-action-policy":
        panel_targets = load_multistep_panel_targets(arguments.validation_panel)
        calibration_targets = build_route_policy_calibration_targets(
            iter_route_action_evaluations(arguments.source_replay),
            tuple(target.target_id for target in panel_targets),
        )
        panel_smiles = {target.target_id: target.smiles for target in panel_targets}
        for target in calibration_targets:
            if panel_smiles[target.route_id] != target.target_smiles:
                raise ValueError(
                    f"calibration panel target contradicts replay: {target.route_id}"
                )
        loaded_library = load_generic_library(arguments.library)
        input_model = load_route_action_policy(arguments.input_model)
        with open_stock_lookup(arguments.stock_index) as stock_index:
            calibrated_model, calibration_report = calibrate_route_action_policy(
                input_model,
                calibration_targets,
                loaded_library,
                stock_index,
                config=RoutePolicyCalibrationConfig(
                    max_depth=arguments.max_depth,
                    molecular_weight_threshold=(
                        arguments.molecular_weight_threshold
                    ),
                    top_k_routes=arguments.top_k_routes,
                    per_step_top_k=arguments.per_step_top_k,
                    beam_width=arguments.beam_width,
                    max_expansions=arguments.max_expansions,
                    max_templates_to_apply=arguments.max_templates,
                    max_candidates_to_validate=(
                        arguments.max_candidates_to_validate
                    ),
                    minimum_validation_targets=(
                        arguments.minimum_validation_targets
                    ),
                ),
                progress=lambda message: print(message, flush=True),
            )
        report_path = arguments.report or f"{arguments.output_model}.calibration.json"
        report = save_route_policy_calibration(
            calibrated_model,
            calibration_report,
            arguments.output_model,
            report_path,
            source_replay=arguments.source_replay,
            input_model=arguments.input_model,
            overwrite=arguments.overwrite,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    if arguments.command == "expand-precedent-routes":
        report = run_precedent_route_expansion_poc(
            arguments.definition,
            output_path=arguments.output_json,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    source = (
        resolve_library_mode(arguments.source, arguments.library_mode)
        if hasattr(arguments, "source") and hasattr(arguments, "library_mode")
        else None
    )
    if arguments.command == "build-library":
        rows = _selected_rows(
            source,
            arguments.include,
            max_rows=arguments.max_rows,
            max_rows_per_include=arguments.max_rows_per_include,
        )
        library, report = build_library(
            rows,
            levels=arguments.level or ("L1", "L2"),
            max_precedents_per_template=arguments.max_precedents_per_template,
        )
        save_library(library, arguments.output)
        print(json.dumps(asdict(report), indent=2, sort_keys=True))
        return 0
    if arguments.command == "compare":
        rows = _selected_rows(
            source,
            arguments.include,
            max_rows=arguments.max_rows,
            max_rows_per_include=arguments.max_rows_per_include,
        )
        report = run_comparison(
            rows,
            arguments.output_directory,
            test_fraction=arguments.test_fraction,
            top_k=arguments.top_k,
            max_test_targets=arguments.max_test_targets,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "compare-diverse":
        rows = load_diverse_rows(
            source,
            max_rows_per_cohort=arguments.max_rows_per_cohort,
        )
        report = run_diverse_benchmark(
            rows,
            arguments.output_directory,
            test_fraction=arguments.test_fraction,
            max_targets_per_transformation=(
                arguments.max_targets_per_transformation
            ),
            top_k=arguments.top_k,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "compare-stress":
        rows = load_stress_rows(
            source,
            max_rows_per_cohort=arguments.max_rows_per_cohort,
        )
        report = run_diverse_benchmark(
            rows,
            arguments.output_directory,
            test_fraction=arguments.test_fraction,
            max_targets_per_transformation=(
                arguments.max_targets_per_transformation
            ),
            top_k=arguments.top_k,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
            include_level_ablations=False,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "compare-operators":
        rows = load_operator_rows(
            source,
            max_rows=arguments.max_rows,
        )
        report = run_operator_coverage_benchmark(
            rows,
            arguments.output_directory,
            test_fraction=arguments.test_fraction,
            max_targets=arguments.max_targets,
            max_targets_per_operator=arguments.max_targets_per_operator,
            top_k=arguments.top_k,
            max_templates_to_apply=arguments.max_templates,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
        )
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "build-operators-full":
        levels = ("L1", "L2") if arguments.skip_l0 else ("L0", "L1", "L2")
        _, build_report = build_full_scale_operator_library(
            source,
            arguments.output_directory,
            config=FullScaleBuildConfig(
                levels=levels,
                max_precedents_per_template=(
                    arguments.max_precedents_per_template
                ),
                max_rows_per_shard=arguments.max_rows_per_shard,
            ),
            max_shards=arguments.max_shards,
            workers=arguments.workers,
            force=arguments.force,
            progress_callback=(
                None if arguments.quiet_progress else _print_build_progress
            ),
            progress_interval_seconds=arguments.progress_interval,
        )
        print(json.dumps(build_report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "audit-operator-coverage":
        rows = islice(iter_library_rows(source), arguments.max_rows)
        audit_report = audit_operator_library_coverage(
            rows,
            load_generic_library(arguments.library),
            arguments.output_directory,
            top_k=arguments.top_k,
            max_templates_to_apply=arguments.max_templates,
            max_candidates_to_validate=(
                arguments.max_candidates_to_validate
            ),
        )
        print(json.dumps(audit_report, indent=2, sort_keys=True))
        return 0

    if arguments.command == "disconnect-generic":
        candidates = disconnect_generic_target(
            arguments.target,
            load_generic_library(arguments.library),
            transformations=arguments.transformation or (),
            levels=arguments.level or (),
            top_k=arguments.top_k,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
            use_context=not arguments.no_context,
            diversify_sites=arguments.diversify_sites,
        )
        if arguments.concise:
            for candidate in candidates:
                print(candidate.proposed_reaction_smiles)
            return 0
        print(
            json.dumps(
                [candidate.to_dict() for candidate in candidates],
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    if arguments.command == "disconnect-operators":
        loaded_library = load_generic_library(arguments.library)
        candidates = disconnect_operator_ladder(
            arguments.target,
            loaded_library,
            top_k=arguments.top_k,
            max_templates_to_apply=arguments.max_templates,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
            use_context=not arguments.no_context,
            include_l0=not arguments.skip_l0,
            diversify=not arguments.no_diversity,
            use_hierarchical_ranking=not arguments.no_hierarchical_ranking,
        )
        condition_ranked = None
        if arguments.condition_index:
            from condition_recommender import GenericConditionRecommender
            from reactive_taxonomy import RxnMapperProvider

            if (
                arguments.condition_use_rxnmapper
                and not RxnMapperProvider.is_available()
            ):
                raise RuntimeError("RXNMAPPER_UNAVAILABLE")
            condition_recommender = GenericConditionRecommender.from_path(
                arguments.condition_index,
                mapping_provider=(
                    RxnMapperProvider()
                    if arguments.condition_use_rxnmapper
                    else None
                ),
                include_review=arguments.condition_unrestricted_fallback,
            )
            condition_ranked = rank_retrosynthesis_candidates_with_conditions(
                candidates,
                condition_recommender,
                condition_top_k=arguments.condition_top_k,
                minimum_pool_size=arguments.condition_minimum_pool_size,
                unrestricted_fallback=(
                    arguments.condition_unrestricted_fallback
                ),
                rerank=not arguments.keep_retrosynthesis_order,
            )
        if arguments.concise:
            for candidate in (
                (value.candidate for value in condition_ranked)
                if condition_ranked is not None
                else candidates
            ):
                print(candidate.proposed_reaction_smiles)
            return 0
        diagnostics = None
        if arguments.diagnostics:
            _, diagnostics = disconnect_generic_target_detailed(
                arguments.target,
                loaded_library,
                levels=("L1", "L2"),
                top_k=arguments.top_k,
                max_templates_to_apply=arguments.max_templates,
                max_candidates_to_validate=(
                    arguments.max_candidates_to_validate
                ),
                use_context=not arguments.no_context,
            )
        serialized_candidates = (
            [value.to_dict() for value in condition_ranked]
            if condition_ranked is not None
            else [candidate.to_dict() for candidate in candidates]
        )
        print(
            json.dumps(
                {
                    "candidates": serialized_candidates,
                    "diagnostics": (
                        diagnostics.to_dict() if diagnostics is not None else None
                    ),
                }
                if arguments.diagnostics
                else serialized_candidates,
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    if arguments.command == "plan-routes":
        loaded_library = load_generic_library(arguments.library)
        with open_stock_lookup(arguments.stock_index) as stock_index:
            result = plan_multistep_routes(
                arguments.target,
                loaded_library,
                stock_index,
                max_depth=arguments.max_depth,
                molecular_weight_threshold=(
                    arguments.molecular_weight_threshold
                ),
                top_k_routes=arguments.top_k_routes,
                per_step_top_k=arguments.per_step_top_k,
                beam_width=arguments.beam_width,
                max_expansions=arguments.max_expansions,
                max_templates_to_apply=arguments.max_templates,
                max_candidates_to_validate=(
                    arguments.max_candidates_to_validate
                ),
                use_context=not arguments.no_context,
                include_l0=not arguments.skip_l0,
                diversify=not arguments.no_diversity,
                use_hierarchical_ranking=not arguments.no_hierarchical_ranking,
                allow_untyped_literature_terminals=(
                    arguments.allow_untyped_literature_terminals
                ),
                route_action_policy=(
                    load_route_action_policy(arguments.route_action_policy)
                    if arguments.route_action_policy
                    else None
                ),
            )
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
        return 0

    if arguments.command == "render-report":
        summary = write_comparison_html(
            arguments.comparison_json,
            arguments.output_html,
            methods=arguments.method or DEFAULT_METHODS,
            top_k=arguments.top_k,
            max_targets=arguments.max_targets,
            title=arguments.title,
        )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if arguments.command == "disconnect-ensemble":
        ensemble_candidates = disconnect_ensemble(
            arguments.target,
            load_baseline_library(arguments.baseline_library),
            load_library(arguments.core_library),
            allowed_bonds=arguments.bond or ("C-N", "C-O", "C-S"),
            top_k=arguments.top_k,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
            validate_forward=not arguments.skip_forward_validation,
        )
        if arguments.concise:
            for candidate in ensemble_candidates:
                print(candidate.proposed_reaction_smiles)
            return 0
        print(
            json.dumps(
                [candidate.to_dict() for candidate in ensemble_candidates],
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    library = load_library(arguments.library)
    candidates = disconnect_target(
        arguments.target,
        library,
        allowed_bonds=arguments.bond or ("C-N", "C-O", "C-S"),
        levels=arguments.level or ("L1", "L2"),
        max_templates_to_apply=arguments.max_templates,
        max_candidates_to_validate=arguments.max_candidates_to_validate,
        top_k=arguments.top_k,
        validate_forward=not arguments.skip_forward_validation,
        use_context=not arguments.no_context,
    )
    if arguments.concise:
        for candidate in candidates:
            print(candidate.proposed_reaction_smiles)
        return 0
    print(
        json.dumps(
            [candidate.to_dict() for candidate in candidates],
            indent=2,
            sort_keys=True,
        )
    )
    return 0


__all__ = ["main"]
