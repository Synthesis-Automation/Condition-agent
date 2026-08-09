"""CLI for core-derived retrosynthesis and paired evaluation."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from itertools import islice
from typing import Any, Dict, Iterable, Iterator, Optional, Sequence

from retrosynthesis_poc.library import iter_rows
from retrosynthesis_poc.library import load_library as load_baseline_library

from .comparison import run_comparison
from .diverse_benchmark import load_diverse_rows, run_diverse_benchmark
from .ensemble import disconnect_ensemble
from .generic_library import load_generic_library
from .generic_search import disconnect_generic_target
from .html_report import DEFAULT_METHODS, write_comparison_html
from .library import build_library, load_library, save_library
from .search import disconnect_target


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reaction-core-derived one-step retrosynthesis POC"
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

    generic = commands.add_parser(
        "disconnect-generic",
        help="apply a structurally diverse generic template library",
    )
    generic.add_argument("library")
    generic.add_argument("target")
    generic.add_argument("--transformation", action="append")
    generic.add_argument("--level", action="append", choices=("RDCHIRAL", "L1", "L2"))
    generic.add_argument("--top-k", type=int, default=20)
    generic.add_argument("--max-candidates-to-validate", type=int, default=50)
    generic.add_argument("--concise", action="store_true")
    generic.add_argument("--no-context", action="store_true")

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
    source: str,
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
                    iter_rows(source, include=(pattern,)),
                    max_rows_per_include,
                )

        values: Iterable[Dict[str, Any]] = balanced()
    else:
        values = iter_rows(source, include=includes)
    return islice(values, max_rows) if max_rows is not None else values


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run library construction, disconnection, or paired evaluation."""

    arguments = _parser().parse_args(argv)
    if arguments.command == "build-library":
        rows = _selected_rows(
            arguments.source,
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
            arguments.source,
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
            arguments.source,
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

    if arguments.command == "disconnect-generic":
        candidates = disconnect_generic_target(
            arguments.target,
            load_generic_library(arguments.library),
            transformations=arguments.transformation or (),
            levels=arguments.level or (),
            top_k=arguments.top_k,
            max_candidates_to_validate=arguments.max_candidates_to_validate,
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
