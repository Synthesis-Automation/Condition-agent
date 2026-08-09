"""Command-line interface for the standalone C-X retrosynthesis POC."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from typing import Optional, Sequence

from .library import build_library, iter_rows, load_library, save_library
from .search import disconnect_target


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Precedent-backed one-step C-X retrosynthesis POC"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build-library", help="extract a template library")
    build.add_argument("source", help="JSONL/JSONL.GZ file or shard directory")
    build.add_argument("output", help="output JSON or JSON.GZ library")
    build.add_argument(
        "--include",
        action="append",
        default=[],
        help="directory glob to include; repeat for multiple source cohorts",
    )
    build.add_argument("--max-rows", type=int)
    build.add_argument("--max-precedents-per-template", type=int, default=8)

    disconnect = subparsers.add_parser(
        "disconnect", help="generate precursor candidates for one target"
    )
    disconnect.add_argument("library", help="template library JSON or JSON.GZ")
    disconnect.add_argument("target", help="target molecule SMILES")
    disconnect.add_argument("--bond", action="append", choices=("C-N", "C-O", "C-S"))
    disconnect.add_argument("--max-templates", type=int, default=250)
    disconnect.add_argument("--top-k", type=int, default=20)
    disconnect.add_argument("--skip-forward-validation", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the requested build or target-disconnection operation."""

    arguments = _parser().parse_args(argv)
    if arguments.command == "build-library":
        library, report = build_library(
            iter_rows(arguments.source, include=arguments.include),
            max_rows=arguments.max_rows,
            max_precedents_per_template=arguments.max_precedents_per_template,
        )
        save_library(library, arguments.output)
        print(json.dumps(asdict(report), indent=2, sort_keys=True))
        return 0
    library = load_library(arguments.library)
    candidates = disconnect_target(
        arguments.target,
        library,
        allowed_bonds=arguments.bond or ("C-N", "C-O", "C-S"),
        max_templates_to_apply=arguments.max_templates,
        top_k=arguments.top_k,
        validate_forward=not arguments.skip_forward_validation,
    )
    print(
        json.dumps(
            [candidate.to_dict() for candidate in candidates],
            indent=2,
            sort_keys=True,
        )
    )
    return 0


__all__ = ["main"]
