"""Command-line interface for forward operator libraries and prediction."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .library import (
    build_forward_library,
    load_forward_library,
    save_forward_library,
)
from .evaluation import evaluate_product_hidden_replay
from .prediction import assess_proposed_step, predict_products


def _read_json(path: str | Path) -> Any:
    source = Path(path)
    if source.suffix.casefold() == ".gz":
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    return json.loads(source.read_text(encoding="utf-8"))


def _recipe(path: str | None) -> Mapping[str, Any] | None:
    if not path:
        return None
    value = _read_json(path)
    if not isinstance(value, Mapping):
        raise ValueError("recipe JSON must contain one object")
    return value


def _rows(path: str | Path) -> tuple[Mapping[str, Any], ...]:
    source = Path(path)
    if source.suffix.casefold() in {".json", ".gz"}:
        value = _read_json(source)
        if isinstance(value, list) and all(isinstance(item, Mapping) for item in value):
            return tuple(value)
        raise ValueError("evaluation JSON must contain a list of row objects")
    rows = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, Mapping):
            raise ValueError("evaluation JSONL rows must be objects")
        rows.append(value)
    return tuple(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate and assess products using graph-derived operators."
    )
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser(
        "build-library",
        help="project a serialized generic operator library forward",
    )
    build.add_argument("generic_library")
    build.add_argument("output")
    build.add_argument("--skip-source-round-trip", action="store_true")

    predict = commands.add_parser("predict", help="blindly predict products")
    predict.add_argument("library")
    predict.add_argument("starting_materials")
    predict.add_argument("--recipe")
    predict.add_argument("--top-k", type=int, default=20)

    assess = commands.add_parser(
        "assess-step",
        help="run targeted replay and blind competition for a route step",
    )
    assess.add_argument("library")
    assess.add_argument("starting_materials")
    assess.add_argument("intended_product")
    assess.add_argument("--recipe")
    assess.add_argument("--operator-hint")
    assess.add_argument("--top-k", type=int, default=20)

    evaluate = commands.add_parser(
        "evaluate",
        help="run reference-disjoint product-hidden replay",
    )
    evaluate.add_argument("library")
    evaluate.add_argument("rows")
    evaluate.add_argument("--top-k", type=int, default=20)
    evaluate.add_argument("--allow-reference-overlap", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    if arguments.command == "build-library":
        source = _read_json(arguments.generic_library)
        library = build_forward_library(
            source,
            require_source_round_trip=not arguments.skip_source_round_trip,
        )
        save_forward_library(library, arguments.output)
        print(json.dumps(library.to_dict(), indent=2, sort_keys=True))
        return 0
    library = load_forward_library(arguments.library)
    if arguments.command == "evaluate":
        result = evaluate_product_hidden_replay(
            _rows(arguments.rows),
            library,
            top_k=arguments.top_k,
            require_reference_disjoint=not arguments.allow_reference_overlap,
        )
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
        return 0
    recipe = _recipe(arguments.recipe)
    if arguments.command == "predict":
        result = predict_products(
            arguments.starting_materials,
            library,
            recipe=recipe,
            top_k=arguments.top_k,
        )
    else:
        result = assess_proposed_step(
            arguments.starting_materials,
            arguments.intended_product,
            library,
            recipe=recipe,
            operator_hint=arguments.operator_hint,
            top_k=arguments.top_k,
        )
    print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
    valid = (
        result.blind_prediction.valid
        if hasattr(result, "blind_prediction")
        else result.valid
    )
    return 0 if valid else 1


__all__ = ["build_parser", "main"]
