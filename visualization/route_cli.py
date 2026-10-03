"""Render a route JSON file with ``python -m visualization.route_cli``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from .annotated_scheme import SchemeAnnotation
from .route_scheme import RouteSchemeStep, render_route_scheme_svg


def _annotation(value: Any) -> SchemeAnnotation:
    if isinstance(value, str):
        return SchemeAnnotation(value, "supplied")
    if (
        not isinstance(value, dict)
        or set(value) - {"text", "basis"}
        or not isinstance(value.get("text"), str)
        or not isinstance(value.get("basis", "supplied"), str)
    ):
        raise ValueError("Annotations require text and optional basis")
    return SchemeAnnotation(value["text"], value.get("basis", "supplied"))


def _step(value: Any) -> str | RouteSchemeStep:
    if isinstance(value, str):
        return value
    if not isinstance(value, dict) or set(value) - set(
        RouteSchemeStep.__dataclass_fields__
    ):
        raise ValueError("Unknown route step fields")
    if not isinstance(value.get("reaction_smiles"), str) or not isinstance(
        value.get("basis", "supplied"), str
    ):
        raise ValueError(
            "Each step requires reaction_smiles and an optional text basis"
        )
    for side in ("above", "below", "reagents"):
        if not isinstance(value.get(side, []), list):
            raise ValueError(f"{side} must be an annotation list")
    return RouteSchemeStep(
        reaction_smiles=value["reaction_smiles"],
        above=tuple(_annotation(item) for item in value.get("above", [])),
        below=tuple(_annotation(item) for item in value.get("below", [])),
        yield_info=_annotation(value["yield_info"])
        if value.get("yield_info") is not None
        else None,
        basis=value.get("basis", "supplied"),
        product_index=value.get("product_index"),
        reagents=tuple(_annotation(item) for item in value.get("reagents", [])),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Write one SVG from an ordered JSON list or annotated route object."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="UTF-8 route JSON")
    parser.add_argument("output", type=Path, help="Output SVG path")
    args = parser.parse_args(argv)
    try:
        payload = json.loads(args.input.read_text(encoding="utf-8-sig"))
        if isinstance(payload, list):
            payload = {"steps": payload}
        if not isinstance(payload, dict) or set(payload) - {
            "steps",
            "title",
            "main_reactant_index",
            "reagents_only",
        }:
            raise ValueError(
                "Expected a route object with steps, title, main_reactant_index, reagents_only"
            )
        if not isinstance(payload.get("steps"), list) or not isinstance(
            payload.get("title", "Reaction route"), str
        ):
            raise ValueError("Route requires a steps list and an optional text title")
        if type(payload.get("reagents_only", False)) is not bool:
            raise ValueError("reagents_only must be a boolean")
        svg = render_route_scheme_svg(
            tuple(_step(value) for value in payload["steps"]),
            title=payload.get("title", "Reaction route"),
            main_reactant_index=payload.get("main_reactant_index"),
            reagents_only=payload.get("reagents_only", False),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(svg)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
