"""Render a provider-backed route workbench JSON report as SVG HTML."""

from __future__ import annotations

import argparse
import json
from typing import Sequence

from core_retrosynthesis.route_workbench_review import (
    write_provider_route_workbench_review_html,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Render one report and print artifact metadata."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_report")
    parser.add_argument("output_html")
    parser.add_argument(
        "--title",
        default="Provider-backed retrosynthesis review",
    )
    arguments = parser.parse_args(argv)
    result = write_provider_route_workbench_review_html(
        arguments.source_report,
        arguments.output_html,
        title=arguments.title,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
