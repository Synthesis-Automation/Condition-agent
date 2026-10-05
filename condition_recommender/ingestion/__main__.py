"""CLI for reproducible preparation of the consolidated raw dataset tree."""

from __future__ import annotations

import argparse
import json
import time

from .artifacts import PreprocessingProgress
from .migration import regenerate_intermediate_datasets, validate_intermediate_manifest


def main() -> None:
    """Regenerate the corpus or validate its prepared artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("regenerate", "validate"))
    parser.add_argument("--raw-root", default="raw_datasets")
    parser.add_argument("--output-root", default="datasets/intermediate_datasets")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    last = 0.0

    def progress(event: PreprocessingProgress) -> None:
        nonlocal last
        now = time.monotonic()
        if event.phase in {"started", "completed", "reused"} or now - last >= 20:
            print(event.message, flush=True)
            last = now

    if args.command == "validate":
        report = validate_intermediate_manifest(args.output_root)
    else:
        report = regenerate_intermediate_datasets(
            args.raw_root,
            args.output_root,
            force=args.force,
            progress_callback=progress,
        )
        report = {key: value for key, value in report.items() if key != "sources"}
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
