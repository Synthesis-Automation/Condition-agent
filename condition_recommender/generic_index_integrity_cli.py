"""CLI for deep validation of a persisted generic reaction index."""

from __future__ import annotations

import argparse
import json

from .generic_indexing import validate_generic_index_artifact


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate a persisted generic reaction index"
    )
    parser.add_argument("index_path")
    parser.add_argument("--output-path")
    args = parser.parse_args()
    report = validate_generic_index_artifact(args.index_path)
    if args.output_path:
        with open(args.output_path, "w", encoding="utf-8") as handle:
            json.dump(report, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
