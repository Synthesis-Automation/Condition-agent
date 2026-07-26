"""CLI for development-only generic ranking calibration."""

from __future__ import annotations

import argparse
import json

from .calibration import calibrate_generic_ranking


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate ranking on development and gate on validation"
    )
    parser.add_argument("development_index_path")
    parser.add_argument("validation_index_path")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    report = calibrate_generic_ranking(
        args.development_index_path,
        args.validation_index_path,
        args.output_dir,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
