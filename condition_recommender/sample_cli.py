"""CLI for deterministic reference-safe pilot dataset sampling."""

from __future__ import annotations

import argparse
import json

from .conversion.sampling import build_reference_safe_samples


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build reference-safe pilot samples from reaction CSVs"
    )
    parser.add_argument("dataset_path", help="CSV file or directory of CSV files")
    parser.add_argument("output_dir", help="Destination for samples and manifest")
    parser.add_argument("--smoke-size", type=int, default=500)
    parser.add_argument("--development-size", type=int, default=5000)
    parser.add_argument("--validation-size", type=int, default=2000)
    parser.add_argument("--test-size", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--max-rows-per-group", type=int, default=20)
    args = parser.parse_args()
    report = build_reference_safe_samples(
        args.dataset_path,
        args.output_dir,
        smoke_size=args.smoke_size,
        development_size=args.development_size,
        validation_size=args.validation_size,
        test_size=args.test_size,
        seed=args.seed,
        max_rows_per_group=args.max_rows_per_group,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
