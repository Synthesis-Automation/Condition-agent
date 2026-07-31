"""CLI for generic type-agnostic reaction-dataset conversion."""

from __future__ import annotations

import argparse
import json

from reactive_taxonomy import RxnMapperProvider

from .conversion.engine import convert_datasets


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert mixed reaction CSVs using chemistry signatures"
    )
    parser.add_argument("dataset_path", help="CSV file or directory of CSV files")
    parser.add_argument("output_dir", help="Directory for converted artifacts")
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional file-order row limit for a conversion smoke test",
    )
    parser.add_argument(
        "--use-rxnmapper",
        action="store_true",
        help=(
            "Use optional RXNMapper evidence for unresolved/ambiguous reactions "
            "and consensus-only core enrichment for resolved reactions"
        ),
    )
    args = parser.parse_args()
    if args.use_rxnmapper and not RxnMapperProvider.is_available():
        parser.error(
            "RXNMapper is not installed; run "
            "'python -m pip install -r requirements-mapping.txt'"
        )
    mapping_provider = RxnMapperProvider() if args.use_rxnmapper else None
    report = convert_datasets(
        args.dataset_path,
        args.output_dir,
        max_rows=args.max_rows,
        mapping_provider=mapping_provider,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
