"""CLI for one-file-to-one-file chemistry-free source preprocessing."""

from __future__ import annotations

import argparse
import json

from .ingestion import adapter_ids, preprocess_files


def main() -> None:
    """Preprocess explicitly selected source files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_files", nargs="+", help="Source CSV file(s)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--adapter",
        default="auto",
        choices=("auto", *adapter_ids()),
        help="Source adapter; auto requires one unambiguous schema match",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild even when the source checksum and adapter are unchanged",
    )
    args = parser.parse_args()
    report = preprocess_files(
        args.source_files,
        args.output_dir,
        adapter_id=args.adapter,
        force=args.force,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
