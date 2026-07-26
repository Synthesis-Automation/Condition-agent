"""CLI for validating a sharded generic conversion artifact."""

from __future__ import annotations

import argparse
import json

from .conversion.sharded import validate_sharded_conversion


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate generic conversion shard provenance and integrity"
    )
    parser.add_argument("manifest_path")
    parser.add_argument(
        "--checksums-only",
        action="store_true",
        help="Skip decompressed row and observation-ID validation",
    )
    args = parser.parse_args()
    report = validate_sharded_conversion(
        args.manifest_path,
        verify_rows=not args.checksums_only,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
