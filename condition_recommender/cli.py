"""Command-line entry points for recommendation data preparation."""

from __future__ import annotations

import argparse
import json

from .conversion.suzuki import convert_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert Suzuki data using reactive_taxonomy features")
    parser.add_argument("input_csv")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    print(json.dumps(convert_file(args.input_csv, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
