from __future__ import annotations

import argparse
import json

from .conversion.heteroatom import convert_file


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert C–O or C–S coupling data")
    parser.add_argument("input_csv"); parser.add_argument("output_dir"); parser.add_argument("--element", required=True, choices=("O", "S"))
    args = parser.parse_args(); print(json.dumps(convert_file(args.input_csv, args.output_dir, element=args.element), indent=2))


if __name__ == "__main__": main()
