from __future__ import annotations

import argparse
import json

from .coverage import DEFAULT_DATASETS, analyze_coverage


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit condition registry coverage across pilot datasets")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    print(json.dumps(analyze_coverage(DEFAULT_DATASETS, args.output_dir), indent=2))


if __name__ == "__main__": main()
