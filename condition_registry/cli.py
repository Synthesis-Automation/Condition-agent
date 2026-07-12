from __future__ import annotations

import argparse
import json

from .migration import run_migration_audit


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit the standalone condition registry")
    parser.add_argument("output_dir"); args = parser.parse_args()
    print(json.dumps(run_migration_audit(args.output_dir), indent=2))


if __name__ == "__main__": main()
