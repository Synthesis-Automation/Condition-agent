#!/usr/bin/env python3
"""Run the standard checklist for substituent fragment edits."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: list[str]) -> int:
    print(f"> {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, cwd=REPO_ROOT)
    return int(proc.returncode)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Checklist for substituent_fragments.v1.json changes."
    )
    parser.add_argument(
        "--full-tests",
        action="store_true",
        help="Run full pytest suite after focused substituent tests.",
    )
    args = parser.parse_args()

    commands: list[list[str]] = [
        [sys.executable, "-m", "chemtools.taxonomy.validate_and_sync", "--check-only"],
        [sys.executable, "-m", "pytest", "-q", "tests/test_substituent_composer.py"],
    ]
    if args.full_tests:
        commands.append([sys.executable, "-m", "pytest", "-q"])

    for cmd in commands:
        rc = _run(cmd)
        if rc != 0:
            print(f"Checklist failed with exit code {rc}.")
            return rc

    print("Checklist passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
