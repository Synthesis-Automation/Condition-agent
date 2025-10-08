#!/usr/bin/env python3
"""
Thin wrapper script for the registry resolver CLI.
Delegates to chemtools.rule_scdb_matcher.cli:main for reuse in console_scripts.
"""
from chemtools.rule_scdb_matcher.cli import main  # type: ignore

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
