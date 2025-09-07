#!/usr/bin/env python3
"""
Thin wrapper script for the registry resolver CLI.
Delegates to chemtools.cli.registry:main for reuse in console_scripts.
"""
from chemtools.cli.registry import main  # type: ignore

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
