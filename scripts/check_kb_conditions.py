#!/usr/bin/env python3
"""Automated checks for KB condition summaries."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


def _run(cmd: list[str], label: str) -> bool:
    print(f"== {label} ==")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout)
    if proc.returncode != 0:
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
        print(f"{label} failed (exit {proc.returncode}).", file=sys.stderr)
        return False
    return True


def _run_pytest() -> bool:
    return _run(
        [sys.executable, "-m", "pytest", "-q", "tests/test_kb_conditions.py"],
        "pytest test_kb_conditions",
    )


def _run_kb_tool_check() -> bool:
    print("== kb_recommend_conditions sanity check ==")
    try:
        repo_root = Path(__file__).resolve().parents[1]
        if str(repo_root) not in sys.path:
            sys.path.insert(0, str(repo_root))
        from chem_assistant.chemtools_wrapper import kb_recommend_conditions

        payload = kb_recommend_conditions.invoke(
            {"query": "5-membered heterocyclic boronates", "top_k": 3}
        )
    except Exception as exc:
        print(f"kb_recommend_conditions failed: {exc}", file=sys.stderr)
        return False

    print(json.dumps(payload, indent=2))
    results = payload.get("results") or []
    if not results:
        print("No KB condition results returned.", file=sys.stderr)
        return False
    return True


def main() -> int:
    ok = True
    ok = _run_pytest() and ok
    ok = _run_kb_tool_check() and ok
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
