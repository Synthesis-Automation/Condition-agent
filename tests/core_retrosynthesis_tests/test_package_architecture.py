"""Regression tests for the canonical retrosynthesis package boundary."""

from __future__ import annotations

from pathlib import Path

import core_retrosynthesis
from core_retrosynthesis.baselines import cx_rdchiral


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_only_canonical_top_level_retrosynthesis_package_exists() -> None:
    assert core_retrosynthesis.__name__ == "core_retrosynthesis"
    assert not (PROJECT_ROOT / "retrosynthesis_poc").exists()
    assert not (PROJECT_ROOT / "core_retrosynthesis_poc").exists()


def test_canonical_package_has_no_legacy_imports() -> None:
    package_root = PROJECT_ROOT / "core_retrosynthesis"
    offenders = []
    for path in package_root.rglob("*.py"):
        source = path.read_text(encoding="utf-8")
        if (
            "from retrosynthesis_poc" in source
            or "import retrosynthesis_poc" in source
        ):
            offenders.append(path.relative_to(PROJECT_ROOT).as_posix())
    assert offenders == []


def test_historical_rdchiral_code_is_explicitly_a_baseline() -> None:
    assert cx_rdchiral.__name__ == (
        "core_retrosynthesis.baselines.cx_rdchiral"
    )
