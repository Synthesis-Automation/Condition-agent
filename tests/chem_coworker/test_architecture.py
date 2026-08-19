"""Architecture boundaries for the rebuilt coworker."""

from __future__ import annotations

import ast
import importlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_legacy_packages_are_removed() -> None:
    assert not (ROOT / "chemtools").exists()
    assert not (ROOT / "chem_coworker" / "tools").exists()
    assert not (ROOT / "chem_coworker" / "skills").exists()


def test_coworker_has_no_legacy_imports() -> None:
    violations = []
    for path in (ROOT / "chem_coworker").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            if any(name == "chemtools" or name.startswith("chemtools.") for name in names):
                violations.append(str(path.relative_to(ROOT)))
    assert violations == []


def test_import_does_not_eagerly_load_legacy_or_llm_modules() -> None:
    before = set(sys.modules)
    importlib.import_module("chem_coworker")
    newly_loaded = set(sys.modules) - before

    assert not any(name == "chemtools" or name.startswith("chemtools.") for name in newly_loaded)
    assert not any(name == "llmtools" or name.startswith("llmtools.") for name in newly_loaded)
