"""Architecture boundaries for the rebuilt coworker."""

from __future__ import annotations

import ast
import importlib
import importlib.util
import subprocess
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


def test_scientific_execution_does_not_depend_on_agent_runtime() -> None:
    """Scientific process control must not bypass identity through a runtime import."""
    prefix = "chem_coworker.scientific_workspace"
    violations = []
    workspace = ROOT / "chem_coworker" / "scientific_workspace"
    for folder in ("core", "adapters"):
        for path in (workspace / folder).glob("*.py"):
            package = f"{prefix}.{folder}"
            for node in ast.walk(ast.parse(path.read_text("utf-8"))):
                if isinstance(node, ast.ImportFrom):
                    name = importlib.util.resolve_name(
                        "." * node.level + (node.module or ""), package,
                    ) if node.level else node.module or ""
                    names = [name]
                    if name == prefix:
                        names.extend(f"{name}.{alias.name}" for alias in node.names)
                elif isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                else:
                    continue
                if any(name == f"{prefix}.runtime" or name.startswith(f"{prefix}.runtime.")
                       for name in names):
                    violations.append(str(path.relative_to(ROOT)))
    assert violations == []


def test_headless_execution_import_keeps_agent_harness_unloaded() -> None:
    """Check a fresh interpreter so earlier conversation tests cannot mask coupling."""
    subprocess.run(
        [sys.executable, "-c", (
            "import sys; import chem_coworker.scientific_workspace.core.execution; "
            "assert not any(name.startswith('chem_coworker.scientific_workspace.runtime') "
            "for name in sys.modules)"
        )],
        cwd=ROOT, check=True, capture_output=True, text=True, timeout=30,
    )
