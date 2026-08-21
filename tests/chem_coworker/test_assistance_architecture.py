"""Architecture boundaries for the assistance application layer."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ASSISTANCE = ROOT / "chem_coworker" / "assistance"


def test_assistance_has_no_dynamic_imports_or_chemistry_implementation() -> None:
    forbidden_calls = {
        "Chem.MolFromSmarts",
        "importlib.import_module",
        "__import__",
    }
    for path in ASSISTANCE.glob("*.py"):
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        assert "chemtools" not in source
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if isinstance(node.func, ast.Name):
                name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                parts = []
                current = node.func
                while isinstance(current, ast.Attribute):
                    parts.append(current.attr)
                    current = current.value
                if isinstance(current, ast.Name):
                    parts.append(current.id)
                name = ".".join(reversed(parts))
            else:
                continue
            assert name not in forbidden_calls


def test_assistance_does_not_create_generic_tool_or_skill_directories() -> None:
    assert not (ROOT / "chem_coworker" / "tools").exists()
    assert not (ROOT / "chem_coworker" / "skills").exists()


def test_provider_request_code_exists_only_in_shared_transport() -> None:
    legacy_modules = (
        ROOT / "chem_coworker" / "review.py",
        ROOT / "chem_coworker" / "retrosynthesis_review.py",
        ROOT / "chem_coworker" / "multistep_review.py",
    )
    for path in legacy_modules:
        source = path.read_text(encoding="utf-8")
        assert "client.responses.create" not in source
        assert "client.chat.completions.create" not in source
        assert "def _legacy_complete" not in source
