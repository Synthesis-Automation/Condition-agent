"""Run the chat controller's DOM-independent interaction regressions with Node."""

from pathlib import Path
import shutil
import subprocess

import pytest


def test_scientific_chat_controller() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node is needed for chat controller tests")
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [node, "--test", str(root / "tests" / "web" / "scientific_chat.test.cjs")],
        cwd=root, capture_output=True, text=True, encoding="utf-8", timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
