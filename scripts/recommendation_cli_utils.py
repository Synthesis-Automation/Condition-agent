"""
Compatibility shim for ``scripts.recommendation_cli_utils``.

The utilities were relocated under ``app``; this module re-exports them lazily.
"""

from __future__ import annotations

import importlib
from typing import Any


_MODULE_NAME = "app.recommendation_cli_utils"


def _load_module():
    return importlib.import_module(_MODULE_NAME)


def __getattr__(name: str) -> Any:
    module = _load_module()
    return getattr(module, name)


def __dir__() -> list[str]:
    module = _load_module()
    return sorted(set(dir(module)))
