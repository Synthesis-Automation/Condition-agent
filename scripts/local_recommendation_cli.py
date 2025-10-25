"""
Compatibility shim for the historical ``scripts.local_recommendation_cli`` path.

The implementation now lives in ``app.local_recommendation_cli``.  This module
simply forwards attribute access lazily to avoid side effects at import time.
"""

from __future__ import annotations

import importlib
from typing import Any


_MODULE_NAME = "app.local_recommendation_cli"


def _load_module():
    return importlib.import_module(_MODULE_NAME)


def __getattr__(name: str) -> Any:
    module = _load_module()
    return getattr(module, name)


def __dir__() -> list[str]:
    module = _load_module()
    return sorted(set(dir(module)))
