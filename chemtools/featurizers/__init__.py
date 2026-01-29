from __future__ import annotations

import importlib
from typing import List

__all__ = [
    "analysis",
    "calculable",
    "molecule",
    "reaction_pair",
    "unified",
]


def __getattr__(name: str):
    if name in __all__:
        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> List[str]:
    return sorted(list(globals().keys()) + __all__)
