"""Model Context Protocol integration utilities for chemtools."""
from __future__ import annotations

from importlib import metadata

DEFAULT_SCHEMA_VERSION = "0.1.0"

try:  # pragma: no cover - package metadata optional during development
    __version__ = metadata.version("chemtools")
except metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

from . import client, tools  # noqa: E402  (re-export for convenience)

__all__ = ["DEFAULT_SCHEMA_VERSION", "__version__", "client", "tools"]
