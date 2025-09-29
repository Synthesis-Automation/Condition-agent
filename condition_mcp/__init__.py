"""Condition MCP package scaffolding.

This package collects the thin wrappers, schemas, and resources used by the
upcoming Model Context Protocol (MCP) façade around the condition recommendation
stack. The goal is to keep these utilities importable without requiring the
full MCP transport layer (which will be added later in the implementation
phases).

"""

from __future__ import annotations

from importlib import metadata

DEFAULT_SCHEMA_VERSION = "0.1.0"

try:
    __version__ = metadata.version("condition_mcp")
except metadata.PackageNotFoundError:  # pragma: no cover - local editable installs
    __version__ = "0.0.0"

__all__ = ["DEFAULT_SCHEMA_VERSION", "__version__"]
