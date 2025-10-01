"""Common helpers for Condition MCP tool wrappers."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List

from pydantic import BaseModel, Field, ValidationError

from chemtools.integrations.mcp import DEFAULT_SCHEMA_VERSION


class ToolError(ValueError):
    """Raised when tool input validation fails."""


def validate_payload(model: type[BaseModel], data: Dict[str, Any]) -> BaseModel:
    """Validate ``data`` against ``model`` and surface consistent errors."""

    try:
        return model.model_validate(data)
    except ValidationError as exc:  # pragma: no cover - exercised via tests
        raise ToolError(str(exc)) from exc


def pick_first(values: Iterable[str | None]) -> str | None:
    """Return the first non-empty string from ``values``."""

    for value in values:
        if isinstance(value, str) and value.strip():
            return value
    return None


class SchemaStamped(BaseModel):
    """Base model that injects the shared schema version field."""

    schema_version: str = Field(default=DEFAULT_SCHEMA_VERSION)

    model_config = {
        "populate_by_name": True,
        "extra": "ignore",
    }


def flatten_smiles_block(items: Iterable[Dict[str, Any]]) -> List[str]:
    """Return canonical SMILES strings from normalized blocks."""

    out: List[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        choice = pick_first(
            (
                item.get("smiles_norm"),
                item.get("largest_smiles"),
                item.get("input"),
            )
        )
        if choice:
            out.append(choice)
    return out
