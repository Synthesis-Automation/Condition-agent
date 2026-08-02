"""Canonical chemist notation shared by molecular and reaction rendering."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping


_PATH = Path(__file__).with_name("definitions") / "chemist_notation.v1.json"
_SUPERSCRIPT_DIGITS = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
_SUBSCRIPT_DIGITS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


@lru_cache(maxsize=1)
def load_chemist_notation() -> dict[str, Any]:
    """Load the versioned canonical chemist-notation definition."""
    with _PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != "1.0":
        raise ValueError("Unsupported chemist-notation schema")
    return dict(payload)


def available_styles() -> tuple[str, ...]:
    """Return supported presentation styles in definition order."""
    return tuple(load_chemist_notation()["styles"])


def notation_style(style: str) -> dict[str, str]:
    """Return one validated notation style."""
    styles = load_chemist_notation()["styles"]
    if style not in styles:
        raise ValueError(f"Unknown chemist-notation style: {style}")
    return dict(styles[style])


@lru_cache(maxsize=1)
def _fragment_symbols() -> dict[str, str]:
    return {
        str(record["id"]): str(record["symbol"])
        for record in load_chemist_notation()["fragment_notations"]
    }


def render_fragment_notation(notation_id: str, *, style: str = "unicode") -> str:
    """Render one registered fragment class."""
    notation_style(style)
    try:
        symbol = _fragment_symbols()[notation_id]
    except KeyError as exc:
        raise ValueError(f"Unknown fragment notation: {notation_id}") from exc
    return format_chemist_text(symbol, style=style)


def _template_values(style: str) -> dict[str, str]:
    styling = notation_style(style)
    return {
        **styling,
        **{
            notation_id: render_fragment_notation(notation_id, style=style)
            for notation_id in _fragment_symbols()
        },
    }


def render_context_notation(context_id: str, *, style: str = "unicode") -> str:
    """Render one registered context classification."""
    templates: Mapping[str, str] = load_chemist_notation()["context_notations"]
    if context_id not in templates:
        raise ValueError(f"Unknown context notation: {context_id}")
    rendered = str(templates[context_id]).format(**_template_values(style))
    return format_chemist_text(rendered, style=style)


def render_remote_class(remote_class: str, *, style: str = "unicode") -> str:
    """Render one reaction-core remote class through the shared registry."""
    bindings: Mapping[str, str] = load_chemist_notation()[
        "remote_class_notations"
    ]
    if remote_class not in bindings:
        raise ValueError(f"Unknown remote fragment class: {remote_class}")
    return render_fragment_notation(bindings[remote_class], style=style)


def format_chemist_text(text: str, *, style: str = "unicode") -> str:
    """Apply canonical formula and fragment-index typography."""
    notation_style(style)
    if style != "unicode":
        return text
    symbols = tuple(load_chemist_notation()["fragment_index_symbols"])
    symbol_pattern = "|".join(
        re.escape(symbol) for symbol in sorted(symbols, key=lambda value: -len(value))
    )
    value = re.sub(
        rf"({symbol_pattern})([0-9]+)",
        lambda match: match.group(1)
        + match.group(2).translate(_SUPERSCRIPT_DIGITS),
        text,
    )
    return re.sub(
        r"([A-Z][a-z]?|\))([0-9]+)",
        lambda match: match.group(1)
        + match.group(2).translate(_SUBSCRIPT_DIGITS),
        value,
    )


__all__ = [
    "available_styles",
    "format_chemist_text",
    "load_chemist_notation",
    "notation_style",
    "render_context_notation",
    "render_fragment_notation",
    "render_remote_class",
]
