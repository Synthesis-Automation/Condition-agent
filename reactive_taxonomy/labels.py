"""Taxonomy-driven deterministic chemist-facing label rendering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List

from .notation import (
    available_styles,
    format_chemist_text,
    notation_style,
    render_context_notation,
    render_fragment_notation,
)


_RENDERING_PATH = Path(__file__).with_name("definitions") / "rendering.v1.json"


@lru_cache(maxsize=1)
def load_rendering_taxonomy() -> Dict[str, Any]:
    with _RENDERING_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _style(style: str) -> Dict[str, str]:
    styling = notation_style(style)
    return {
        "bond": styling["single"],
        "double": styling["double"],
        "triple": styling["triple"],
        "negative_charge": styling["negative_charge"],
        "positive_charge": styling["positive_charge"],
    }


def render_context(token: str, *, style: str = "unicode") -> str:
    return render_context_notation(token, style=style)


def _render_context_or_fragment(token: str, *, style: str) -> str:
    try:
        return render_context_notation(token, style=style)
    except ValueError:
        return render_fragment_notation(token, style=style)


def _render_contexts(tokens: List[str], *, style: str) -> List[str]:
    """Render classification contexts, excluding non-display activation facets."""
    rendered = []
    for token in tokens:
        try:
            rendered.append(render_context_notation(token, style=style))
        except ValueError:
            continue
    return rendered


def render_anion(
    center: str,
    contexts: List[str],
    charge: int,
    *,
    style: str = "unicode",
) -> str:
    """Render an explicitly charged heteroatom nucleophile."""
    styling = _style(style)
    bond = styling["bond"]
    rendered_contexts = _render_contexts(contexts, style=style)
    context = rendered_contexts[0] if rendered_contexts else ""
    magnitude = abs(int(charge))
    sign = styling["negative_charge"] if charge < 0 else styling["positive_charge"]
    charge_label = (str(magnitude) if magnitude > 1 else "") + sign
    prefix = f"{context}{bond}" if context else ""
    return format_chemist_text(f"{prefix}{center}{charge_label}", style=style)


def render_edge(context: str, handle: str, *, style: str = "unicode") -> str:
    bond = _style(style)["bond"]
    template = load_rendering_taxonomy().get("edge_template", "{context}{bond}{handle}")
    return format_chemist_text(
        str(template).format(
            context=_render_context_or_fragment(context, style=style),
            bond=bond,
            handle=handle,
        ),
        style=style,
    )


def render_named_handle(template_id: str, *, context: str = "", style: str = "unicode") -> str:
    """Render a taxonomy-owned label for a center-, atom-, or bond-site."""
    bond = _style(style)["bond"]
    templates = load_rendering_taxonomy().get("named_handle_templates", {})
    template = templates.get(template_id)
    if template is None:
        raise ValueError(f"Unknown named handle template: {template_id}")
    return format_chemist_text(str(template).format(
        bond=bond,
        double=_style(style)["double"],
        triple=_style(style)["triple"],
        context=render_context_notation(context, style=style) if context else "",
    ), style=style)


def render_unsaturated_bond(
    bond_order: int,
    endpoint_h_counts: List[int],
    endpoint_substituent_counts: List[int],
    stereochemistry: str = "",
    *,
    style: str = "unicode",
) -> str:
    """Render a C=C or C≡C site with explicit H/R substitution topology."""
    styling = _style(style)
    templates = load_rendering_taxonomy().get("unsaturated_bond_templates", {})
    if bond_order == 2:
        rules = templates["alkene"]
        h_left, h_right = (int(value) for value in endpoint_h_counts)
        substituent_left, substituent_right = (
            int(value) for value in endpoint_substituent_counts
        )
        if h_left + substituent_left != 2 or h_right + substituent_right != 2:
            return f"C{styling['double']}C"
        next_r = 1

        def endpoint(side: str, h_count: int, substituent_count: int) -> str:
            nonlocal next_r
            indices = list(range(next_r, next_r + substituent_count))
            next_r += substituent_count
            values = {
                "r1": indices[0] if indices else "",
                "r2": indices[1] if len(indices) > 1 else "",
            }
            return str(rules[f"{side}_endpoints"][str(h_count)]).format(**values)

        left = endpoint("left", h_left, substituent_left)
        right = endpoint("right", h_right, substituent_right)
        stereo_suffix = f" ({stereochemistry})" if stereochemistry in {"E", "Z"} else ""
        return format_chemist_text(str(rules["template"]).format(
            left=left,
            right=right,
            double=styling["double"],
            stereo_suffix=stereo_suffix,
        ), style=style)
    if bond_order == 3:
        rules = templates["alkyne"]
        substituent_count = sum(int(value) for value in endpoint_substituent_counts)
        key = "acetylene" if substituent_count == 0 else (
            "terminal" if substituent_count == 1 else "internal"
        )
        return format_chemist_text(str(rules[key]).format(
            bond=styling["bond"],
            triple=styling["triple"],
            r1=1,
            r2=2,
        ), style=style)
    raise ValueError(f"Unsupported unsaturated bond order: {bond_order}")


def _rule_matches(
    rule: Dict[str, Any],
    center: str,
    h_count: int,
    contexts: List[str],
    features: Dict[str, Any],
) -> bool:
    if rule.get("center") != center:
        return False
    if "h_count" in rule and int(rule["h_count"]) != h_count:
        return False
    if "contexts_exact" in rule and list(rule["contexts_exact"]) != contexts:
        return False
    if "contexts_set" in rule and sorted(rule["contexts_set"]) != sorted(contexts):
        return False
    if "contexts_contains" in rule and rule["contexts_contains"] not in contexts:
        return False
    if "context_first" in rule and (not contexts or contexts[0] != rule["context_first"]):
        return False
    if "context_first_prefix" in rule and (not contexts or not contexts[0].startswith(rule["context_first_prefix"])):
        return False
    if any(
        features.get(key) != expected
        for key, expected in (rule.get("feature_equals") or {}).items()
    ):
        return False
    return True


def render_xh(
    center: str,
    h_count: int,
    contexts: List[str],
    features: Dict[str, Any] | None = None,
    *,
    style: str = "unicode",
) -> str:
    styling = _style(style)
    bond = styling["bond"]
    feature_values = dict(features or {})
    rules = sorted(load_rendering_taxonomy().get("xh_rules", []), key=lambda item: -int(item.get("priority", 0)))
    rule = next(
        (
            item
            for item in rules
            if _rule_matches(item, center, h_count, contexts, feature_values)
        ),
        None,
    )
    if rule is None:
        return f"{center}{bond}H"
    rendered_contexts = _render_contexts(contexts, style=style)
    context = rendered_contexts[0] if rendered_contexts else "H"
    suffix = "2" if h_count == 2 else ("R" if h_count == 1 and len(contexts) > 1 else ("" if h_count == 1 else str(h_count)))
    return format_chemist_text(str(rule["template"]).format(
        bond=bond,
        context=context,
        contexts=bond.join(rendered_contexts) if rendered_contexts else "H",
        suffix=suffix,
        triple=styling["triple"],
    ), style=style)


__all__ = ["available_styles", "load_rendering_taxonomy", "render_context", "render_edge", "render_named_handle", "render_unsaturated_bond", "render_xh"]
