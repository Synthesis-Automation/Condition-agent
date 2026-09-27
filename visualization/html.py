"""HTML embedding for trusted, locally generated molecular SVG drawings."""

from html import escape
import math
import re
import xml.etree.ElementTree as ET


def svg_html(svg: str, *, label: str = "Chemical structure") -> str:
    """Keep an SVG's native scale inside a keyboard-accessible scroll region.

    This embeds trusted renderer output, not arbitrary or sanitized user HTML.
    Inline dimensions prevent report styles from stretching individual drawings.
    """
    svg = svg[svg.index("<svg"):]
    root = ET.fromstring(svg)
    dimensions = []
    for key in ("width", "height"):
        value = float(root.attrib[key].removesuffix("px"))
        if not math.isfinite(value) or value <= 0:
            raise ValueError("SVG requires positive pixel dimensions")
        dimensions.append(value)
    width, height = dimensions
    opening, content = svg.split(">", 1)
    opening = re.sub(r"\sstyle=(['\"]).*?\1", "", opening, flags=re.S)
    style = (
        f"display:block;width:{width:g}px;height:{height:g}px;"
        "min-width:0;min-height:0;max-width:none;max-height:none;flex:none"
    )
    return (
        '<div class="chemical-svg-scroll" tabindex="0" role="region" '
        f'aria-label="{escape(label, quote=True)}" '
        'style="max-width:100%;min-width:0;overflow:auto;overscroll-behavior:contain">'
        f'{opening} style="{style}">{content}</div>'
    )
