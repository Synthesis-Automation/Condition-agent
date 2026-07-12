"""Taxonomy-driven chemist-friendly reaction rendering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict

from .labels import render_context, render_xh
from .reaction_models import ReactionSiteReference


_PATH = Path(__file__).with_name("definitions") / "reaction_rendering.v1.json"
_PRECEDENCE = {"SO2R": 0, "C(O)NR": 1, "C(O)OR": 2, "C(O)R": 3, "HeteroAr": 4, "Ar": 5, "Alkenyl": 6, "Alkyl": 7, "N": 8}


@lru_cache(maxsize=1)
def load_reaction_rendering() -> Dict[str, Dict[str, Any]]:
    with _PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle).get("rules") or {}


def _nitrogen_product(anchor_context: str, nitrogen: ReactionSiteReference, bond: str, style: str) -> str:
    if nitrogen.details.get("center_token") == "N_aromatic":
        return f"{render_context(anchor_context, style=style)}{bond}AromN"
    contexts = [anchor_context, *list(nitrogen.details.get("contexts") or [])]
    contexts.sort(key=lambda token: (_PRECEDENCE.get(token, 99), token))
    return render_xh("N", max(0, int(nitrogen.details["h_count"]) - 1), contexts, style=style)


def render_reactant_label(assignment: Dict[str, ReactionSiteReference], *, style: str = "unicode") -> str:
    """Render only the observed reactant handles for a candidate assignment."""
    bond = "–" if style == "unicode" else "-"
    return " + ".join(site.chemist_label.replace("–", bond) for site in assignment.values())


def render_reaction_label(grammar: Dict[str, Any], assignment: Dict[str, ReactionSiteReference], *, style: str = "unicode") -> str:
    bond = "–" if style == "unicode" else "-"
    reactant_labels = render_reactant_label(assignment, style=style)
    rule = load_reaction_rendering().get(str(grammar["id"]), {})
    kind = rule.get("product_kind")
    if kind == "join_contexts":
        left, right = assignment[rule["left_role"]], assignment[rule["right_role"]]
        product = f"{render_context(left.details['anchor_context'], style=style)}{bond}{render_context(right.details['anchor_context'], style=style)}"
    elif kind == "nitrogen_substitution":
        product = _nitrogen_product(assignment[rule["anchor_role"]].details["anchor_context"], assignment[rule["nitrogen_role"]], bond, style)
    elif kind == "heteroatom_substitution":
        anchor, partner = assignment[rule["anchor_role"]], assignment[rule["partner_role"]]
        retained = list(partner.details.get("contexts") or [])
        right = render_context(retained[0], style=style) if retained else "H"
        product = f"{render_context(anchor.details['anchor_context'], style=style)}{bond}{rule['element']}{bond}{right}"
    elif kind == "terminal_alkyne":
        anchor = assignment[rule["anchor_role"]]
        product = f"{render_context(anchor.details['anchor_context'], style=style)}{bond}C≡C{bond}R"
    elif kind == "chan_lam":
        anchor, partner = assignment[rule["anchor_role"]], assignment[rule["partner_role"]]
        center = partner.details.get("center_token")
        if center in {"N", "N_aromatic"}:
            product = _nitrogen_product(anchor.details["anchor_context"], partner, bond, style)
        else:
            retained = list(partner.details.get("contexts") or [])
            right = render_context(retained[0], style=style) if retained else "H"
            product = f"{render_context(anchor.details['anchor_context'], style=style)}{bond}{center}{bond}{right}"
    elif kind == "amide":
        nitrogen = assignment[rule["nitrogen_role"]]
        product = render_xh("N", max(0, int(nitrogen.details["h_count"]) - 1), ["C(O)R", *list(nitrogen.details.get("contexts") or [])], style=style)
    elif kind == "acyl_heteroatom":
        acyl, partner = assignment[rule["acyl_role"]], assignment[rule["partner_role"]]
        retained = list(partner.details.get("contexts") or [])
        right = render_context(retained[0], style=style) if retained else "H"
        product = f"{render_context(acyl.details['retained_context'], style=style)}{bond}C(O){bond}{rule['element']}{bond}{right}"
    elif kind == "sulfonamide":
        sulfonyl, nitrogen = assignment[rule["sulfonyl_role"]], assignment[rule["nitrogen_role"]]
        suffix = "2" if int(nitrogen.details["h_count"]) == 2 else "R"
        product = f"{render_context(sulfonyl.details['retained_context'], style=style)}{bond}S(O)2{bond}NH{suffix}"
    else:
        product = str(grammar.get("generic_label") or grammar["id"])
    return f"{reactant_labels} → {product}"


__all__ = ["load_reaction_rendering", "render_reactant_label", "render_reaction_label"]
