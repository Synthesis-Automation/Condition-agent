"""Taxonomy-driven chemist-friendly reaction rendering."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict

from .labels import render_context, render_xh
from .reaction_models import ReactionSiteReference


_PATH = Path(__file__).with_name("definitions") / "reaction_rendering.v1.json"
_PRECEDENCE = {
    "SO2R": 0,
    "C(O)NR": 1,
    "C(O)OR": 2,
    "C(O)R": 3,
    "HeteroAr": 4,
    "Ar": 5,
    "Alkenyl": 6,
    "Alkyl": 7,
    "N": 8,
}


@lru_cache(maxsize=1)
def load_reaction_rendering() -> Dict[str, Dict[str, Any]]:
    with _PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle).get("rules") or {}


@lru_cache(maxsize=1)
def load_product_context_precedence() -> tuple[str, ...]:
    with _PATH.open("r", encoding="utf-8") as handle:
        return tuple(json.load(handle).get("product_context_precedence") or ())


def render_context_connection(
    left_context: str, right_context: str, *, style: str = "unicode"
) -> str:
    """Render a symmetric C–C connection in canonical display order."""
    bond = "–" if style == "unicode" else "-"
    precedence = {
        token: index for index, token in enumerate(load_product_context_precedence())
    }
    contexts = sorted(
        (left_context, right_context),
        key=lambda token: (precedence.get(token, 999), token),
    )
    return bond.join(render_context(context, style=style) for context in contexts)


def _nitrogen_product(
    anchor_context: str, nitrogen: ReactionSiteReference, bond: str, style: str
) -> str:
    if nitrogen.details.get("center_token") == "N_aromatic":
        return f"{render_context(anchor_context, style=style)}{bond}AromN"
    anchor = render_context(anchor_context, style=style)

    def n_substituent(token: str) -> str:
        attachment_oriented = {
            "C(O)R": f"C(O){bond}R",
            "C(O)OR": f"C(O){bond}OR",
            "C(O)NR": f"C(O){bond}NR2",
            "SO2R": f"S(O)2{bond}R",
            "P(O)R2": "P(O)R2",
        }
        return attachment_oriented.get(token, render_context(token, style=style))

    retained = [
        n_substituent(token) for token in nitrogen.details.get("contexts") or []
    ]
    remaining_h = max(0, int(nitrogen.details["h_count"]) - 1)
    if remaining_h >= 2:
        return f"{anchor}{bond}NH2"
    if remaining_h == 1:
        return f"{anchor}{bond}NH" + (f"{bond}{retained[0]}" if retained else "")
    if not retained:
        return f"{anchor}{bond}N"
    if len(retained) == 2 and retained[0] == retained[1] == "R":
        return f"{anchor}{bond}NR1R2"
    return f"{anchor}{bond}N" + "".join(f"({context})" for context in retained)


def render_reactant_label(
    assignment: Dict[str, ReactionSiteReference], *, style: str = "unicode"
) -> str:
    """Render only the observed reactant handles for a candidate assignment."""
    bond = "–" if style == "unicode" else "-"
    return " + ".join(
        site.chemist_label.replace("–", bond) for site in assignment.values()
    )


def render_heteroatom_product(
    anchor_context: str,
    partner: ReactionSiteReference,
    element: str,
    *,
    style: str = "unicode",
) -> str:
    """Render a C–O/C–S product from anchor and retained partner context."""
    bond = "–" if style == "unicode" else "-"
    retained = list(partner.details.get("contexts") or [])
    right = render_context(retained[0], style=style) if retained else "H"
    return f"{render_context(anchor_context, style=style)}{bond}{element}{bond}{right}"


def _ordered_context_labels(site: ReactionSiteReference, *, style: str) -> list[str]:
    tokens = sorted(
        (str(token) for token in site.details.get("contexts") or ()),
        key=lambda token: (_PRECEDENCE.get(token, 999), token),
    )
    return [render_context(token, style=style) for token in tokens]


def _carbonyl_fragment(
    carbonyl: ReactionSiteReference, *, reduced: bool, style: str
) -> str:
    bond = "–" if style == "unicode" else "-"
    subtype = str(carbonyl.details.get("carbonyl_subtype") or "ketone")
    contexts = _ordered_context_labels(carbonyl, style=style)
    if subtype == "formaldehyde":
        return "H3C" if reduced else "H2C=O"
    if subtype == "aldehyde":
        context = contexts[0] if contexts else "R"
        center = "CH2" if reduced else "CH=O"
        return f"{context}{bond}{center}"
    while len(contexts) < 2:
        contexts.append("R")
    left, right = contexts[:2]
    center = "CH" if reduced else "C"
    suffix = "" if reduced else "=O"
    return f"{left}{bond}{center}({right}){suffix}"


def _reduced_amine_fragment(amine: ReactionSiteReference, *, style: str) -> str:
    bond = "–" if style == "unicode" else "-"
    contexts = _ordered_context_labels(amine, style=style)
    remaining_h = max(0, int(amine.details.get("h_count", 0)) - 1)
    if remaining_h >= 2:
        return "NH2"
    if remaining_h == 1:
        return "NH" + (f"{bond}{contexts[0]}" if contexts else "")
    if len(contexts) == 2 and contexts[0] == contexts[1] == "R":
        return "NR2"
    return "N" + "".join(f"({context})" for context in contexts)


def render_reaction_label(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    *,
    style: str = "unicode",
) -> str:
    bond = "–" if style == "unicode" else "-"
    reactant_labels = render_reactant_label(assignment, style=style)
    rule = load_reaction_rendering().get(str(grammar["id"]), {})
    kind = rule.get("product_kind")
    if kind == "join_contexts":
        left, right = assignment[rule["left_role"]], assignment[rule["right_role"]]
        product = render_context_connection(
            left.details["anchor_context"], right.details["anchor_context"], style=style
        )
    elif kind == "nitrogen_substitution":
        product = _nitrogen_product(
            assignment[rule["anchor_role"]].details["anchor_context"],
            assignment[rule["nitrogen_role"]],
            bond,
            style,
        )
    elif kind == "heteroatom_substitution":
        anchor, partner = (
            assignment[rule["anchor_role"]],
            assignment[rule["partner_role"]],
        )
        product = render_heteroatom_product(
            anchor.details["anchor_context"], partner, str(rule["element"]), style=style
        )
    elif kind == "terminal_alkyne":
        anchor = assignment[rule["anchor_role"]]
        product = f"{render_context(anchor.details['anchor_context'], style=style)}{bond}C≡C{bond}R"
    elif kind == "chan_lam":
        anchor, partner = (
            assignment[rule["anchor_role"]],
            assignment[rule["partner_role"]],
        )
        center = partner.details.get("center_token")
        if center in {"N", "N_aromatic"}:
            product = _nitrogen_product(
                anchor.details["anchor_context"], partner, bond, style
            )
        else:
            retained = list(partner.details.get("contexts") or [])
            right = render_context(retained[0], style=style) if retained else "H"
            product = f"{render_context(anchor.details['anchor_context'], style=style)}{bond}{center}{bond}{right}"
    elif kind == "reductive_amination":
        carbonyl = assignment[rule["carbonyl_role"]]
        amine = assignment[rule["amine_role"]]
        reactant_labels = (
            f"{_carbonyl_fragment(carbonyl, reduced=False, style=style)} + "
            f"{amine.chemist_label.replace('–', bond)}"
        )
        product = (
            f"{_carbonyl_fragment(carbonyl, reduced=True, style=style)}{bond}"
            f"{_reduced_amine_fragment(amine, style=style)}"
        )
    elif kind == "amide":
        nitrogen = assignment[rule["nitrogen_role"]]
        product = render_xh(
            "N",
            max(0, int(nitrogen.details["h_count"]) - 1),
            ["C(O)R", *list(nitrogen.details.get("contexts") or [])],
            style=style,
        )
    elif kind == "acyl_heteroatom":
        acyl, partner = assignment[rule["acyl_role"]], assignment[rule["partner_role"]]
        retained = list(partner.details.get("contexts") or [])
        right = render_context(retained[0], style=style) if retained else "H"
        product = f"{render_context(acyl.details['retained_context'], style=style)}{bond}C(O){bond}{rule['element']}{bond}{right}"
    elif kind == "sulfonamide":
        sulfonyl, nitrogen = (
            assignment[rule["sulfonyl_role"]],
            assignment[rule["nitrogen_role"]],
        )
        suffix = "2" if int(nitrogen.details["h_count"]) == 2 else "R"
        product = f"{render_context(sulfonyl.details['retained_context'], style=style)}{bond}S(O)2{bond}NH{suffix}"
    else:
        product = str(grammar.get("generic_label") or grammar["id"])
    return f"{reactant_labels} → {product}"


__all__ = [
    "load_product_context_precedence",
    "load_reaction_rendering",
    "render_context_connection",
    "render_heteroatom_product",
    "render_reactant_label",
    "render_reaction_label",
]
