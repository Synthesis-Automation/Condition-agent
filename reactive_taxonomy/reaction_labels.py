"""Taxonomy-driven chemist-friendly reaction rendering."""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional

from .labels import render_context
from .reaction_models import ReactionSiteReference
from .reaction_topology import assignment_component_scope


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


@dataclass(frozen=True)
class _FragmentAliasIndex:
    """Reaction-local aliases for graph-derived generic fragments."""

    aliases_by_fragment: Dict[tuple[Any, ...], str]
    fragments_by_role: Dict[str, tuple[tuple[Any, ...], ...]]
    symbols_by_role: Dict[str, tuple[str, ...]]

    def alias_for(self, role: str, context_index: int) -> Optional[str]:
        fragments = self.fragments_by_role.get(role, ())
        if context_index >= len(fragments):
            return None
        return self.aliases_by_fragment.get(fragments[context_index])

    def aliases_for_symbol(self, role: str, symbol: str) -> tuple[str, ...]:
        fragments = self.fragments_by_role.get(role, ())
        symbols = self.symbols_by_role.get(role, ())
        return tuple(
            self.aliases_by_fragment[fragment]
            for fragment, fragment_symbol in zip(fragments, symbols)
            if fragment_symbol == symbol
        )


def _fragment_alias_index(
    assignment: Dict[str, ReactionSiteReference],
) -> _FragmentAliasIndex:
    """Assign one deterministic display namespace across all selected roles."""
    definition = _load_reaction_rendering_definition()
    indexing = definition.get("fragment_indexing") or {}
    token_symbols = {
        str(token): str(symbol)
        for token, symbol in (indexing.get("context_symbols") or {}).items()
    }
    alias_template = str(indexing.get("alias_template") or "{symbol}{index}")
    ordered_fragments: list[tuple[Any, ...]] = []
    fragment_symbols: Dict[tuple[Any, ...], str] = {}
    fragments_by_role: Dict[str, tuple[tuple[Any, ...], ...]] = {}
    symbols_by_role: Dict[str, tuple[str, ...]] = {}
    for role, site in assignment.items():
        role_fragments: list[tuple[Any, ...]] = []
        role_symbols: list[str] = []
        for context_index, record in enumerate(
            site.details.get("context_records") or ()
        ):
            token = str(record.get("token") or "")
            symbol = token_symbols.get(token)
            atom_indices = tuple(
                int(index) for index in record.get("fragment_atom_indices") or ()
            )
            fragment = (
                (
                    site.component_index,
                    atom_indices,
                    token,
                )
                if atom_indices
                else (
                    site.component_index,
                    role,
                    context_index,
                    token,
                )
            )
            role_fragments.append(fragment)
            role_symbols.append(symbol or "")
            if symbol and fragment not in fragment_symbols:
                fragment_symbols[fragment] = symbol
                ordered_fragments.append(fragment)
        fragments_by_role[role] = tuple(role_fragments)
        symbols_by_role[role] = tuple(role_symbols)
    counts = Counter(fragment_symbols.values())
    next_index: Counter[str] = Counter()
    aliases: Dict[tuple[Any, ...], str] = {}
    for fragment in ordered_fragments:
        symbol = fragment_symbols[fragment]
        if counts[symbol] > 1:
            next_index[symbol] += 1
            aliases[fragment] = alias_template.format(
                symbol=symbol,
                index=next_index[symbol],
            )
        else:
            aliases[fragment] = symbol
    return _FragmentAliasIndex(aliases, fragments_by_role, symbols_by_role)


def _render_site_label(
    role: str,
    site: ReactionSiteReference,
    aliases: _FragmentAliasIndex,
    *,
    style: str,
) -> str:
    """Render one handle using reaction-global rather than site-local indices."""
    if site.details.get("center_family") == "Carbonyl":
        return _carbonyl_fragment(
            site,
            reduced=False,
            style=style,
            role=role,
            aliases=aliases,
        )
    bond = "–" if style == "unicode" else "-"
    label = site.chemist_label.replace("–", bond)
    symbols = sorted(
        set(aliases.symbols_by_role.get(role, ())) - {""}, key=len, reverse=True
    )
    for symbol in symbols:
        replacements = iter(aliases.aliases_for_symbol(role, symbol))
        pattern = re.compile(rf"(?<![A-Za-z]){re.escape(symbol)}(?:\d+)?")
        label = pattern.sub(lambda match: next(replacements, match.group(0)), label)
    return label


def _context_label(
    role: str,
    site: ReactionSiteReference,
    context_index: int,
    token: str,
    aliases: _FragmentAliasIndex,
    *,
    style: str,
) -> str:
    return aliases.alias_for(role, context_index) or render_context(token, style=style)


def _nitrogen_retained_labels(
    role: str,
    nitrogen: ReactionSiteReference,
    aliases: _FragmentAliasIndex,
    *,
    style: str,
) -> list[str]:
    """Render N substituents in attachment-oriented form with global aliases."""
    bond = "–" if style == "unicode" else "-"
    attachment_oriented = {
        "C(O)R": f"C(O){bond}R",
        "C(O)OR": f"C(O){bond}OR",
        "C(O)NR": f"C(O){bond}NR2",
        "SO2R": f"S(O)2{bond}R",
        "P(O)R2": "P(O)R2",
    }
    return [
        aliases.alias_for(role, index)
        or attachment_oriented.get(str(token), render_context(str(token), style=style))
        for index, token in enumerate(nitrogen.details.get("contexts") or ())
    ]


@lru_cache(maxsize=1)
def _load_reaction_rendering_definition() -> Dict[str, Any]:
    with _PATH.open("r", encoding="utf-8") as handle:
        return dict(json.load(handle))


@lru_cache(maxsize=1)
def load_reaction_rendering() -> Dict[str, Dict[str, Any]]:
    return _load_reaction_rendering_definition().get("rules") or {}


@lru_cache(maxsize=1)
def load_product_context_precedence() -> tuple[str, ...]:
    return tuple(
        _load_reaction_rendering_definition().get("product_context_precedence") or ()
    )


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
    anchor_context: str,
    nitrogen: ReactionSiteReference,
    bond: str,
    style: str,
    *,
    anchor_label: Optional[str] = None,
    retained_labels: Optional[list[str]] = None,
) -> str:
    if nitrogen.details.get("center_token") == "N_aromatic":
        return (
            f"{anchor_label or render_context(anchor_context, style=style)}{bond}AromN"
        )
    anchor = anchor_label or render_context(anchor_context, style=style)

    def n_substituent(token: str) -> str:
        attachment_oriented = {
            "C(O)R": f"C(O){bond}R",
            "C(O)OR": f"C(O){bond}OR",
            "C(O)NR": f"C(O){bond}NR2",
            "SO2R": f"S(O)2{bond}R",
            "P(O)R2": "P(O)R2",
        }
        return attachment_oriented.get(token, render_context(token, style=style))

    context_tokens = list(nitrogen.details.get("contexts") or [])
    retained = retained_labels or [n_substituent(token) for token in context_tokens]
    remaining_h = max(0, int(nitrogen.details["h_count"]) - 1)
    if remaining_h >= 2:
        return f"{anchor}{bond}NH2"
    if remaining_h == 1:
        return f"{anchor}{bond}NH" + (f"{bond}{retained[0]}" if retained else "")
    if not retained:
        return f"{anchor}{bond}N"
    if len(retained) == 2 and context_tokens == ["Alkyl", "Alkyl"]:
        return f"{anchor}{bond}N{retained[0]}{retained[1]}"
    return f"{anchor}{bond}N" + "".join(f"({context})" for context in retained)


def render_reactant_label(
    assignment: Dict[str, ReactionSiteReference], *, style: str = "unicode"
) -> str:
    """Render only the observed reactant handles for a candidate assignment."""
    aliases = _fragment_alias_index(assignment)
    scope = assignment_component_scope(assignment)
    separator = " / " if scope == "intramolecular" else " + "
    rendered = separator.join(
        _render_site_label(role, site, aliases, style=style)
        for role, site in assignment.items()
    )
    return f"intramolecular {rendered}" if scope == "intramolecular" else rendered


def render_heteroatom_product(
    anchor_context: str,
    partner: ReactionSiteReference,
    element: str,
    *,
    style: str = "unicode",
    anchor_label: Optional[str] = None,
    retained_label: Optional[str] = None,
) -> str:
    """Render a C–O/C–S product from anchor and retained partner context."""
    bond = "–" if style == "unicode" else "-"
    retained = list(partner.details.get("contexts") or [])
    right = retained_label or (
        render_context(retained[0], style=style) if retained else "H"
    )
    left = anchor_label or render_context(anchor_context, style=style)
    return f"{left}{bond}{element}{bond}{right}"


def _ordered_context_labels(
    site: ReactionSiteReference,
    *,
    style: str,
    role: Optional[str] = None,
    aliases: Optional[_FragmentAliasIndex] = None,
) -> list[str]:
    tokens = list(str(token) for token in site.details.get("contexts") or ())
    entries = sorted(
        enumerate(tokens),
        key=lambda entry: (_PRECEDENCE.get(entry[1], 999), entry[1], entry[0]),
    )
    return [
        (
            _context_label(role, site, index, token, aliases, style=style)
            if aliases is not None and role is not None
            else render_context(token, style=style)
        )
        for index, token in entries
    ]


def _carbonyl_fragment(
    carbonyl: ReactionSiteReference,
    *,
    reduced: bool,
    style: str,
    role: Optional[str] = None,
    aliases: Optional[_FragmentAliasIndex] = None,
) -> str:
    bond = "–" if style == "unicode" else "-"
    subtype = str(carbonyl.details.get("carbonyl_subtype") or "ketone")
    contexts = _ordered_context_labels(
        carbonyl, style=style, role=role, aliases=aliases
    )
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


def _reduced_amine_fragment(
    amine: ReactionSiteReference,
    *,
    style: str,
    role: Optional[str] = None,
    aliases: Optional[_FragmentAliasIndex] = None,
) -> str:
    bond = "–" if style == "unicode" else "-"
    contexts = _ordered_context_labels(amine, style=style, role=role, aliases=aliases)
    remaining_h = max(0, int(amine.details.get("h_count", 0)) - 1)
    if remaining_h >= 2:
        return "NH2"
    if remaining_h == 1:
        return "NH" + (f"{bond}{contexts[0]}" if contexts else "")
    if list(amine.details.get("contexts") or ()) == ["Alkyl", "Alkyl"]:
        return "N" + "".join(contexts)
    return "N" + "".join(f"({context})" for context in contexts)


def render_product_label(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    *,
    style: str = "unicode",
) -> str:
    bond = "–" if style == "unicode" else "-"
    aliases = _fragment_alias_index(assignment)
    rule = load_reaction_rendering().get(str(grammar["id"]), {})
    kind = rule.get("product_kind")
    if kind == "join_contexts":
        left_role, right_role = str(rule["left_role"]), str(rule["right_role"])
        left, right = assignment[left_role], assignment[right_role]
        pairs = [
            (
                str(left.details["anchor_context"]),
                _context_label(
                    left_role,
                    left,
                    0,
                    str(left.details["anchor_context"]),
                    aliases,
                    style=style,
                ),
            ),
            (
                str(right.details["anchor_context"]),
                _context_label(
                    right_role,
                    right,
                    0,
                    str(right.details["anchor_context"]),
                    aliases,
                    style=style,
                ),
            ),
        ]
        precedence = {
            token: index
            for index, token in enumerate(load_product_context_precedence())
        }
        pairs.sort(key=lambda pair: (precedence.get(pair[0], 999), pair[0]))
        product = bond.join(label for _, label in pairs)
    elif kind == "nitrogen_substitution":
        anchor_role, nitrogen_role = (
            str(rule["anchor_role"]),
            str(rule["nitrogen_role"]),
        )
        anchor, nitrogen = assignment[anchor_role], assignment[nitrogen_role]
        retained = _nitrogen_retained_labels(
            nitrogen_role, nitrogen, aliases, style=style
        )
        product = _nitrogen_product(
            str(anchor.details["anchor_context"]),
            nitrogen,
            bond,
            style,
            anchor_label=_context_label(
                anchor_role,
                anchor,
                0,
                str(anchor.details["anchor_context"]),
                aliases,
                style=style,
            ),
            retained_labels=retained,
        )
    elif kind == "heteroatom_substitution":
        anchor_role, partner_role = str(rule["anchor_role"]), str(rule["partner_role"])
        anchor, partner = (
            assignment[anchor_role],
            assignment[partner_role],
        )
        retained = list(partner.details.get("contexts") or [])
        product = render_heteroatom_product(
            str(anchor.details["anchor_context"]),
            partner,
            str(rule["element"]),
            style=style,
            anchor_label=_context_label(
                anchor_role,
                anchor,
                0,
                str(anchor.details["anchor_context"]),
                aliases,
                style=style,
            ),
            retained_label=(
                _context_label(
                    partner_role,
                    partner,
                    0,
                    str(retained[0]),
                    aliases,
                    style=style,
                )
                if retained
                else None
            ),
        )
    elif kind == "terminal_alkyne":
        anchor_role = str(rule["anchor_role"])
        anchor = assignment[anchor_role]
        anchor_label = _context_label(
            anchor_role,
            anchor,
            0,
            str(anchor.details["anchor_context"]),
            aliases,
            style=style,
        )
        product = f"{anchor_label}{bond}C≡C{bond}R"
    elif kind == "heck_alkene":
        anchor = assignment[rule["anchor_role"]]
        alkene = assignment[rule["alkene_role"]]
        right = (
            "CH2" if int(alkene.details.get("substitution_degree", 0)) == 0 else "CHR1"
        )
        product = (
            f"{render_context(anchor.details['anchor_context'], style=style)}"
            f"{bond}CH={right}"
        )
    elif kind == "chan_lam":
        anchor_role, partner_role = str(rule["anchor_role"]), str(rule["partner_role"])
        anchor, partner = (
            assignment[anchor_role],
            assignment[partner_role],
        )
        anchor_label = _context_label(
            anchor_role,
            anchor,
            0,
            str(anchor.details["anchor_context"]),
            aliases,
            style=style,
        )
        center = partner.details.get("center_token")
        if center in {"N", "N_aromatic"}:
            retained = _nitrogen_retained_labels(
                partner_role, partner, aliases, style=style
            )
            product = _nitrogen_product(
                str(anchor.details["anchor_context"]),
                partner,
                bond,
                style,
                anchor_label=anchor_label,
                retained_labels=retained,
            )
        else:
            retained = list(partner.details.get("contexts") or [])
            right = (
                _context_label(
                    partner_role,
                    partner,
                    0,
                    str(retained[0]),
                    aliases,
                    style=style,
                )
                if retained
                else "H"
            )
            product = f"{anchor_label}{bond}{center}{bond}{right}"
    elif kind == "reductive_amination":
        carbonyl_role, amine_role = str(rule["carbonyl_role"]), str(rule["amine_role"])
        carbonyl = assignment[carbonyl_role]
        amine = assignment[amine_role]
        product = (
            f"{_carbonyl_fragment(carbonyl, reduced=True, style=style, role=carbonyl_role, aliases=aliases)}{bond}"
            f"{_reduced_amine_fragment(amine, style=style, role=amine_role, aliases=aliases)}"
        )
    elif kind == "amide":
        acyl_role, nitrogen_role = str(rule["acyl_role"]), str(rule["nitrogen_role"])
        acyl, nitrogen = assignment[acyl_role], assignment[nitrogen_role]
        acyl_context = _context_label(
            acyl_role,
            acyl,
            0,
            str(acyl.details["retained_context"]),
            aliases,
            style=style,
        )
        nitrogen_contexts = _ordered_context_labels(
            nitrogen, style=style, role=nitrogen_role, aliases=aliases
        )
        remaining_h = max(0, int(nitrogen.details["h_count"]) - 1)
        nitrogen_label = "N" + ("H" if remaining_h else "")
        if list(nitrogen.details.get("contexts") or ()) == ["Alkyl", "Alkyl"]:
            nitrogen_label += "".join(nitrogen_contexts)
        elif len(nitrogen_contexts) == 1:
            separator = bond if remaining_h else ""
            nitrogen_label += separator + nitrogen_contexts[0]
        else:
            nitrogen_label += "".join(f"({value})" for value in nitrogen_contexts)
        product = f"{acyl_context}{bond}C(O){bond}{nitrogen_label}"
    elif kind == "acyl_heteroatom":
        acyl_role, partner_role = str(rule["acyl_role"]), str(rule["partner_role"])
        acyl, partner = assignment[acyl_role], assignment[partner_role]
        retained = list(partner.details.get("contexts") or [])
        left = _context_label(
            acyl_role,
            acyl,
            0,
            str(acyl.details["retained_context"]),
            aliases,
            style=style,
        )
        right = (
            _context_label(
                partner_role, partner, 0, str(retained[0]), aliases, style=style
            )
            if retained
            else "H"
        )
        product = f"{left}{bond}C(O){bond}{rule['element']}{bond}{right}"
    elif kind == "aryl_acylation":
        acyl_role, aromatic_role = str(rule["acyl_role"]), str(rule["aromatic_role"])
        acyl = assignment[acyl_role]
        aromatic = assignment[aromatic_role]
        aromatic_label = _context_label(
            aromatic_role,
            aromatic,
            0,
            str(aromatic.details["ring_context"]),
            aliases,
            style=style,
        )
        acyl_label = _context_label(
            acyl_role,
            acyl,
            0,
            str(acyl.details["retained_context"]),
            aliases,
            style=style,
        )
        product = f"{aromatic_label}{bond}C(O){bond}{acyl_label}"
    elif kind == "sulfonamide":
        sulfonyl_role, nitrogen_role = (
            str(rule["sulfonyl_role"]),
            str(rule["nitrogen_role"]),
        )
        sulfonyl, nitrogen = (
            assignment[sulfonyl_role],
            assignment[nitrogen_role],
        )
        left = _context_label(
            sulfonyl_role,
            sulfonyl,
            0,
            str(sulfonyl.details["retained_context"]),
            aliases,
            style=style,
        )
        contexts = _ordered_context_labels(
            nitrogen, style=style, role=nitrogen_role, aliases=aliases
        )
        remaining_h = max(0, int(nitrogen.details["h_count"]) - 1)
        suffix = "H" if remaining_h else ""
        if len(contexts) == 1:
            suffix += (bond if remaining_h else "") + contexts[0]
        elif contexts:
            suffix += "".join(f"({context})" for context in contexts)
        product = f"{left}{bond}S(O)2{bond}N{suffix}"
    else:
        return str(grammar.get("generic_label") or grammar["id"])
    return product


def render_reaction_label(
    grammar: Dict[str, Any],
    assignment: Dict[str, ReactionSiteReference],
    *,
    style: str = "unicode",
) -> str:
    arrow = "→" if style == "unicode" else "->"
    reactant_labels = render_reactant_label(assignment, style=style)
    product = render_product_label(grammar, assignment, style=style)
    return f"{reactant_labels} {arrow} {product}"


__all__ = [
    "load_product_context_precedence",
    "load_reaction_rendering",
    "render_context_connection",
    "render_heteroatom_product",
    "render_product_label",
    "render_reactant_label",
    "render_reaction_label",
]
