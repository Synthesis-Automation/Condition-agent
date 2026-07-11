"""Deterministic chemist-facing label rendering."""

from __future__ import annotations

from typing import List


_PREFIX = {"Ar": "Ar", "HeteroAr": "HeteroAr", "Alkenyl": "Alkenyl", "Alkynyl": "Alkynyl", "Alkyl": "R"}


def render_edge(context: str, handle: str) -> str:
    return f"{_PREFIX.get(context, context)}–{handle}"


def render_xh(center: str, h_count: int, contexts: List[str]) -> str:
    if center == "Csp":
        return "R–C≡C–H"
    if center == "O":
        return f"{_PREFIX.get(contexts[0], contexts[0]) if contexts else 'H'}–OH"
    if center == "S":
        return f"{_PREFIX.get(contexts[0], contexts[0]) if contexts else 'H'}–SH"
    if center != "N":
        return f"{center}–H"
    if h_count == 3 and not contexts:
        return "NH3"
    if h_count == 2:
        prefix = _PREFIX.get(contexts[0], contexts[0]) if contexts else "H"
        return f"{prefix}–NH2"
    if contexts == ["Alkyl", "Alkyl"]:
        return "R1R2–NH"
    if contexts == ["Ar", "Ar"]:
        return "Ar1Ar2–NH"
    if sorted(contexts) == ["Alkyl", "Ar"]:
        return "Ar–NH–R"
    if contexts and contexts[0] == "SO2R":
        return "RSO2–NHR" if len(contexts) > 1 else "RSO2–NH2"
    if contexts and contexts[0].startswith("C(O)"):
        return "RC(O)–NHR" if len(contexts) > 1 else "RC(O)–NH2"
    if contexts == ["N"]:
        return "N–NH2" if h_count == 2 else "N–NH–R"
    return "–".join([*(_PREFIX.get(c, c) for c in contexts), "NH"])


__all__ = ["render_edge", "render_xh"]
