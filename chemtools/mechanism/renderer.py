"""
Utilities to convert mechanism classification output into agent-friendly narratives.

The renderer keeps Option-1 lightweight: no expensive computations, just structured
text that cites detection evidence, bond changes, and experimental context so LLM
agents can speak confidently.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional


def build_mechanism_narrative(
    mechanism: Dict[str, Any],
    detection: Optional[Dict[str, Any]] = None,
    bond_changes: Optional[Dict[str, Any]] = None,
    context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Convert mechanism outputs into a natural-language summary and highlight list.
    """

    detection = detection or {}
    bond_changes = bond_changes or {}
    context = context or {}

    family = detection.get("family") or "unknown family"
    mechanism_type = mechanism.get("mechanism_type") or "unknown mechanism"
    confidence = mechanism.get("confidence")

    lines: List[str] = []
    intro = f"Reaction classified as {family} following a {mechanism_type} pathway."
    if isinstance(confidence, (int, float)):
        intro += f" Confidence ~{confidence:.2f}."
    lines.append(intro)

    if context:
        ctx_tokens = _format_context(context)
        if ctx_tokens:
            lines.append(f"Conditions: {', '.join(ctx_tokens)}.")

    bond_bits = _summarize_bonds_for_narrative(bond_changes)
    if bond_bits:
        lines.append(bond_bits)

    highlights = _build_highlights(mechanism.get("steps") or [])
    warnings = _collect_renderer_warnings(detection, bond_changes)

    if mechanism.get("warnings"):
        warnings.extend([w for w in mechanism["warnings"] if w and w not in warnings])

    return {
        "narrative": " ".join(lines).strip(),
        "highlights": highlights,
        "warnings": warnings,
    }


def _format_context(context: Dict[str, Any]) -> List[str]:
    tokens: List[str] = []
    solvent = context.get("solvent")
    if solvent:
        tokens.append(f"solvent {solvent}")
    temperature = context.get("temperature_c")
    if isinstance(temperature, (int, float)):
        tokens.append(f"{temperature} C")
    catalyst = context.get("catalyst")
    if catalyst:
        tokens.append(f"catalyst {catalyst}")
    atmosphere = context.get("atmosphere")
    if atmosphere:
        tokens.append(f"{atmosphere} atmosphere")
    return tokens


def _summarize_bonds_for_narrative(bond_changes: Dict[str, Any]) -> Optional[str]:
    broken = bond_changes.get("broken_bonds") or []
    formed = bond_changes.get("formed_bonds") or []
    if not broken and not formed:
        return None
    return f"Bond evidence: {len(broken)} broken / {len(formed)} formed (RXNMapper/MCS)."


def _build_highlights(steps: List[Dict[str, Any]]) -> List[str]:
    highlights: List[str] = []
    for idx, step in enumerate(steps, start=1):
        title = step.get("title") or f"Step {idx}"
        desc = step.get("description") or ""
        highlights.append(f"{idx}. {title}: {desc}".strip())
        if len(highlights) >= 5:
            break
    return highlights


def _collect_renderer_warnings(
    detection: Dict[str, Any],
    bond_changes: Dict[str, Any],
) -> List[str]:
    warnings: List[str] = []
    status = detection.get("status")
    if status == "conflict":
        warnings.append("Rule/ML detection disagreement - review classification.")
    if bond_changes and not (bond_changes.get("broken_bonds") or bond_changes.get("formed_bonds")):
        warnings.append("Bond analysis returned empty lists.")
    return warnings


__all__ = ["build_mechanism_narrative"]
