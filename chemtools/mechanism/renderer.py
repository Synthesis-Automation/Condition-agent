"""
Utilities to convert mechanism classification output into agent-friendly narratives.

The renderer keeps Option-1 lightweight: no expensive computations, just structured
text that cites detection evidence, bond changes, and experimental context so LLM
agents can speak confidently.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple


def build_mechanism_narrative(
    mechanism: Dict[str, Any],
    detection: Optional[Dict[str, Any]] = None,
    bond_changes: Optional[Dict[str, Any]] = None,
    context: Optional[Dict[str, Any]] = None,
    electron_flow: Optional[Dict[str, Any]] = None,
    intermediates: Optional[List[Dict[str, Any]]] = None,
    precedents: Optional[List[Dict[str, Any]]] = None,
    electron_balance: Optional[Dict[str, Any]] = None,
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

    if electron_flow:
        flow_line, flow_warnings = _summarize_electron_flow(electron_flow)
        if flow_line:
            lines.append(flow_line)
        flow_notes = electron_flow.get("notes") or []
        for note in flow_notes:
            if "No electron flow" in note:
                flow_warnings.append(note)
            else:
                lines.append(note)
    else:
        flow_warnings = []

    extra_highlights: List[str] = []
    if intermediates:
        names = ", ".join(entry.get("name", "intermediate") for entry in intermediates[:3])
        if names:
            lines.append(f"Representative intermediates: {names}.")
            extra_highlights.append(f"Intermediates: {names}")

    precedent_line = _summarize_precedents(precedents or [])
    if precedent_line:
        lines.append(precedent_line)
        extra_highlights.append(precedent_line)

    if electron_balance:
        balance_line = _summarize_electron_balance(electron_balance)
        if balance_line:
            lines.append(balance_line)

    highlights = _build_highlights(mechanism.get("steps") or [])
    highlights.extend(extra_highlights)
    warnings = _collect_renderer_warnings(detection, bond_changes)
    warnings.extend(flow_warnings)

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


def _summarize_electron_flow(electron_flow: Dict[str, Any]) -> Tuple[Optional[str], List[str]]:
    arrows = electron_flow.get("arrows") or []
    if not arrows:
        return None, []
    first = arrows[0]
    description = first.get("description") or ""
    return f"Electron flow: {description}", []


def _summarize_precedents(precedents: List[Dict[str, Any]]) -> Optional[str]:
    if not precedents:
        return None
    snippets: List[str] = []
    for entry in precedents[:2]:
        parts: List[str] = []
        sim = entry.get("similarity")
        if isinstance(sim, (int, float)):
            parts.append(f"{sim:.2f} sim")
        yld = entry.get("yield")
        if yld not in (None, ""):
            parts.append(f"{yld}% yield")
        catalyst = (entry.get("conditions") or {}).get("catalyst")
        if catalyst:
            parts.append(str(catalyst))
        if parts:
            snippets.append("; ".join(parts))
    if not snippets:
        return None
    return "Closest precedents: " + " | ".join(snippets)


def _summarize_electron_balance(balance: Dict[str, Any]) -> Optional[str]:
    atoms = balance.get("atom_balance") or []
    if not atoms:
        return None
    preview = ", ".join(
        f"{entry['atom']}: {entry['delta_lone_pairs']:+}"
        for entry in atoms[:4]
    )
    return f"Electron balance Δ(lone pairs) [reactant→product]: {preview}."


__all__ = ["build_mechanism_narrative"]
