"""
LLM-assisted review for uncertain reaction featurization outputs.

This module is intentionally advisory-only: callers decide if and how to apply
LLM suggestions after deterministic taxonomy validation.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional

from chemtools.taxonomy.reaction_catalog import list_reaction_type_ids
from llmtools.clients import LLMClient
from llmtools.prompts import get_template


def _strip_markdown_fences(content: str) -> str:
    """Remove optional markdown fences from model output."""
    text = (content or "").strip()
    if text.startswith("```json") or text.startswith("```JSON"):
        text = text[7:].lstrip()
    elif text.startswith("```"):
        text = text[3:].lstrip()
    if text.endswith("```"):
        text = text[:-3].rstrip()
    return text.strip()


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _format_values(values: Iterable[Any], limit: int = 60) -> str:
    items = [str(v).strip() for v in values if str(v).strip()]
    if len(items) > limit:
        items = items[:limit]
    return json.dumps(items, ensure_ascii=True)


def _format_text(value: Any, max_len: int = 1200) -> str:
    text = str(value or "").strip()
    if not text:
        return "None"
    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


def _format_json(value: Any, max_len: int = 2000) -> str:
    if value is None:
        return "None"
    try:
        text = json.dumps(value, ensure_ascii=True, sort_keys=True)
    except Exception:
        return _format_text(value, max_len=max_len)
    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


@dataclass
class LLMReactionFeaturizationOptions:
    """Runtime configuration for LLM featurization review."""

    enabled: bool = False
    provider: Optional[str] = None
    model: Optional[str] = None
    temperature: float = 0.0
    max_tokens: int = 700
    timeout: int = 60

    def is_ready(self) -> bool:
        return bool(self.enabled and self.provider and self.model)


def build_reaction_featurization_prompt(context: Dict[str, Any]) -> str:
    """Build the review prompt from deterministic featurization context."""
    try:
        candidates = list_reaction_type_ids()
    except Exception:
        candidates = []
    candidate_text = "\n".join(candidates) if candidates else "Unknown"

    template = get_template("reaction_featurization_review")
    return template.format(
        reaction_smiles=_format_text(context.get("reaction_smiles")),
        normalized_reaction=_format_text(context.get("normalized")),
        deterministic_reaction_type=_format_text(
            context.get("deterministic_reaction_type") or "Unknown"
        ),
        deterministic_confidence=str(_safe_float(context.get("deterministic_confidence"), 0.0)),
        mapping_warning=_format_text(context.get("mapping_warning")),
        reaction_key_raw=_format_text(context.get("reaction_key_raw")),
        reaction_key=_format_text(context.get("reaction_key")),
        reacted_motifs=_format_values(context.get("reacted_motifs", [])),
        formed_motifs=_format_values(context.get("formed_motifs", [])),
        spectator_motifs=_format_values(context.get("spectator_motifs", [])),
        reacted_formed_overlap=_format_values(context.get("reacted_formed_overlap", [])),
        product_broad_tags=_format_values(context.get("product_broad_tags", [])),
        product_motifs_reactive=_format_values(context.get("product_motifs_reactive", [])),
        event_kinds=_format_values(context.get("event_kinds", [])),
        reaction_key_quality=_format_json(context.get("reaction_key_quality")),
        stoichiometry_delta=_format_json(context.get("stoichiometry_delta")),
        reacted_motif_counts=_format_json(context.get("reacted_motif_counts")),
        formed_motif_counts=_format_json(context.get("formed_motif_counts")),
        spectator_motif_counts=_format_json(context.get("spectator_motif_counts")),
        reaction_type_candidates=candidate_text,
    )


def review_reaction_featurization(
    context: Dict[str, Any],
    options: LLMReactionFeaturizationOptions,
    client: Optional[LLMClient] = None,
) -> Dict[str, Any]:
    """
    Request a structured LLM review for uncertain reaction featurization output.

    Returns structured diagnostics even when parsing fails.
    """
    if not options.is_ready():
        return {"enabled": False, "status": "disabled"}

    prompt = build_reaction_featurization_prompt(context)

    llm_client = client
    if llm_client is None:
        llm_client = LLMClient(
            provider=options.provider or "openai",
            model=options.model,
            temperature=options.temperature,
            max_tokens=options.max_tokens,
            timeout=options.timeout,
        )

    try:
        response = llm_client.chat(
            prompt=prompt,
            temperature=options.temperature,
            max_tokens=options.max_tokens,
        )
    except Exception as exc:
        return {
            "enabled": True,
            "status": "error",
            "error": str(exc),
            "provider": options.provider,
            "model": options.model,
        }

    raw_content = (response.content or "").strip()
    parsed: Optional[Dict[str, Any]] = None
    parse_error: Optional[str] = None
    if raw_content:
        cleaned = _strip_markdown_fences(raw_content)
        try:
            parsed = json.loads(cleaned)
        except json.JSONDecodeError as exc:
            parse_error = f"{exc.__class__.__name__}: {exc}"

    output: Dict[str, Any] = {
        "enabled": True,
        "status": "ok" if parsed else "parse_error",
        "provider": options.provider,
        "model": response.model,
        "raw_response": raw_content,
        "prompt_tokens": response.prompt_tokens,
        "completion_tokens": response.completion_tokens,
        "total_tokens": response.total_tokens,
        "latency_ms": response.latency_ms,
        "analysis": None,
    }

    if parse_error:
        output["error"] = parse_error

    if isinstance(parsed, dict):
        flags = parsed.get("uncertainty_flags")
        if isinstance(flags, list):
            norm_flags = [str(x) for x in flags if str(x).strip()]
        else:
            norm_flags = []
        output["analysis"] = {
            "suggested_reaction_type": str(
                parsed.get("suggested_reaction_type") or "Unknown"
            ).strip(),
            "confidence": max(0.0, min(1.0, _safe_float(parsed.get("confidence"), 0.0))),
            "rationale": str(parsed.get("rationale") or "").strip(),
            "requires_human_review": bool(parsed.get("requires_human_review", False)),
            "uncertainty_flags": norm_flags,
            "mechanistic_family": str(parsed.get("mechanistic_family") or "").strip(),
            "mechanistic_rationale": str(parsed.get("mechanistic_rationale") or "").strip(),
            "tautomer_or_representation_issue": bool(
                parsed.get("tautomer_or_representation_issue", False)
            ),
            "taxonomy_gap_suspected": bool(parsed.get("taxonomy_gap_suspected", False)),
            "taxonomy_gap_note": str(parsed.get("taxonomy_gap_note") or "").strip(),
            "deterministic_checks_used": [
                str(v).strip()
                for v in (parsed.get("deterministic_checks_used") or [])
                if str(v).strip()
            ]
            if isinstance(parsed.get("deterministic_checks_used"), list)
            else [],
        }

    return output


__all__ = [
    "LLMReactionFeaturizationOptions",
    "build_reaction_featurization_prompt",
    "review_reaction_featurization",
]
