"""
Utilities for canonical reaction-key splitting and event payload extraction.

The canonical minimal key keeps only:
  |reacted motifs -> formed motifs | spectators: ...

All extra sections (bond/event/reaction-type diagnostics) are collected into a
separate structured payload for optional downstream use.
"""

from __future__ import annotations

import json
import re
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple


_EMPTY_KEY_MARKERS = {
    "",
    "none",
    "[]->[]||[]",
    "crk-v1|[]->[]",
    "|[]->[]",
}


def _split_motif_tokens(value: str) -> List[str]:
    if not value:
        return []
    return [tok.strip() for tok in re.split(r"[|,;/]", str(value)) if tok.strip()]


def _split_event_tokens(value: str) -> List[str]:
    if not value:
        return []
    return [tok.strip() for tok in str(value).split("+") if tok.strip()]


def _split_bond_tokens(value: str) -> List[str]:
    if not value:
        return []
    return [tok.strip() for tok in str(value).split(";") if tok.strip()]


def _dedupe_sorted(values: Iterable[str]) -> List[str]:
    cleaned = {str(value).strip() for value in values if str(value).strip()}
    return sorted(cleaned)


def _format_motif_tokens(values: Iterable[str]) -> str:
    tokens = _dedupe_sorted(values)
    return "|".join(tokens) if tokens else "[]"


def _split_sections(reaction_key: str) -> Tuple[str, List[str]]:
    parts = [part.strip() for part in str(reaction_key or "").split(" | ") if part.strip()]
    if not parts:
        return "", []
    summary = parts[0]
    if summary.startswith("CRK-v1"):
        summary = summary[len("CRK-v1"):].strip()
    if summary.startswith("|"):
        summary = summary[1:].strip()
    return summary, parts[1:]


def _extract_named_sections(parts: Iterable[str]) -> Dict[str, str]:
    named: Dict[str, str] = {}
    for part in parts:
        if ":" not in part:
            continue
        label, payload = part.split(":", 1)
        key = str(label).strip().lower()
        value = str(payload).strip()
        if key and value and key not in named:
            named[key] = value
    return named


def canonicalize_reaction_key_minimal(reaction_key: Any) -> str:
    """
    Convert CRK text into minimal canonical form used by recommendation matching.

    Output shape:
      |<reacted> -> <formed> | spectators: <spectators>
    """
    text = str(reaction_key or "").strip()
    if not text:
        return ""
    compact = text.replace(" ", "").lower()
    if compact in _EMPTY_KEY_MARKERS:
        return ""
    if "->" not in text:
        return text

    summary, rest = _split_sections(text)
    if "->" not in summary:
        return text
    left, right = [chunk.strip() for chunk in summary.split("->", 1)]

    reacted = _split_motif_tokens(left)
    formed = _split_motif_tokens(right)
    base = f"|{_format_motif_tokens(reacted)} -> {_format_motif_tokens(formed)}"

    named = _extract_named_sections(rest)
    spectators = _split_motif_tokens(named.get("spectators", ""))
    if spectators:
        base = f"{base} | spectators: {_format_motif_tokens(spectators)}"
    return base


def build_reaction_events_payload(
    reaction_key: Any,
    reaction_events: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Build a standardized structured payload for non-core reaction key sections.
    """
    text = str(reaction_key or "").strip()
    if not text:
        text = ""
    _, rest = _split_sections(text)
    named = _extract_named_sections(rest)

    payload: Dict[str, Any] = {}

    bond_formed = _split_bond_tokens(named.get("bond_formed", ""))
    bond_formed_labeled = _split_bond_tokens(named.get("bond_formed_labeled", ""))
    bond_broken = _split_bond_tokens(named.get("bond_broken", ""))
    event_signature = str(named.get("events") or "").strip()
    reaction_types = _split_event_tokens(named.get("reaction_types", ""))

    if bond_formed:
        payload["bond_formed"] = _dedupe_sorted(bond_formed)
    if bond_formed_labeled:
        payload["bond_formed_labeled"] = _dedupe_sorted(bond_formed_labeled)
    if bond_broken:
        payload["bond_broken"] = _dedupe_sorted(bond_broken)
    if event_signature:
        payload["event_signature"] = event_signature
    if reaction_types:
        payload["reaction_types"] = _dedupe_sorted(reaction_types)

    if isinstance(reaction_events, Mapping):
        events = reaction_events.get("events")
        if isinstance(events, list) and events:
            kinds = []
            for event in events:
                if not isinstance(event, Mapping):
                    continue
                kind = str(event.get("kind") or "").strip()
                if kind:
                    kinds.append(kind)
            if kinds:
                payload["event_kinds"] = _dedupe_sorted(kinds)

        quality = reaction_events.get("reaction_key_quality")
        if isinstance(quality, Mapping):
            quality_payload: Dict[str, Any] = {}
            level = str(quality.get("level") or "").strip()
            if level:
                quality_payload["level"] = level
            try:
                score = quality.get("score_0_1")
                if score is not None:
                    quality_payload["score_0_1"] = round(float(score), 3)
            except Exception:
                pass
            reasons = quality.get("reasons")
            if isinstance(reasons, list):
                cleaned_reasons = _dedupe_sorted([str(reason) for reason in reasons])
                if cleaned_reasons:
                    quality_payload["reasons"] = cleaned_reasons
            if quality_payload:
                payload["reaction_key_quality"] = quality_payload

        redox = reaction_events.get("redox_assessment")
        if isinstance(redox, Mapping):
            classification = str(redox.get("classification") or "").strip()
            if classification:
                payload["redox_classification"] = classification
                payload["redox_neutral"] = classification == "redox_neutral"
            try:
                confidence = redox.get("confidence")
                if confidence is not None:
                    payload["redox_confidence"] = round(float(confidence), 3)
            except Exception:
                pass

    return payload


def serialize_reaction_events_payload(payload: Mapping[str, Any]) -> str:
    """
    Serialize standardized event payload to a concise single-line format.

    Format:
      sig:<...> | form:<...> | break:<...> | redox:<...> | q:<level>(<score>)
    """
    if not payload:
        return ""
    parts: List[str] = []
    event_sig = str(payload.get("event_signature") or "").strip()
    if event_sig:
        parts.append(f"sig:{event_sig}")

    formed = _dedupe_sorted(payload.get("bond_formed") or [])
    if formed:
        parts.append("form:" + ";".join(formed))

    formed_labeled = _dedupe_sorted(payload.get("bond_formed_labeled") or [])
    if formed_labeled:
        parts.append("form_labeled:" + ";".join(formed_labeled))

    broken = _dedupe_sorted(payload.get("bond_broken") or [])
    if broken:
        parts.append("break:" + ";".join(broken))

    reaction_types = _dedupe_sorted(payload.get("reaction_types") or [])
    if reaction_types:
        parts.append("types:" + "+".join(reaction_types))

    redox = str(payload.get("redox_classification") or "").strip()
    if redox:
        parts.append(f"redox:{redox}")
    redox_conf = payload.get("redox_confidence")
    try:
        if redox_conf is not None:
            parts.append(f"redox_conf:{float(redox_conf):.2f}")
    except Exception:
        pass

    quality = payload.get("reaction_key_quality")
    if isinstance(quality, Mapping):
        q_level = str(quality.get("level") or "").strip()
        q_score = quality.get("score_0_1")
        q_chunk = ""
        if q_level:
            q_chunk = f"q:{q_level}"
        if q_score is not None:
            try:
                q_chunk = f"{q_chunk}({float(q_score):.2f})" if q_chunk else f"q:({float(q_score):.2f})"
            except Exception:
                pass
        if q_chunk:
            parts.append(q_chunk)
        reasons = _dedupe_sorted(quality.get("reasons") or [])
        if reasons:
            parts.append("q_reason:" + ",".join(reasons))

    kinds = _dedupe_sorted(payload.get("event_kinds") or [])
    if kinds:
        parts.append("kinds:" + "+".join(kinds))

    return " | ".join(parts)


def deserialize_reaction_events_text(text: Any) -> Dict[str, Any]:
    """
    Parse event text from either legacy JSON or compact key-value format.
    """
    raw = str(text or "").strip()
    if not raw:
        return {}
    if raw.startswith("{"):
        try:
            parsed = json.loads(raw)
            return parsed if isinstance(parsed, dict) else {}
        except Exception:
            return {}

    payload: Dict[str, Any] = {}
    parts = [part.strip() for part in raw.split(" | ") if part.strip()]
    for part in parts:
        if ":" not in part:
            continue
        key, value = part.split(":", 1)
        label = key.strip().lower()
        data = value.strip()
        if not data:
            continue
        if label == "sig":
            payload["event_signature"] = data
        elif label == "form":
            payload["bond_formed"] = _split_bond_tokens(data)
        elif label == "form_labeled":
            payload["bond_formed_labeled"] = _split_bond_tokens(data)
        elif label == "break":
            payload["bond_broken"] = _split_bond_tokens(data)
        elif label == "types":
            payload["reaction_types"] = _split_event_tokens(data)
        elif label == "redox":
            payload["redox_classification"] = data
            payload["redox_neutral"] = data == "redox_neutral"
        elif label == "redox_conf":
            try:
                payload["redox_confidence"] = float(data)
            except Exception:
                pass
        elif label == "kinds":
            payload["event_kinds"] = _split_event_tokens(data)
        elif label == "q":
            # q:high(0.85) or q:(0.85)
            match = re.match(r"(?P<level>[a-zA-Z_]+)?(?:\((?P<score>[^)]+)\))?$", data)
            q_payload: Dict[str, Any] = {}
            if match:
                level = str(match.group("level") or "").strip()
                score = str(match.group("score") or "").strip()
                if level:
                    q_payload["level"] = level
                if score:
                    try:
                        q_payload["score_0_1"] = float(score)
                    except Exception:
                        pass
            if q_payload:
                payload["reaction_key_quality"] = q_payload
        elif label == "q_reason":
            quality = payload.setdefault("reaction_key_quality", {})
            if isinstance(quality, dict):
                quality["reasons"] = [tok.strip() for tok in data.split(",") if tok.strip()]
    return payload


def normalize_reaction_events_text(value: Any) -> str:
    payload = deserialize_reaction_events_text(value)
    if not payload:
        return ""
    return serialize_reaction_events_payload(payload)


__all__ = [
    "canonicalize_reaction_key_minimal",
    "build_reaction_events_payload",
    "serialize_reaction_events_payload",
    "deserialize_reaction_events_text",
    "normalize_reaction_events_text",
]
