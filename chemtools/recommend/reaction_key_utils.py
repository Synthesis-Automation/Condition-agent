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


def _split_csv_tokens(value: str) -> List[str]:
    if not value:
        return []
    return [tok.strip() for tok in re.split(r"[,+]", str(value)) if tok.strip()]


def _canonical_bond_class_token(value: Any) -> str:
    token = str(value or "").strip()
    if "-" not in token:
        return ""
    left, right = [part.strip() for part in token.split("-", 1)]

    def _element(part: str) -> str:
        text = str(part or "").strip()
        if "(" in text:
            text = text.split("(", 1)[0].strip()
        # Keep atom symbol only (e.g., "C(ar)" -> "C").
        match = re.search(r"[A-Z][a-z]?", text)
        return match.group(0) if match else ""

    left_el = _element(left)
    right_el = _element(right)
    if not left_el or not right_el:
        return ""
    return "-".join(sorted((left_el, right_el)))


def _bond_classes_from_pairs(values: Any) -> List[str]:
    if not isinstance(values, (list, tuple, set)):
        return []
    out: List[str] = []
    for pair in values:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            continue
        left = str(pair[0]).strip()
        right = str(pair[1]).strip()
        if not left or not right:
            continue
        out.append("-".join(sorted((left, right))))
    return _dedupe_sorted(out)


def _event_families_from_kinds(values: Any) -> List[str]:
    kinds = {str(v).strip() for v in (values or []) if str(v).strip()}
    if not kinds:
        return []
    families: set[str] = set()
    if "benzyl_o_alkylation_like" in kinds:
        families.add("o_alkylation")
    if "ester_hydrolysis_like" in kinds:
        families.add("hydrolysis")
    if "amidation_like" in kinds:
        families.add("acyl_transfer")
    if "ring_closure_or_annulation" in kinds:
        families.add("annulation")
    bond_formation = {
        "c_n_bond_formation",
        "c_o_bond_formation",
        "c_s_bond_formation",
        "c_c_bond_formation",
    }
    has_bond = bool(kinds & bond_formation)
    has_disp = "leaving_group_displacement" in kinds
    if has_disp and has_bond:
        families.add("substitution")
    elif has_disp:
        families.add("displacement")
    elif has_bond:
        families.add("bond_formation")
    return _dedupe_sorted(families)


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
        payload["formed_bond_classes"] = _dedupe_sorted(
            _canonical_bond_class_token(token) for token in bond_formed
        )
    if bond_formed_labeled:
        payload["bond_formed_labeled"] = _dedupe_sorted(bond_formed_labeled)
    if bond_broken:
        payload["bond_broken"] = _dedupe_sorted(bond_broken)
        payload["broken_bond_classes"] = _dedupe_sorted(
            _canonical_bond_class_token(token) for token in bond_broken
        )
    if event_signature:
        payload["event_signature"] = event_signature
    if reaction_types:
        payload["reaction_types"] = _dedupe_sorted(reaction_types)

    if isinstance(reaction_events, Mapping):
        events = reaction_events.get("events")
        if isinstance(events, list) and events:
            kinds = []
            leaving_groups: List[str] = []
            nucleophile_elements: List[str] = []
            molecularity = ""
            for event in events:
                if not isinstance(event, Mapping):
                    continue
                kind = str(event.get("kind") or "").strip()
                if kind:
                    kinds.append(kind)
                if kind == "intramolecular_likely":
                    molecularity = "intramolecular"
                elif kind == "intermolecular_or_multi_component" and not molecularity:
                    molecularity = "intermolecular_or_multi_component"
                if kind != "leaving_group_displacement":
                    continue
                details = event.get("details")
                if not isinstance(details, Mapping):
                    continue
                lg = str(details.get("leaving_group") or "").strip()
                if lg:
                    leaving_groups.append(lg)
                nuc = str(details.get("nucleophile_element") or "").strip()
                if nuc:
                    nucleophile_elements.append(nuc)
            if kinds:
                payload["event_kinds"] = _dedupe_sorted(kinds)
                if "event_families" not in payload:
                    inferred = _event_families_from_kinds(kinds)
                    if inferred:
                        payload["event_families"] = inferred
            if molecularity:
                payload["molecularity"] = molecularity
            if leaving_groups:
                payload["leaving_groups"] = _dedupe_sorted(leaving_groups)
            if nucleophile_elements:
                payload["nucleophile_elements"] = _dedupe_sorted(nucleophile_elements)

        event_families = reaction_events.get("event_families")
        if isinstance(event_families, list):
            cleaned_families = _dedupe_sorted(event_families)
            if cleaned_families:
                payload["event_families"] = cleaned_families

        anomalies = reaction_events.get("anomalies")
        if isinstance(anomalies, list):
            cleaned_anomalies = _dedupe_sorted(anomalies)
            if cleaned_anomalies:
                payload["anomalies"] = cleaned_anomalies

        ring_change = reaction_events.get("ring_change")
        if isinstance(ring_change, Mapping):
            try:
                ring_delta = int(ring_change.get("delta"))
                payload["ring_delta"] = ring_delta
            except Exception:
                pass

        bond_pairs = reaction_events.get("bond_pairs")
        if isinstance(bond_pairs, Mapping):
            formed_classes = _bond_classes_from_pairs(bond_pairs.get("formed"))
            broken_classes = _bond_classes_from_pairs(bond_pairs.get("broken"))
            if formed_classes:
                payload.setdefault("formed_bond_classes", formed_classes)
            if broken_classes:
                payload.setdefault("broken_bond_classes", broken_classes)

        profile = reaction_events.get("transformation_profile")
        if isinstance(profile, Mapping):
            molecularity = str(profile.get("molecularity") or "").strip()
            if molecularity:
                payload["molecularity"] = molecularity
            for key in ("formed_bond_classes", "broken_bond_classes", "leaving_groups", "nucleophile_elements"):
                values = profile.get(key)
                if isinstance(values, list):
                    cleaned_values = _dedupe_sorted(values)
                    if cleaned_values:
                        payload[key] = cleaned_values
            try:
                ring_delta = profile.get("ring_delta")
                if ring_delta is not None:
                    payload["ring_delta"] = int(ring_delta)
            except Exception:
                pass

        electrophile_profile = reaction_events.get("electrophile_profile")
        if isinstance(electrophile_profile, Mapping):
            hybridization = str(electrophile_profile.get("hybridization_guess") or "").strip()
            if hybridization:
                payload["electrophile_hybridization"] = hybridization
            environment_tags = electrophile_profile.get("environment_tags")
            if isinstance(environment_tags, list):
                cleaned_env = _dedupe_sorted(environment_tags)
                if cleaned_env:
                    payload["electrophile_environment_tags"] = cleaned_env

        nucleophile_profile = reaction_events.get("nucleophile_profile")
        if isinstance(nucleophile_profile, Mapping):
            candidate_classes = nucleophile_profile.get("candidate_classes")
            if isinstance(candidate_classes, list):
                cleaned_classes = _dedupe_sorted(candidate_classes)
                if cleaned_classes:
                    payload["nucleophile_candidate_classes"] = cleaned_classes
            ambident = nucleophile_profile.get("ambident_possible")
            if isinstance(ambident, bool):
                payload["ambident_possible"] = ambident

        mechanism_shortlist = reaction_events.get("mechanism_shortlist")
        if isinstance(mechanism_shortlist, list):
            mechanism_names: List[str] = []
            for item in mechanism_shortlist:
                if isinstance(item, Mapping):
                    name = str(item.get("name") or "").strip()
                else:
                    name = str(item).strip()
                if name:
                    mechanism_names.append(name)
            cleaned_mechanisms = _dedupe_sorted(mechanism_names)
            if cleaned_mechanisms:
                payload["mechanism_shortlist"] = cleaned_mechanisms

        selectivity_risks = reaction_events.get("selectivity_risks")
        if isinstance(selectivity_risks, list):
            cleaned_risks = _dedupe_sorted(selectivity_risks)
            if cleaned_risks:
                payload["selectivity_risks"] = cleaned_risks

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

    families = _dedupe_sorted(payload.get("event_families") or [])
    if families:
        parts.append("fam:" + "+".join(families))

    molecularity = str(payload.get("molecularity") or "").strip()
    if molecularity:
        parts.append(f"mol:{molecularity}")

    formed_classes = _dedupe_sorted(payload.get("formed_bond_classes") or [])
    if formed_classes:
        parts.append("form_cls:" + ";".join(formed_classes))

    broken_classes = _dedupe_sorted(payload.get("broken_bond_classes") or [])
    if broken_classes:
        parts.append("break_cls:" + ";".join(broken_classes))

    leaving_groups = _dedupe_sorted(payload.get("leaving_groups") or [])
    if leaving_groups:
        parts.append("lg:" + "+".join(leaving_groups))

    nucleophile_elements = _dedupe_sorted(payload.get("nucleophile_elements") or [])
    if nucleophile_elements:
        parts.append("nuc:" + "+".join(nucleophile_elements))

    ring_delta = payload.get("ring_delta")
    try:
        if ring_delta is not None:
            parts.append(f"ringd:{int(ring_delta)}")
    except Exception:
        pass

    anomalies = _dedupe_sorted(payload.get("anomalies") or [])
    if anomalies:
        parts.append("anom:" + ",".join(anomalies))

    mechanisms = _dedupe_sorted(payload.get("mechanism_shortlist") or [])
    if mechanisms:
        parts.append("mech:" + "+".join(mechanisms))

    risks = _dedupe_sorted(payload.get("selectivity_risks") or [])
    if risks:
        parts.append("risk:" + ",".join(risks))

    electrophile_hyb = str(payload.get("electrophile_hybridization") or "").strip()
    if electrophile_hyb:
        parts.append(f"ehyb:{electrophile_hyb}")

    electrophile_env = _dedupe_sorted(payload.get("electrophile_environment_tags") or [])
    if electrophile_env:
        parts.append("eenv:" + "+".join(electrophile_env))

    nucleophile_classes = _dedupe_sorted(payload.get("nucleophile_candidate_classes") or [])
    if nucleophile_classes:
        parts.append("nclass:" + "+".join(nucleophile_classes))

    ambident_possible = payload.get("ambident_possible")
    if isinstance(ambident_possible, bool):
        parts.append("amb:1" if ambident_possible else "amb:0")

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
        elif label == "fam":
            payload["event_families"] = _split_event_tokens(data)
        elif label == "mol":
            payload["molecularity"] = data
        elif label == "form_cls":
            payload["formed_bond_classes"] = _split_bond_tokens(data)
        elif label == "break_cls":
            payload["broken_bond_classes"] = _split_bond_tokens(data)
        elif label == "lg":
            payload["leaving_groups"] = _split_event_tokens(data)
        elif label == "nuc":
            payload["nucleophile_elements"] = _split_event_tokens(data)
        elif label == "ringd":
            try:
                payload["ring_delta"] = int(data)
            except Exception:
                pass
        elif label == "anom":
            payload["anomalies"] = [tok.strip() for tok in data.split(",") if tok.strip()]
        elif label == "mech":
            payload["mechanism_shortlist"] = _split_event_tokens(data)
        elif label == "risk":
            payload["selectivity_risks"] = _split_csv_tokens(data)
        elif label == "ehyb":
            payload["electrophile_hybridization"] = data
        elif label == "eenv":
            payload["electrophile_environment_tags"] = _split_event_tokens(data)
        elif label == "nclass":
            payload["nucleophile_candidate_classes"] = _split_event_tokens(data)
        elif label == "amb":
            payload["ambident_possible"] = data in {"1", "true", "yes"}
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
