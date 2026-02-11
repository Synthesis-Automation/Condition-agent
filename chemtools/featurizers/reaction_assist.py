"""
Optional LLM post-processing for reaction featurization outputs.

This module is intentionally separate from deterministic reaction formatting.
It consumes the deterministic payload and may apply a taxonomy-validated LLM
override only when configured.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple

from .formatters.molecule import to_bool

_LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD = 0.60


def _to_float_or_default(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _normalize_llm_assist_options(options: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    cfg = (options or {}).get("llm_assist") if isinstance(options, dict) else None
    defaults = {
        "enabled": False,
        "provider": None,
        "model": None,
        "temperature": 0.0,
        "max_tokens": 700,
        "timeout": 60,
        "only_on_uncertain": True,
        "confidence_threshold": _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
        "require_crk_validation": True,
    }
    if isinstance(cfg, bool):
        defaults["enabled"] = cfg
        return defaults
    if not isinstance(cfg, dict):
        return defaults
    defaults["enabled"] = to_bool(cfg.get("enabled"), default=True)
    defaults["provider"] = (
        str(cfg.get("provider")).strip() if cfg.get("provider") is not None else None
    )
    defaults["model"] = (
        str(cfg.get("model")).strip() if cfg.get("model") is not None else None
    )
    defaults["temperature"] = _to_float_or_default(cfg.get("temperature"), 0.0)
    defaults["max_tokens"] = int(cfg.get("max_tokens") or 700)
    defaults["timeout"] = int(cfg.get("timeout") or 60)
    defaults["only_on_uncertain"] = to_bool(cfg.get("only_on_uncertain"), default=True)
    defaults["confidence_threshold"] = _to_float_or_default(
        cfg.get("confidence_threshold"),
        _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
    )
    defaults["require_crk_validation"] = to_bool(
        cfg.get("require_crk_validation"),
        default=True,
    )
    return defaults


def _get_reaction_type_id_and_confidence(payload: Dict[str, Any]) -> Tuple[Optional[str], float]:
    reaction_type = payload.get("reaction_type")
    if isinstance(reaction_type, dict):
        rid = reaction_type.get("reaction_type") or reaction_type.get("name")
        conf = _to_float_or_default(
            reaction_type.get("confidence"),
            _to_float_or_default(payload.get("confidence"), 0.0),
        )
    else:
        rid = reaction_type
        conf = _to_float_or_default(payload.get("confidence"), 0.0)

    rid_text = str(rid or "").strip()
    if not rid_text or rid_text.lower() == "unknown":
        return None, conf
    return rid_text, conf


def _set_reaction_type_payload(
    payload: Dict[str, Any],
    reaction_type_id: str,
    confidence: float,
) -> None:
    payload["reaction_type"] = reaction_type_id
    payload["confidence"] = round(
        max(0.0, min(1.0, _to_float_or_default(confidence, 0.0))),
        3,
    )


def is_reaction_uncertain_for_llm_assist(
    reaction_type: Any,
    detection_payload: Dict[str, Any],
    reaction_key: Optional[str],
    confidence_threshold: float,
) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    if isinstance(reaction_type, dict):
        rid = reaction_type.get("reaction_type")
        confidence = _to_float_or_default(reaction_type.get("confidence"), 0.0)
    else:
        rid = reaction_type
        confidence = 0.0
    rid_text = str(rid or "").strip()
    if not rid_text or rid_text.lower() == "unknown":
        reasons.append("unknown_reaction_type")
    if confidence < max(0.0, confidence_threshold):
        reasons.append("low_confidence")
    if isinstance(detection_payload.get("mapping_warning"), dict):
        reasons.append("mapping_warning")
    if not str(reaction_key or "").strip():
        reasons.append("missing_reaction_key")
    key_quality = detection_payload.get("reaction_key_quality")
    if isinstance(key_quality, dict):
        level = str(key_quality.get("level") or "").strip().lower()
        score = _to_float_or_default(key_quality.get("score_0_1"), 1.0)
        if level == "low" or score < 0.45:
            reasons.append("low_reaction_key_quality")
        elif level == "medium":
            reasons.append("medium_reaction_key_quality")
    validation_payload = detection_payload.get("validation")
    if isinstance(validation_payload, dict):
        validated = str(validation_payload.get("validated_detection") or "").strip()
        if not validated or validated.lower() == "unknown":
            reasons.append("unknown_validated_detection")
    return bool(reasons), reasons


@lru_cache(maxsize=1)
def _load_reaction_catalog_data() -> Tuple[Dict[str, Any], Dict[str, str]]:
    try:
        from chemtools.taxonomy.reaction_catalog import load_reaction_catalog
    except Exception:
        return {}, {}
    try:
        definitions, alias_map = load_reaction_catalog()
    except Exception:
        return {}, {}
    return definitions, alias_map


def _resolve_reaction_type_id(reaction_type: Optional[str]) -> Optional[str]:
    if not reaction_type:
        return None
    definitions, alias_map = _load_reaction_catalog_data()
    if not definitions:
        return None
    label = str(reaction_type).strip()
    if not label:
        return None
    if label in definitions:
        return label
    resolved = alias_map.get(label.lower())
    if resolved and resolved in definitions:
        return resolved
    return None


def _run_llm_reaction_assist(
    context: Dict[str, Any],
    llm_assist: Dict[str, Any],
) -> Dict[str, Any]:
    try:
        from llmtools.reaction_featurization_review import (
            LLMReactionFeaturizationOptions,
            review_reaction_featurization,
        )
    except Exception as exc:
        return {
            "enabled": True,
            "status": "error",
            "error": f"llmtools unavailable: {exc}",
        }

    opts = LLMReactionFeaturizationOptions(
        enabled=True,
        provider=llm_assist.get("provider"),
        model=llm_assist.get("model"),
        temperature=_to_float_or_default(llm_assist.get("temperature"), 0.0),
        max_tokens=int(llm_assist.get("max_tokens") or 700),
        timeout=int(llm_assist.get("timeout") or 60),
    )
    return review_reaction_featurization(context, opts)


def _validate_llm_suggested_type_with_crk(
    suggested_reaction_type: str,
    reaction_key_raw: Optional[str],
) -> Tuple[bool, str]:
    if not str(reaction_key_raw or "").strip():
        return True, "no_reaction_key_for_validation"
    try:
        from .formatters.detection_validation import validate_detection_with_crk_key
    except Exception as exc:
        return False, f"validation_import_error: {exc}"
    try:
        validated = validate_detection_with_crk_key(
            initial_detection=suggested_reaction_type,
            initial_confidence=0.9,
            reaction_key=str(reaction_key_raw),
            include_evidence=False,
        )
    except Exception as exc:
        return False, f"validation_error: {exc}"
    validated_type = _resolve_reaction_type_id(validated.get("reaction_type"))
    if validated_type == suggested_reaction_type:
        return True, "validated_with_crk"
    return False, f"validation_mismatch:{validated_type or 'Unknown'}"


def apply_llm_reaction_assist(
    payload: Dict[str, Any],
    *,
    reaction_smiles: str,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    llm_assist = _normalize_llm_assist_options(options)
    if not llm_assist.get("enabled"):
        return payload

    out = dict(payload)
    detection_payload = out.get("detection")
    if not isinstance(detection_payload, dict):
        detection_payload = {}
    else:
        detection_payload = dict(detection_payload)
    aggregates = out.get("aggregates")
    if not isinstance(aggregates, dict):
        aggregates = {}

    llm_assist_meta: Dict[str, Any] = {
        "enabled": True,
        "used": False,
        "status": "skipped",
        "decision": "none",
    }
    if not llm_assist.get("provider") or not llm_assist.get("model"):
        llm_assist_meta["status"] = "config_error"
        llm_assist_meta["decision"] = "missing_provider_or_model"
    else:
        current_rt_id, current_conf = _get_reaction_type_id_and_confidence(out)
        uncertain, uncertainty_reasons = is_reaction_uncertain_for_llm_assist(
            reaction_type={
                "reaction_type": current_rt_id or "Unknown",
                "confidence": current_conf,
            },
            detection_payload=detection_payload,
            reaction_key=out.get("reaction_key"),
            confidence_threshold=_to_float_or_default(
                llm_assist.get("confidence_threshold"),
                _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
            ),
        )
        llm_assist_meta["uncertainty_reasons"] = uncertainty_reasons
        only_on_uncertain = bool(llm_assist.get("only_on_uncertain", True))
        if only_on_uncertain and not uncertain:
            llm_assist_meta["status"] = "skipped_not_uncertain"
            llm_assist_meta["decision"] = "deterministic_kept"
        else:
            llm_assist_meta["used"] = True
            reaction_key_raw = (
                (detection_payload.get("validation") or {}).get("reaction_key_raw")
                if isinstance(detection_payload.get("validation"), dict)
                else None
            )
            review_context = {
                "reaction_smiles": reaction_smiles,
                "normalized": out.get("normalized") or "",
                "deterministic_reaction_type": current_rt_id or "Unknown",
                "deterministic_confidence": current_conf,
                "reaction_key_raw": reaction_key_raw or "",
                "reaction_key": out.get("reaction_key") or "",
                "mapping_warning": detection_payload.get("mapping_warning"),
                "reacted_motifs": aggregates.get("reacted_motifs") or [],
                "formed_motifs": aggregates.get("formed_motifs_all")
                or aggregates.get("formed_motifs")
                or [],
                "spectator_motifs": aggregates.get("spectator_motifs") or [],
                "product_broad_tags": out.get("product_broad_tags") or [],
                "product_motifs_reactive": out.get("product_motifs_reactive") or [],
            }
            review = _run_llm_reaction_assist(review_context, llm_assist)
            llm_assist_meta["status"] = str(review.get("status") or "unknown")
            llm_assist_meta["provider"] = review.get("provider") or llm_assist.get(
                "provider"
            )
            llm_assist_meta["model"] = review.get("model") or llm_assist.get("model")
            if review.get("total_tokens") is not None:
                llm_assist_meta["total_tokens"] = review.get("total_tokens")
            if review.get("latency_ms") is not None:
                llm_assist_meta["latency_ms"] = review.get("latency_ms")
            if review.get("error"):
                llm_assist_meta["error"] = review.get("error")

            analysis = review.get("analysis")
            if review.get("status") == "ok" and isinstance(analysis, dict):
                suggested_label = str(analysis.get("suggested_reaction_type") or "").strip()
                llm_assist_meta["suggested_reaction_type"] = suggested_label or "Unknown"
                llm_assist_meta["suggested_confidence"] = _to_float_or_default(
                    analysis.get("confidence"),
                    0.0,
                )
                llm_assist_meta["requires_human_review"] = bool(
                    analysis.get("requires_human_review", False)
                )
                if analysis.get("uncertainty_flags"):
                    llm_assist_meta["uncertainty_flags"] = list(
                        analysis.get("uncertainty_flags") or []
                    )

                if not suggested_label or suggested_label.lower() == "unknown":
                    llm_assist_meta["decision"] = "no_valid_suggestion"
                else:
                    suggested_rt_id = _resolve_reaction_type_id(suggested_label)
                    if not suggested_rt_id:
                        llm_assist_meta["decision"] = "invalid_taxonomy_label"
                    elif suggested_rt_id == current_rt_id:
                        llm_assist_meta["decision"] = "no_change"
                    else:
                        validation_ok = True
                        validation_reason = "validation_bypassed"
                        if bool(llm_assist.get("require_crk_validation", True)):
                            validation_ok, validation_reason = (
                                _validate_llm_suggested_type_with_crk(
                                    suggested_rt_id,
                                    reaction_key_raw or out.get("reaction_key"),
                                )
                            )
                        if validation_ok:
                            llm_conf = _to_float_or_default(
                                analysis.get("confidence"),
                                current_conf,
                            )
                            merged_conf = max(current_conf, llm_conf)
                            _set_reaction_type_payload(out, suggested_rt_id, merged_conf)
                            llm_assist_meta["decision"] = "applied"
                            llm_assist_meta["validation"] = validation_reason
                        else:
                            llm_assist_meta["decision"] = "rejected_validation_mismatch"
                            llm_assist_meta["validation"] = validation_reason
            elif llm_assist_meta.get("decision") == "none":
                llm_assist_meta["decision"] = "llm_review_failed"

    detection_payload["llm_assist"] = {
        "decision": llm_assist_meta.get("decision"),
        "status": llm_assist_meta.get("status"),
        "used": llm_assist_meta.get("used", False),
    }
    out["detection"] = detection_payload

    meta_payload = out.get("meta")
    if not isinstance(meta_payload, dict):
        meta_payload = {}
    else:
        meta_payload = dict(meta_payload)
    meta_payload["llm_assist"] = llm_assist_meta
    out["meta"] = meta_payload
    return out


__all__ = [
    "apply_llm_reaction_assist",
    "is_reaction_uncertain_for_llm_assist",
]
