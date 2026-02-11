from __future__ import annotations

from chemtools.featurizers.formatters import reaction as reaction_formatter
from chemtools.featurizers.formatters import detection_validation
from chemtools.featurizers import reaction_assist


def test_get_crk_options_disables_primary_detection() -> None:
    opts = reaction_formatter.get_crk_options()
    assert opts.get("include_reaction_type") is False


def test_primary_detection_is_exposed_and_applied(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        reaction_formatter,
        "_run_primary_reaction_type_detection",
        lambda _: {
            "status": "ok",
            "error": None,
            "summary": {
                "reaction_type": "C_N_Coupling",
                "name": "C_N_Coupling",
                "confidence": 0.91,
                "slot_evidence": {"electrophile": ["Ar-I"], "nucleophile": ["Ar-NH2"]},
            },
            "matches": [
                {
                    "reaction_type": "C_N_Coupling",
                    "name": "C_N_Coupling",
                    "confidence": 0.91,
                    "slot_evidence": {"electrophile": ["Ar-I"], "nucleophile": ["Ar-NH2"]},
                }
            ],
        },
    )
    monkeypatch.setattr(
        detection_validation,
        "validate_detection_with_crk_key",
        lambda **kwargs: {
            "reaction_type": kwargs.get("initial_detection", "Unknown"),
            "confidence": kwargs.get("initial_confidence", 0.0),
            "validation_method": "test_stub",
            "reason": "stubbed",
            "corrected_from": None,
            "evidence": {},
        },
    )

    result = reaction_formatter.featurize_reaction(
        "Ic1ccccc1.N>>Nc1ccccc1",
        options={"include_reaction_type": True},
    )
    assert result.get("reaction_type") == "C_N_Coupling"
    detection = result.get("detection") or {}
    assert (detection.get("primary_detection") or {}).get("status") == "ok"
    assert detection.get("matches")


def test_low_reaction_key_quality_adds_meta_warning(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        reaction_formatter,
        "summarize_reaction_events",
        lambda **_: {
            "events": [],
            "anomalies": ["mapping_unreliable_fallback_used"],
            "reaction_key_quality": {
                "score_0_1": 0.2,
                "level": "low",
                "reasons": ["forced_test_case"],
            },
        },
    )

    result = reaction_formatter.featurize_reaction(
        "Ic1ccccc1.N>>Nc1ccccc1",
        options={"include_reaction_type": False},
    )
    errors = (result.get("meta") or {}).get("errors") or []
    assert "low_reaction_key_quality" in errors
    assert "reaction_key_anomaly:mapping_unreliable_fallback_used" in errors


def test_llm_uncertainty_includes_key_quality_reason() -> None:
    uncertain, reasons = reaction_assist.is_reaction_uncertain_for_llm_assist(
        reaction_type={"reaction_type": "C_N_Coupling", "confidence": 0.92},
        detection_payload={
            "reaction_key_quality": {
                "score_0_1": 0.3,
                "level": "low",
            }
        },
        reaction_key="CRK-v1 | Ar-I -> Ar-NR2",
        confidence_threshold=0.60,
    )
    assert uncertain is True
    assert "low_reaction_key_quality" in reasons
