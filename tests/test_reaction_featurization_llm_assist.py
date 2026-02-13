from __future__ import annotations

import pytest

from chemtools.featurizers import reaction_assist
from chemtools.featurizers import unified as unified_featurizer


def _llm_assist_options() -> dict:
    return {
        "llm_assist": {
            "enabled": True,
            "provider": "openai",
            "model": "gpt-test",
            "only_on_uncertain": True,
        }
    }


def test_llm_assist_applies_taxonomy_validated_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        reaction_assist,
        "is_reaction_uncertain_for_llm_assist",
        lambda **_: (True, ["unknown_reaction_type"]),
    )
    monkeypatch.setattr(
        reaction_assist,
        "_run_llm_reaction_assist",
        lambda context, llm_assist: {
            "status": "ok",
            "provider": llm_assist.get("provider"),
            "model": llm_assist.get("model"),
            "total_tokens": 42,
            "latency_ms": 8.0,
            "analysis": {
                "suggested_reaction_type": "Amide_formation",
                "confidence": 0.9,
                "rationale": "motif pattern supports amide formation",
                "requires_human_review": False,
                "uncertainty_flags": [],
                "mechanistic_family": "acyl_transfer",
                "mechanistic_rationale": "C-N formation adjacent to carbonyl",
                "tautomer_or_representation_issue": False,
                "taxonomy_gap_suspected": False,
                "deterministic_checks_used": ["reaction_key", "event_kinds"],
            },
        },
    )
    monkeypatch.setattr(
        reaction_assist,
        "_validate_llm_suggested_type_with_crk",
        lambda suggested_reaction_type, reaction_key_raw: (True, "validated_with_crk"),
    )

    result = unified_featurizer.featurize_reaction(
        "CCO.CN>>CCN",
        options=_llm_assist_options(),
    )

    assert result.get("reaction_type") == "Amide_formation"
    meta = result.get("meta", {}).get("llm_assist", {})
    assert meta.get("decision") == "applied"
    assert meta.get("status") == "ok"
    assert meta.get("mechanistic_family") == "acyl_transfer"
    assert meta.get("deterministic_checks_used") == ["reaction_key", "event_kinds"]


def test_llm_assist_rejects_override_on_validation_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        reaction_assist,
        "is_reaction_uncertain_for_llm_assist",
        lambda **_: (True, ["unknown_reaction_type"]),
    )
    monkeypatch.setattr(
        reaction_assist,
        "_run_llm_reaction_assist",
        lambda context, llm_assist: {
            "status": "ok",
            "provider": llm_assist.get("provider"),
            "model": llm_assist.get("model"),
            "analysis": {
                "suggested_reaction_type": "Amide_formation",
                "confidence": 0.95,
                "rationale": "suggested correction",
                "requires_human_review": False,
                "uncertainty_flags": [],
            },
        },
    )
    monkeypatch.setattr(
        reaction_assist,
        "_validate_llm_suggested_type_with_crk",
        lambda suggested_reaction_type, reaction_key_raw: (
            False,
            "validation_mismatch:C_N_Coupling",
        ),
    )

    result = unified_featurizer.featurize_reaction(
        "CCO.CN>>CCN",
        options=_llm_assist_options(),
    )

    assert result.get("reaction_type") != "Amide_formation"
    meta = result.get("meta", {}).get("llm_assist", {})
    assert meta.get("decision") == "rejected_validation_mismatch"


def test_llm_assist_skips_when_not_uncertain(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    called = {"count": 0}

    monkeypatch.setattr(
        reaction_assist,
        "is_reaction_uncertain_for_llm_assist",
        lambda **_: (False, []),
    )

    def _never_called(context, llm_assist):
        called["count"] += 1
        return {"status": "ok", "analysis": {}}

    monkeypatch.setattr(reaction_assist, "_run_llm_reaction_assist", _never_called)

    result = unified_featurizer.featurize_reaction(
        "CCO.CN>>CCN",
        options=_llm_assist_options(),
    )

    meta = result.get("meta", {}).get("llm_assist", {})
    assert meta.get("decision") == "deterministic_kept"
    assert meta.get("used") is False
    assert called["count"] == 0


def test_llm_assist_review_context_includes_richer_deterministic_signals(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict = {}

    monkeypatch.setattr(
        reaction_assist,
        "is_reaction_uncertain_for_llm_assist",
        lambda **_: (True, ["low_confidence"]),
    )

    def _capture_context(context, llm_assist):
        captured.update(context)
        return {
            "status": "ok",
            "provider": llm_assist.get("provider"),
            "model": llm_assist.get("model"),
            "analysis": {
                "suggested_reaction_type": "Unknown",
                "confidence": 0.2,
                "rationale": "no override",
                "requires_human_review": True,
                "uncertainty_flags": ["taxonomy_boundary_case"],
            },
        }

    monkeypatch.setattr(reaction_assist, "_run_llm_reaction_assist", _capture_context)

    _ = unified_featurizer.featurize_reaction(
        "Clc1ncc(-c2ccccc2)cn1>>NN=c1ncc(-c2ccccc2)c[nH]1",
        options=_llm_assist_options(),
    )

    assert isinstance(captured.get("stoichiometry_delta"), dict)
    assert "element_delta" in captured.get("stoichiometry_delta", {})
    assert isinstance(captured.get("event_kinds"), list)
    assert "reaction_key_quality" in captured
