from __future__ import annotations

from chemtools.recommend.reaction_key_utils import (
    build_reaction_events_payload,
    canonicalize_reaction_key_minimal,
    serialize_reaction_events_payload,
)


def test_canonicalize_reaction_key_minimal_keeps_only_summary_and_spectators() -> None:
    raw = (
        "CRK-v1 |Ar-Br|HeteroAr-NH2 -> HeteroAr-NHR "
        "| bond_formed: C(ar)-N "
        "| bond_formed_labeled: C(ar)-N "
        "| bond_broken: Br-C(ar) "
        "| spectators: Alkyl-AromN|Ar-COR|HeteroAr-H "
        "| events: LGDisp+C-N"
    )
    out = canonicalize_reaction_key_minimal(raw)
    assert out == (
        "|Ar-Br|HeteroAr-NH2 -> HeteroAr-NHR "
        "| spectators: Alkyl-AromN|Ar-COR|HeteroAr-H"
    )
    assert "bond_formed" not in out
    assert "events:" not in out


def test_build_reaction_events_payload_includes_redox_and_event_details() -> None:
    raw = (
        "CRK-v1 |Ar-Br|HeteroAr-NH2 -> HeteroAr-NHR "
        "| bond_formed: C(ar)-N "
        "| bond_broken: Br-C(ar) "
        "| events: LGDisp+C-N"
    )
    reaction_events = {
        "events": [
            {"kind": "c_n_bond_formation", "confidence": 0.9},
            {
                "kind": "leaving_group_displacement",
                "confidence": 0.92,
                "details": {"leaving_group": "Br", "nucleophile_element": "N"},
            },
        ],
        "event_families": ["substitution"],
        "ring_change": {"delta": 0},
        "transformation_profile": {
            "molecularity": "intermolecular_or_multi_component",
            "formed_bond_classes": ["C-N"],
            "broken_bond_classes": ["Br-C"],
            "leaving_groups": ["Br"],
            "nucleophile_elements": ["N"],
            "ring_delta": 0,
        },
        "reaction_key_quality": {
            "score_0_1": 0.85,
            "level": "high",
            "reasons": ["multi_transform_record_possible"],
        },
        "redox_assessment": {
            "classification": "redox_neutral",
            "confidence": 0.65,
        },
    }
    payload = build_reaction_events_payload(raw, reaction_events)
    assert payload["bond_formed"] == ["C(ar)-N"]
    assert payload["bond_broken"] == ["Br-C(ar)"]
    assert payload["event_signature"] == "LGDisp+C-N"
    assert payload["redox_neutral"] is True
    assert payload["redox_classification"] == "redox_neutral"
    assert payload["event_families"] == ["substitution"]
    assert payload["molecularity"] == "intermolecular_or_multi_component"
    assert payload["formed_bond_classes"] == ["C-N"]
    assert payload["broken_bond_classes"] == ["Br-C"]
    assert payload["leaving_groups"] == ["Br"]
    assert payload["nucleophile_elements"] == ["N"]
    assert payload["ring_delta"] == 0
    assert "reaction_key_quality" in payload

    text = serialize_reaction_events_payload(payload)
    assert "sig:LGDisp+C-N" in text
    assert "form:C(ar)-N" in text
    assert "break:Br-C(ar)" in text
    assert "redox:redox_neutral" in text
    assert "fam:substitution" in text
    assert "mol:intermolecular_or_multi_component" in text
    assert "form_cls:C-N" in text
    assert "lg:Br" in text
    assert "nuc:N" in text
