from __future__ import annotations

from scripts import triage_no_motif_evidence as triage


def test_classify_missing_side_bucket() -> None:
    bucket, detail = triage.classify_no_motif_evidence_case(
        "CCO>>",
        {
            "reaction_type": "Unresolved:NoMotifEvidence",
            "reaction_key": "",
            "reaction_events": {
                "reaction_key_quality": {
                    "level": "low",
                    "score_0_1": 0.2,
                    "reasons": [
                        "missing_formed_bond_and_product_motif_evidence",
                        "missing_bond_key",
                    ],
                }
            },
            "aggregates": {"reacted_motifs": [], "formed_motifs": []},
        },
    )
    assert bucket == "missing_side_or_malformed_reaction_smiles"
    assert detail["product_count"] == 0


def test_classify_organometallic_bucket() -> None:
    bucket, detail = triage.classify_no_motif_evidence_case(
        "CC[P](CC)(CC)[Pt+2]([Br-])>>CC[P](CC)(CC)[Pt+2]([Br-])",
        {
            "reaction_type": "Unresolved:NoMotifEvidence",
            "reaction_key": "",
            "reaction_events": {
                "reaction_key_quality": {
                    "level": "low",
                    "score_0_1": 0.2,
                    "reasons": [
                        "missing_formed_bond_and_product_motif_evidence",
                        "missing_bond_key",
                    ],
                }
            },
            "aggregates": {"reacted_motifs": [], "formed_motifs": []},
        },
    )
    assert bucket == "organometallic_or_coordination_complex"
    assert "Pt" in (detail.get("metal_tokens") or [])


def test_no_motif_case_gate_requires_empty_motif_sets() -> None:
    assert not triage._is_no_motif_evidence_case(  # noqa: SLF001
        {
            "reaction_type": "Unknown",
            "aggregates": {"reacted_motifs": ["Ar-Br"], "formed_motifs": []},
            "reaction_events": {
                "reaction_key_quality": {
                    "reasons": [
                        "missing_formed_bond_and_product_motif_evidence",
                        "missing_bond_key",
                    ]
                }
            },
        }
    )
