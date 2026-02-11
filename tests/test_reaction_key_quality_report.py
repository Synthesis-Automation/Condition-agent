from __future__ import annotations

from scripts import reaction_key_quality_report as report_mod


def test_build_report_aggregates_quality_metrics(monkeypatch) -> None:
    payloads = {
        "rxn_ok": {
            "reaction_type": "C_N_Coupling",
            "reaction_key": "CRK-v1 | Ar-I -> Ar-NR2",
            "reaction_events": {
                "reaction_key_quality": {"level": "high", "score_0_1": 0.9},
                "anomalies": [],
            },
            "aggregates": {"reacted_motifs": ["Ar-I"], "formed_motifs": ["Ar-NR2"]},
        },
        "rxn_unknown_low": {
            "reaction_type": "Unknown",
            "reaction_key": "",
            "reaction_events": {
                "reaction_key_quality": {"level": "low", "score_0_1": 0.2},
                "anomalies": ["mapping_unreliable_fallback_used"],
            },
            "aggregates": {"reacted_motifs": ["Ar-I"], "formed_motifs": []},
        },
    }

    monkeypatch.setattr(
        report_mod,
        "featurize_reaction",
        lambda reaction_smiles, options=None: payloads[reaction_smiles],
    )

    report = report_mod.build_report(["rxn_ok", "rxn_unknown_low"])
    assert report["total_reactions"] == 2
    assert report["unknown_reaction_type"]["count"] == 1
    assert report["empty_reaction_key"]["count"] == 1
    assert report["low_reaction_key_quality"]["count"] == 1
    assert report["quality_levels"]["high"] == 1
    assert report["quality_levels"]["low"] == 1
    top_anomalies = dict(report["top_anomalies"])
    assert top_anomalies.get("mapping_unreliable_fallback_used") == 1
