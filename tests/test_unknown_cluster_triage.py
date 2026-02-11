from __future__ import annotations

from scripts import triage_unknown_clusters as triage


def test_build_unknown_clusters_groups_unknown_rows(monkeypatch) -> None:
    def fake_featurize(reaction_smiles: str, options=None):
        if reaction_smiles == "A>>B":
            return {
                "reaction_type": "Unknown",
                "aggregates": {"reacted_motifs": ["Ar-I"], "formed_motifs": ["Ar-NR2"]},
                "reaction_events": {
                    "events": [{"kind": "c_n_bond_formation"}],
                    "reaction_key_quality": {"level": "high", "score_0_1": 0.9},
                    "anomalies": [],
                },
                "reaction_key": "CRK-1",
            }
        if reaction_smiles == "C>>D":
            return {
                "reaction_type": "Unknown",
                "aggregates": {"reacted_motifs": ["Ar-I"], "formed_motifs": ["Ar-NR2"]},
                "reaction_events": {
                    "events": [{"kind": "c_n_bond_formation"}],
                    "reaction_key_quality": {"level": "medium", "score_0_1": 0.6},
                    "anomalies": ["mapping_unreliable_fallback_used"],
                },
                "reaction_key": "CRK-2",
            }
        return {
            "reaction_type": "C_N_Coupling",
            "aggregates": {"reacted_motifs": ["Ar-Br"], "formed_motifs": ["Ar-NR2"]},
        }

    monkeypatch.setattr(triage, "featurize_reaction", fake_featurize)
    rows = [
        (2, {"reaction_smiles": "A>>B", "reaction_id": "r1", "source_file": "f1.csv", "reaction_type": "X"}),
        (3, {"reaction_smiles": "C>>D", "reaction_id": "r2", "source_file": "f1.csv", "reaction_type": "X"}),
        (4, {"reaction_smiles": "E>>F", "reaction_id": "r3", "source_file": "f2.csv", "reaction_type": "Y"}),
    ]
    out = triage.build_unknown_clusters(rows, sample_per_cluster=2, seed=1)
    assert out["processed_reactions"] == 3
    assert out["unknown_reactions"] == 2
    assert out["cluster_count"] == 1
    cluster = out["clusters"][0]
    assert cluster["count"] == 2
    assert cluster["motif_signature"] == "Ar-I -> Ar-NR2"
