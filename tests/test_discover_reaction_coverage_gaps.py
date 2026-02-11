from __future__ import annotations

from scripts import discover_reaction_coverage_gaps as discovery


class _Slot:
    def __init__(self, allowed):
        self.allowed = list(allowed or [])


class _Defn:
    def __init__(self, name: str, category: str, reactants: dict, products: dict):
        self.name = name
        self.category = category
        self.reactants = reactants
        self.products = products


def test_build_discovery_report_counts_gaps(monkeypatch) -> None:
    defs = {
        "c_n_cross_coupling": _Defn(
            name="C-N cross coupling",
            category="coupling",
            reactants={
                "electrophile": _Slot(["Ar-Br", "Ar-I"]),
                "nucleophile": _Slot(["R-NH2"]),
            },
            products={"product": _Slot(["Ar-NR2"])},
        )
    }

    def fake_catalog():
        return defs, {"c-n cross coupling": "c_n_cross_coupling"}

    def fake_featurize(reaction_smiles: str, options=None):
        if reaction_smiles == "A>>B":
            return {
                "reaction_type": "Unknown",
                "reaction_key": "CRK-A",
                "aggregates": {"reacted_motifs": ["Ar-I"], "formed_motifs": ["Ar-NR2"]},
                "reactants": [{"motifs": [{"id": "Unclassified-Reactant"}]}],
                "products": [{"motifs": [{"id": "Ar-NR2"}]}],
                "reaction_events": {
                    "events": [{"kind": "c_n_bond_formation"}],
                    "reaction_key_quality": {"level": "medium", "score_0_1": 0.6},
                    "anomalies": [],
                },
                "detection": {"mapping_warning": {"reason": "x"}},
            }
        if reaction_smiles == "C>>D":
            return {
                "reaction_type": "c_n_cross_coupling",
                "reaction_key": "",
                "aggregates": {"reacted_motifs": ["Ar-Br"], "formed_motifs": []},
                "reactants": [{"motifs": [{"id": "Ar-Br"}]}],
                "products": [{"motifs": [{"id": "unknown"}]}],
                "reaction_events": {
                    "events": [{"kind": "c_n_bond_formation"}],
                    "reaction_key_quality": {"level": "low", "score_0_1": 0.2},
                    "anomalies": ["mapping_unreliable_fallback_used"],
                },
            }
        return {
            "reaction_type": "c_n_cross_coupling",
            "reaction_key": "CRK-E",
            "aggregates": {"reacted_motifs": ["Ar-X"], "formed_motifs": ["Ar-NR2"]},
            "reactants": [{"motifs": [{"id": "Ar-X"}]}],
            "products": [{"motifs": [{"id": "Ar-NR2"}]}],
            "reaction_events": {
                "events": [{"kind": "c_n_bond_formation"}],
                "reaction_key_quality": {"level": "high", "score_0_1": 0.9},
                "anomalies": [],
            },
        }

    monkeypatch.setattr(discovery, "load_reaction_catalog", fake_catalog)
    monkeypatch.setattr(discovery, "featurize_reaction", fake_featurize)

    rows = [
        {"reaction_smiles": "A>>B", "source_file": "f1.csv", "source_reaction_label": "X"},
        {"reaction_smiles": "C>>D", "source_file": "f1.csv", "source_reaction_label": "X"},
        {"reaction_smiles": "E>>F", "source_file": "f2.csv", "source_reaction_label": "Y"},
    ]
    report = discovery.build_discovery_report(rows, sample_per_cluster=2, progress_every=0)

    summary = report["summary"]
    assert summary["processed_reactions"] == 3
    assert summary["unknown_reaction_type"]["count"] == 1
    assert summary["empty_reaction_key"]["count"] == 1
    assert summary["low_reaction_key_quality"]["count"] == 1
    assert summary["missing_formed_motifs"]["count"] == 1
    assert summary["with_mapping_warning"]["count"] == 1
    assert summary["unresolved_reactions"]["count"] == 3

    assert report["top_unclassified_motifs"][0][0] in {"Unclassified-Reactant", "unknown"}
    outside = dict(report["taxonomy_coverage"]["top_motifs_outside_reaction_taxonomy"])
    assert outside.get("Ar-X", 0) == 1
    assert report["clusters"]
