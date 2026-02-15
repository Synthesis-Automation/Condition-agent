from __future__ import annotations

import pandas as pd

from app import A_convert_to_hte_format as hte_convert


def test_cleanup_reaction_smiles_removes_coordination_and_counterions() -> None:
    raw = "Cl.Clc1ccncc1.CN(C)->[Zn+2]([Cl-])([Cl-])<-N(C)C>>CN(C)c1ccncc1.Cl"
    cleaned, stats = hte_convert._cleanup_reaction_smiles_for_featurization(raw)

    assert cleaned == "Clc1ccncc1>>CN(C)c1ccncc1"
    assert stats["cleanup_applied"] == 1
    assert stats["coordination_removed"] == 1
    assert stats["counterion_removed"] == 2


def test_scope_filter_marks_explicit_cc_as_out_of_scope_for_cn_dataset() -> None:
    out_of_scope, reason = hte_convert._is_out_of_scope_for_dataset(
        "C_N_Coupling",
        {
            "formed_bond_classes": ["C-C"],
            "event_kinds": ["c_c_bond_formation"],
        },
    )
    assert out_of_scope is True
    assert reason == "formed_bond_class_conflict"


def test_scope_filter_keeps_explicit_cn_in_scope_for_cn_dataset() -> None:
    out_of_scope, reason = hte_convert._is_out_of_scope_for_dataset(
        "C_N_Coupling",
        {
            "formed_bond_classes": ["C-N"],
            "event_kinds": ["c_n_bond_formation"],
        },
    )
    assert out_of_scope is False
    assert reason in {"expected_bond_class_present", "expected_event_kind_present"}


def test_scope_filter_uses_formed_motifs_for_cn_guardrails() -> None:
    out_of_scope, reason = hte_convert._is_out_of_scope_for_dataset(
        "C_N_Coupling",
        {},
        {"Ar-Ar", "HeteroAr-H"},
    )
    assert out_of_scope is True
    assert reason == "formed_motif_conflict_ar_ar"

    kept, kept_reason = hte_convert._is_out_of_scope_for_dataset(
        "C_N_Coupling",
        {"formed_bond_classes": ["C-C"]},
        {"HeteroAr-NHR"},
    )
    assert kept is False
    assert kept_reason == "expected_n_formed_motif_present"


def test_process_dataset_drops_explicit_out_of_scope_rows(tmp_path, monkeypatch) -> None:
    input_path = tmp_path / "C_N_Coupling.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame(
        [
            {
                "reaction_smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
                "yield_pct": 52,
            }
        ]
    ).to_csv(input_path, index=False)

    monkeypatch.setattr(
        hte_convert,
        "cached_featurize",
        lambda smiles: {"motifs": [{"id": "Ar-Br"}], "context_motifs": []},
    )
    monkeypatch.setattr(
        hte_convert,
        "cached_featurize_reaction",
        lambda smiles, llm_assist_signature="": {
            "reaction_type": "Unknown",
            "reaction_key": "|Ar-B(OH)2|Ar-Br -> Ar-Ar | bond_formed: C-C | bond_broken: C-Br",
            "reaction_events": {
                "events": [{"kind": "c_c_bond_formation", "confidence": 0.95}],
                "bond_pairs": {"formed": [["C", "C"]], "broken": [["C", "Br"]]},
            },
            "aggregates": {
                "reacted_motifs": ["Ar-B(OH)2", "Ar-Br"],
                "formed_motifs": ["Ar-Ar"],
                "spectator_motifs": [],
            },
        },
    )

    hte_convert.process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
    )

    assert not output_path.exists()


def test_process_dataset_keeps_cn_when_n_formed_motif_present(tmp_path, monkeypatch) -> None:
    input_path = tmp_path / "C_N_Coupling.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame(
        [
            {
                "reaction_smiles": "Clc1ccncc1.CNC>>CNCc1ccncc1",
                "yield_pct": 48,
            }
        ]
    ).to_csv(input_path, index=False)

    monkeypatch.setattr(
        hte_convert,
        "cached_featurize",
        lambda smiles: {"motifs": [{"id": "HeteroAr-Cl"}], "context_motifs": []},
    )
    monkeypatch.setattr(
        hte_convert,
        "cached_featurize_reaction",
        lambda smiles, llm_assist_signature="": {
            "reaction_type": "Unknown",
            "reaction_key": "|HeteroAr-Cl -> HeteroAr-NHR | bond_formed: C(ar)-C(ar) | bond_broken: C(ar)-Cl",
            "reaction_events": {
                "events": [{"kind": "c_c_bond_formation", "confidence": 0.9}],
                "bond_pairs": {"formed": [["C", "C"]], "broken": [["C", "Cl"]]},
            },
            "aggregates": {
                "reacted_motifs": ["HeteroAr-Cl"],
                "formed_motifs": ["HeteroAr-NHR"],
                "spectator_motifs": [],
            },
        },
    )

    hte_convert.process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
    )

    out = pd.read_csv(output_path)
    assert len(out) == 1


def test_scope_filter_for_co_dataset() -> None:
    blocked, blocked_reason = hte_convert._is_out_of_scope_for_dataset(
        "C_O_Coupling",
        {"formed_bond_classes": ["C-C"], "event_kinds": ["c_c_bond_formation"]},
    )
    assert blocked is True
    assert blocked_reason == "formed_bond_class_conflict"

    kept, kept_reason = hte_convert._is_out_of_scope_for_dataset(
        "C_O_Coupling",
        {"formed_bond_classes": ["C-O"]},
    )
    assert kept is False
    assert kept_reason == "expected_bond_class_present"


def test_scope_filter_for_cs_dataset() -> None:
    blocked, blocked_reason = hte_convert._is_out_of_scope_for_dataset(
        "C_S_Coupling",
        {"event_kinds": ["c_n_bond_formation"]},
    )
    assert blocked is True
    assert blocked_reason == "event_kind_conflict"

    kept, kept_reason = hte_convert._is_out_of_scope_for_dataset(
        "C_S_Coupling",
        {"event_kinds": ["c_s_bond_formation"]},
    )
    assert kept is False
    assert kept_reason == "expected_event_kind_present"
