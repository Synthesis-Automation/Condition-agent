from __future__ import annotations

import json

import pandas as pd

from app import A_convert_to_hte_format as hte_convert


def _fake_reaction_bundle(*, taxonomy_gap_suspected: bool) -> dict:
    return {
        "reaction_type": "Unknown",
        "reaction_key": "|Ar-Cl -> AromN-H | bond_formed: C-N | bond_broken: C-Cl",
        "aggregates": {
            "reacted_motifs": ["Ar-Cl"],
            "formed_motifs": ["AromN-H"],
            "spectator_motifs": [],
        },
        "meta": {
            "llm_assist": {
                "enabled": True,
                "used": True,
                "status": "ok",
                "decision": "no_valid_suggestion",
                "provider": "openai",
                "model": "gpt-test",
                "suggested_reaction_type": "Unknown",
                "suggested_confidence": 0.2,
                "requires_human_review": True,
                "uncertainty_flags": ["taxonomy_boundary_case"],
                "mechanistic_family": "SNAr",
                "mechanistic_rationale": "C-Cl displacement with C-N formation",
                "tautomer_or_representation_issue": True,
                "taxonomy_gap_suspected": taxonomy_gap_suspected,
                "taxonomy_gap_note": "Hydrazine-like outcome not explicitly covered",
                "deterministic_checks_used": ["reaction_key", "event_kinds"],
                "uncertainty_reasons": ["unknown_validated_detection"],
            }
        },
    }


def test_process_dataset_exports_taxonomy_gap_sidecars_when_flagged(
    tmp_path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame([{"reaction_smiles": "Clc1ccccc1>>Nc1ccccc1", "yield_pct": 55}]).to_csv(
        input_path, index=False
    )

    monkeypatch.setattr(
        hte_convert,
        "cached_featurize",
        lambda smiles: {"motifs": [{"id": "Ar-Cl"}], "context_motifs": []},
    )
    monkeypatch.setattr(
        hte_convert,
        "cached_featurize_reaction",
        lambda smiles, llm_assist_signature="": _fake_reaction_bundle(
            taxonomy_gap_suspected=True
        ),
    )

    hte_convert.process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
        llm_assist_options={"enabled": True, "provider": "openai", "model": "gpt-test"},
    )

    json_sidecar = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.json")
    csv_sidecar = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.csv")
    assert json_sidecar.exists()
    assert csv_sidecar.exists()

    payload = json.loads(json_sidecar.read_text(encoding="utf-8"))
    assert len(payload) == 1
    assert payload[0]["taxonomy_gap_suspected"] is True
    assert payload[0]["mechanistic_family"] == "SNAr"

    csv_payload = pd.read_csv(csv_sidecar)
    assert len(csv_payload) == 1
    assert bool(csv_payload.iloc[0]["taxonomy_gap_suspected"]) is True


def test_process_dataset_skips_taxonomy_gap_sidecars_when_not_flagged(
    tmp_path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "input.csv"
    output_path = tmp_path / "output.csv"
    pd.DataFrame([{"reaction_smiles": "Clc1ccccc1>>Nc1ccccc1", "yield_pct": 55}]).to_csv(
        input_path, index=False
    )

    monkeypatch.setattr(
        hte_convert,
        "cached_featurize",
        lambda smiles: {"motifs": [{"id": "Ar-Cl"}], "context_motifs": []},
    )
    monkeypatch.setattr(
        hte_convert,
        "cached_featurize_reaction",
        lambda smiles, llm_assist_signature="": _fake_reaction_bundle(
            taxonomy_gap_suspected=False
        ),
    )

    hte_convert.process_reaction_dataset(
        str(input_path),
        str(output_path),
        drop_no_catalyst=False,
        llm_assist_options={"enabled": True, "provider": "openai", "model": "gpt-test"},
    )

    json_sidecar = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.json")
    csv_sidecar = output_path.with_name(f"{output_path.stem}.taxonomy_gaps.csv")
    assert not json_sidecar.exists()
    assert not csv_sidecar.exists()
