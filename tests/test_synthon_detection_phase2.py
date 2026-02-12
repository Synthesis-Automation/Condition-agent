from __future__ import annotations

import chemtools.detection as detection
from chemtools.featurizers import detection as legacy_detection


def _patch_minimal_reaction_defs(monkeypatch) -> None:
    monkeypatch.setattr(
        detection,
        "_load_reaction_types",
        lambda: [
            {
                "id": "C_N_Coupling",
                "reactants": {
                    "electrophile": "@sp2_electrophiles",
                    "nucleophile": "@amines_nh",
                },
                "products": {},
            }
        ],
    )
    monkeypatch.setattr(
        detection,
        "_load_motif_sets",
        lambda: {
            "sp2_electrophiles": {"Ar-Br", "Ar-Cl", "Ar-I"},
            "amines_nh": {"Ar-NH2", "Ar-NHR"},
        },
    )


def test_detect_reaction_type_uses_synthon_fallback_for_missing_slot(monkeypatch) -> None:
    _patch_minimal_reaction_defs(monkeypatch)
    monkeypatch.setattr(
        detection,
        "extract_reaction_key",
        lambda _rxn: (["Ar-Br"], [], [], "CRK-v1 | Ar-Br ->"),
    )

    result = detection.detect_reaction_type("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1")
    assert result.top_match is not None
    assert result.top_match.reaction_type == "C_N_Coupling"
    assert result.top_match.electrophile == ["Ar-Br"]
    assert result.top_match.nucleophile == ["Ar-NH2"]
    assert result.top_match.slot_sources.get("nucleophile") == "synthon"
    synthon = result.synthon_evidence
    assert "Ar-NH2" in ((synthon.get("role_motifs") or {}).get("nucleophile") or [])


def test_detect_reaction_type_does_not_over_match_without_synthon_partner(monkeypatch) -> None:
    _patch_minimal_reaction_defs(monkeypatch)
    monkeypatch.setattr(
        detection,
        "extract_reaction_key",
        lambda _rxn: (["Ar-Br"], [], [], "CRK-v1 | Ar-Br ->"),
    )

    # No amine partner; should not satisfy nucleophile slot.
    result = detection.detect_reaction_type("Brc1ccccc1.CC>>Cc1ccccc1")
    assert result.matches == []


def test_detect_reaction_type_prefers_motif_when_present(monkeypatch) -> None:
    _patch_minimal_reaction_defs(monkeypatch)
    monkeypatch.setattr(
        detection,
        "extract_reaction_key",
        lambda _rxn: (["Ar-Br", "Ar-NH2"], [], [], "CRK-v1 | Ar-Br;Ar-NH2 ->"),
    )

    result = detection.detect_reaction_type("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1")
    assert result.top_match is not None
    assert result.top_match.slot_sources.get("nucleophile") == "motif"
    assert "nucleophile" not in result.top_match.synthon_slot_evidence


def test_legacy_detection_api_exposes_synthon_slot_metadata(monkeypatch) -> None:
    _patch_minimal_reaction_defs(monkeypatch)
    monkeypatch.setattr(
        detection,
        "extract_reaction_key",
        lambda _rxn: (["Ar-Br"], [], [], "CRK-v1 | Ar-Br ->"),
    )

    result = legacy_detection.detect_reaction_types("Brc1ccccc1.Nc1ccccc1>>Nc1ccccc1c1ccccc1")
    assert result.top_match is not None
    slot_evidence = result.top_match.slot_evidence
    assert (slot_evidence.get("slot_sources") or {}).get("nucleophile") == "synthon"
