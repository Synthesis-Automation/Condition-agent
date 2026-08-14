"""Tests for the Qt-free web featurization presentation service."""

from __future__ import annotations

from app.web_api.features import analyze_features, detect_input_kind


MAPPED_INTERCOMPONENT_ANNULATION = (
    "CS(=O)(=O)O[CH2:1][CH2:2][c:3]1[c:9](Cl)"
    "[cH:8][cH:7][cH:6][cH:4]1.[NH2:28][CH3:26]>>"
    "[CH3:26][N:28]1[CH2:1][CH2:2][c:3]2[c:9]1"
    "[cH:8][cH:7][cH:6][cH:4]2"
)


def test_detect_input_kind_supports_molecules_and_reactions() -> None:
    assert detect_input_kind("Brc1ccccc1") == "molecule"
    assert detect_input_kind("Brc1ccccc1.B(O)O>>c1ccccc1") == "reaction"


def test_molecule_feature_analysis_returns_compact_and_full_results() -> None:
    result = analyze_features("Brc1ccc(N)cc1C#N")

    assert result["input_kind"] == "molecule"
    assert result["valid"] is True
    assert result["overview"]["atom_count"] == 10
    assert result["overview"]["motif_count"] == len(result["motifs"])
    assert result["analysis"]["structure"]["canonical_smiles"]


def test_reaction_feature_analysis_includes_core_projection() -> None:
    result = analyze_features(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    assert result["input_kind"] == "reaction"
    assert result["valid"] is True
    assert result["reaction_core"] is not None
    assert result["reaction_core"]["event_count"] >= 1
    assert result["core_projection"] is not None
    assert "<svg" in result["core_graphic_svg"]
    assert result["reactive_sites"]
    assert all(site["side"] == "reactant" for site in result["reactive_sites"])


def test_web_core_graphic_preserves_intercomponent_annulation_ring() -> None:
    result = analyze_features(MAPPED_INTERCOMPONENT_ANNULATION)

    assert result["valid"] is True
    assert result["core_graphic_error"] is None
    assert "<svg" in result["core_graphic_svg"]
    assert result["core_projection"]["reaction_scope"] == "intermolecular"
    assert result["core_projection"]["minimum_reaction_smiles"] == (
        "*COS(C)(=O)=O.*c1ccccc1Cl.*N>>*CN(*)c1ccccc1*"
    )
    assert result["core_graphic_svg"].count("stroke-dasharray") == 2
