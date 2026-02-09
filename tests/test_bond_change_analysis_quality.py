import pytest

import chemtools._atom_mapping as atom_mapping
from chemtools.featurizers.formatters import reaction as reaction_formatter
from chemtools.util.rdkit_helpers import rdkit_available
from chemtools.util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_identify_changed_atoms_dedupes_bidirectional_mapped_bonds() -> None:
    mapped_rxn = "[Cl:1][CH3:2].[OH2:3]>>[OH:3][CH3:2].[Cl-:1]"
    result = identify_changed_atoms_from_mapped_smiles(mapped_rxn)

    assert result.get("success")
    assert result.get("broken_bonds") == [(1, 2)]
    assert result.get("formed_bonds") == [(2, 3)]


def test_get_bond_change_analysis_uses_mcs_validation(monkeypatch: pytest.MonkeyPatch) -> None:
    called = {"use_mcs": None}

    def fake_hybrid(_reaction_smiles: str, use_mcs: bool = True, **_kwargs):
        called["use_mcs"] = use_mcs
        return {
            "success": True,
            "method": "hybrid",
            "agreement": {"rxnmapper_vs_mcs": True},
            "combined_confidence": 0.9,
            "recommended_result": {
                "success": True,
                "broken_bonds": [(1, "Br (leaving group)")],
                "formed_bonds": [(1, 2)],
            },
        }

    monkeypatch.setattr(reaction_formatter.rdkit_helpers, "rdkit_available", lambda: True)
    monkeypatch.setattr(atom_mapping, "analyze_bond_changes_hybrid", fake_hybrid)

    result = reaction_formatter._get_bond_change_analysis("CCBr.CS>>CCSC")

    assert called["use_mcs"] is True
    assert result is not None
    assert result.get("success") is True


def test_get_bond_change_analysis_rejects_low_conf_disagreement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_hybrid(_reaction_smiles: str, use_mcs: bool = True, **_kwargs):
        return {
            "success": True,
            "method": "hybrid",
            "agreement": {"rxnmapper_vs_mcs": False},
            "combined_confidence": 0.4,
            "recommended_result": {
                "success": True,
                "broken_bonds": [(7, "O (leaving group)")],
                "formed_bonds": [(7, 8), (8, 7)],
            },
        }

    monkeypatch.setattr(reaction_formatter.rdkit_helpers, "rdkit_available", lambda: True)
    monkeypatch.setattr(atom_mapping, "analyze_bond_changes_hybrid", fake_hybrid)

    result = reaction_formatter._get_bond_change_analysis(
        "Brc1cncnc1.O=C(O)C1CCC(F)(F)C1>>FC1(F)CCC(c2ncncc2Br)C1"
    )

    assert result is None


def test_decarboxylative_rescue_without_bond_key_selects_product_motif() -> None:
    product_motifs = [
        {"compound_id": "Ar-COR"},
        {"compound_id": "Ar-F"},
        {"compound_id": "AromN-H"},
    ]
    reactive = reaction_formatter._select_reactive_product_motifs(
        product_motifs,
        bond_key=None,
        formed_motifs=[],
        reacted_motifs=["Acyl-CO2H", "HeteroAr-H"],
        reaction_type=None,
    )
    assert "Ar-COR" in reactive


def test_fallback_mapping_used_for_product_projection_even_when_strict_bond_key_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured = {"bond_key": None}

    def fake_with_quality(_reaction_smiles: str):
        return (
            {
                "success": True,
                "broken_bonds": [(1, "O (leaving group)")],
                "formed_bonds": [(1, 2)],
            },
            True,  # unreliable => strict mapping disabled
            {
                "method": "hybrid",
                "combined_confidence": 0.4,
                "agreement": {"rxnmapper_vs_mcs": False},
                "validation": "Methods disagree - review recommended",
            },
        )

    def fake_infer_broad(*, bond_key, product_smiles):
        captured["bond_key"] = bond_key
        return []

    monkeypatch.setattr(
        reaction_formatter,
        "_get_bond_change_analysis_with_quality",
        fake_with_quality,
    )
    monkeypatch.setattr(
        reaction_formatter,
        "_infer_product_broad_tags_with_validation",
        fake_infer_broad,
    )

    result = reaction_formatter.featurize_reaction(
        "c1cc[nH]c1.O=C(O)C(=O)c1ccccc1F>>O=C(c1ccc[nH]1)c1ccccc1F",
        options={"detailed": True},
    )
    reaction_key = result.get("reaction_key") or ""

    assert "-> []" not in reaction_key
    assert "bond_formed:" not in reaction_key
    assert "bond_broken:" not in reaction_key
    assert captured["bond_key"] is not None
