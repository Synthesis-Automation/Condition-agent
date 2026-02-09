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
