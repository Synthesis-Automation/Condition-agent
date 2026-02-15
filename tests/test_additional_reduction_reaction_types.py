import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reduction_imine_to_amine_detected() -> None:
    rxn = "c1ccccc1C=Nc1ccccc1>>c1ccccc1CNc1ccccc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Reduction_imine_to_amine"
    assert "Ar-C=N" in reaction_key
    assert "Bn-NHR" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reduction_nitrile_to_amine_detected() -> None:
    rxn = "N#Cc1ccccc1>>NCc1ccccc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Reduction_nitrile_to_amine"
    assert "Ar-CN" in reaction_key
    assert "Bn-NH2" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reduction_ester_to_alcohol_detected() -> None:
    rxn = "COC(=O)c1ccccc1>>OCc1ccccc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Reduction_ester_to_alcohol"
    assert "Ar-COOR" in reaction_key
    assert "Bn-OH" in reaction_key


