import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_alkyl_nucleophilic_substitution_alcohol() -> None:
    rxn = "CCO.CBr>>CCOC"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Alkyl_Nucleophilic_Substitution"
    assert "CH3-Br" in reaction_key
    assert "RCH2-OH" in reaction_key
    assert "RCH2-OR" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_alkyl_nucleophilic_substitution_thiol() -> None:
    rxn = "CS.CBr>>CSC"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Alkyl_Nucleophilic_Substitution"
    assert "CH3-Br" in reaction_key
    assert "CH3-SH" in reaction_key
    assert "CH3-SR" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_alkyl_nucleophilic_substitution_amine() -> None:
    rxn = "CN.CBr>>CNC"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Alkyl_Nucleophilic_Substitution"
    assert "CH3-Br" in reaction_key
    assert "CH3-NH2" in reaction_key
    assert "CH3-NHR" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_alkyl_nucleophilic_substitution_with_aryl_alcohol_nucleophile() -> None:
    rxn = "Oc1ccccc1.CBr>>COc1ccccc1"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Alkyl_Nucleophilic_Substitution"
    assert "Ar-OH" in reaction_key
    assert "Ar-OR" in reaction_key
