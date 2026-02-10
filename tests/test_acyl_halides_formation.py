import pytest

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_detects_acyl_halides_formation_from_carboxylic_acid() -> None:
    rxn = "Cc1cc(C)c([13C](=O)O)c(C)c1>>Cc1cc(C)c([13C](=O)Cl)c(C)c1"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Acyl_Halides_formation"
    assert "|Ar-CO2H -> Ar-COCl" in reaction_key


def test_validate_detection_with_crk_key_acyl_halides_formation() -> None:
    crk_raw = "|Ar-CO2H -> Ar-COCl | bond_formed: C-Cl | bond_broken: C-O | spectators: Ar-Alkyl"
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Acyl_Halides_formation"
    assert result["validation_method"] == "crk_pattern"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_detects_aliphatic_acid_chloride_as_alkyl_cocl() -> None:
    rxn = "O=C(O)CN1C(=O)c2ccccc2C1=O>>O=C(Cl)CN1C(=O)c2ccccc2C1=O"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""
    aggregates = result.get("aggregates") or {}

    assert result.get("reaction_type") == "Acyl_Halides_formation"
    assert "-> Alkyl-COCl" in reaction_key
    assert "Alkyl-COCl" in (aggregates.get("formed_motifs") or [])
