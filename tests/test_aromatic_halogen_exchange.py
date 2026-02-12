import pytest

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def test_validate_detection_with_crk_key_aromatic_halogen_exchange() -> None:
    crk_raw = "|Ar-Cl -> Ar-F | bond_formed: C(ar)-F | bond_broken: C(ar)-Cl"
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Aromatic_Halogen_Exchange"
    assert result["validation_method"] == "crk_pattern"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_detects_aromatic_halogen_exchange_from_reaction_smiles() -> None:
    rxn = "O=[N+]([O-])c1ccc(Cl)cc1>>O=[N+]([O-])c1ccc(F)cc1"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert result.get("reaction_type") == "Aromatic_Halogen_Exchange"
    assert "|Ar-Cl -> Ar-F" in reaction_key


def test_halogenation_aromatic_not_overwritten_by_exchange_rule() -> None:
    crk_raw = "|Ar-H -> Ar-Cl | bond_formed: C(ar)-Cl"
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key=crk_raw,
    )
    assert result["reaction_type"] == "Halogenation_aromatic"
