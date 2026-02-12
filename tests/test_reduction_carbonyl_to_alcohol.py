import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reduction_carbonyl_to_alcohol_detected_for_aldehyde_reduction() -> None:
    rxn = "O=CC1=CC=CC=C1>>OCC1=CC=CC=C1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Reduction_carbonyl_to_alcohol"
    assert "Ar-CHO" in reaction_key
    assert "Bn-OH" in reaction_key
