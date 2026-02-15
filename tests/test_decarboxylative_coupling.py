import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_decarboxylative_coupling_detected_for_ar_x_partner() -> None:
    rxn = "O=C(O)C(=O)c1ccccc1.N#Cc1ccc(Br)cc1>>N#Cc1ccc(C(=O)c2ccccc2)cc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Decarboxylative_Coupling"
    assert "|Acyl-COOH|Ar-Br" in reaction_key


