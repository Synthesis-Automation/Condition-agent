import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_kumada_detects_aryl_alkyl_product() -> None:
    rxn = "Brc1ccccc1.CCCCCCCC[Mg]Br>>CCCCCCCCc1ccccc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Kumada"
    assert "|Ar-Br|RCH2-Mg*" in reaction_key
    assert "-> Ar-Alkyl" in reaction_key

