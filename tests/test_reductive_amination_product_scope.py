import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reductive_amination_detects_benzylic_amine_product() -> None:
    rxn = "O=Cc1ccc(CO)o1.CCCCCCCCN>>CCCCCCCCNCc1ccc(CO)o1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_reaction_type": True,
            "confirm_coupling_products": True,
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert result.get("reaction_type") == "Reductive_amination"
    assert "|Ar-CHO|RCH2-NH2" in reaction_key
    assert "-> Bn-NHR|RCH2-NHR" in reaction_key

