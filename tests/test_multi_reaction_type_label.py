import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_multi_reaction_type_label_for_suzuki_plus_sonogashira_record() -> None:
    rxn = (
        "OB(O)c1ccc(F)cc1.Brc1sc(Br)c(Br)c1Br.C#Cc1ccc(CC)cc1"
        ">>CCc1ccc(C#Cc2c(-c3ccc(F)cc3)sc(-c3ccc(F)cc3)c2C#Cc2ccc(CC)cc2)cc1"
    )
    result = featurize_reaction(
        rxn,
        options={"detailed": True, "include_reaction_type": True, "confirm_coupling_products": True},
    )
    reaction_type = str(result.get("reaction_type") or "")
    reaction_key = str(result.get("reaction_key") or "")
    detection = result.get("detection") or {}

    assert reaction_type == "Suzuki_miyaura / Sonogashira"
    assert "reaction_types: Suzuki_miyaura+Sonogashira" in reaction_key
    assert detection.get("primary_reaction_type") == "Suzuki_miyaura"
    assert detection.get("co_reaction_types") == ["Sonogashira"]
