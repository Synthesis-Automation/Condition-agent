import pytest

from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_crk_reacted_motifs_excludes_background_alkyl_h_for_miyaura_borylation() -> None:
    rxn = (
        "CCOC(=O)c1cnc2ccccc2c1Cl.CC1(C)COB(B2OCC(C)(C)CO2)OC1"
        ">>"
        "CCOC(=O)c1cnc2ccccc2c1B1OCC(C)(C)CO1"
    )
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "confirm_coupling_products": True,
            "motif_site_filter": "substituent",
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert "|Alkyl-H|" not in reaction_key
    assert "|HeteroAr-Cl" in reaction_key
