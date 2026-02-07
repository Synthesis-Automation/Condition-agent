import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_chan_lam_cn_detects_aryl_sulfonamide_nucleophile() -> None:
    rxn = "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "discovery_mode": False,
            "include_ar_h": False,
            "motif_site_filter": "substituent",
            "confirm_coupling_products": True,
        },
    )
    reaction_key = result.get("reaction_key") or ""
    assert result.get("reaction_type") == "Chan_Lam_C_N_Coupling"
    assert "Ar-SO2NH2" in reaction_key
    assert "Ar-N(R)SO2R" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_chan_lam_cn_detects_alkyl_sulfonamide_nucleophile() -> None:
    rxn = "CS(=O)(=O)N.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)C)cc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "discovery_mode": False,
            "include_ar_h": False,
            "motif_site_filter": "substituent",
            "confirm_coupling_products": True,
        },
    )
    reaction_key = result.get("reaction_key") or ""
    assert result.get("reaction_type") == "Chan_Lam_C_N_Coupling"
    assert "Alkyl-SO2NH2" in reaction_key
    assert "Ar-N(R)SO2R" in reaction_key

