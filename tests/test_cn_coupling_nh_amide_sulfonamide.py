import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_cn_coupling_detects_primary_amide_nucleophile() -> None:
    rxn = "CC(N)=O.O=[N+]([O-])c1ccc(Oc2ccnc(Cl)c2)c(F)c1>>CC(=O)Nc1cc(Oc2ccc([N+](=O)[O-])cc2F)ccn1"
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
    assert result.get("reaction_type") == "C_N_Coupling"
    assert "Alkyl-CONH2" in reaction_key
    assert "-> Ar-NHCOR" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_cn_coupling_detects_tosylamide_nucleophile() -> None:
    rxn = "Cc1ccc(S(N)(=O)=O)cc1.Clc1ncccn1>>Cc1ccc(S(=O)(=O)Nc2ncccn2)cc1"
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
    assert result.get("reaction_type") == "C_N_Coupling"
    assert "Ar-SO2NH2" in reaction_key
    assert "SO2NHR" in reaction_key
    assert "N(R)SO2R" not in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_cn_coupling_detects_r2ch_primary_amide_nucleophile() -> None:
    rxn = (
        "NC(=O)C1CC1.O=[N+]([O-])c1ccc(Oc2ccnc(Cl)c2)c(F)c1"
        ">>O=C(Nc1cc(Oc2ccc([N+](=O)[O-])cc2F)ccn1)C1CC1"
    )
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
    assert result.get("reaction_type") == "C_N_Coupling"
    assert "R2CH-CONH2" in reaction_key
    assert "-> Ar-NHCOR" in reaction_key
