import pytest

from chemtools.reaction.featurize import _select_reactive_product_motifs
from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


def test_select_reactive_product_motifs_prefers_reacted_heteroaryl_scaffold() -> None:
    product_motifs = [
        {"id": "Ar-OR", "atoms": {0, 1}},
        {"id": "HeteroAr-OR", "atoms": {0, 1}, "priority": 7},
    ]

    selected = _select_reactive_product_motifs(
        product_motifs,
        bond_key="break: C(ar)-Cl | form: C(ar)-O",
        formed_motifs=["Ar-OR"],
        reacted_motifs=["HeteroAr-Cl", "RCH2-OH"],
        reaction_type="C_O_Coupling",
    )

    assert selected == ["HeteroAr-OR"]


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_key_prefers_heteroaryl_or_for_heteroaryl_chloride_coupling() -> None:
    rxn = (
        "C=C(CO)COCc1ccccc1.COc1cc(-c2nnc(Cl)c3ccccc23)cc(OC)c1OC"
        ">>"
        "C=C(COCc1ccccc1)COc1nnc(-c2cc(OC)c(OC)c(OC)c2)c2ccccc12"
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
    assert "|HeteroAr-Cl|RCH2-OH -> HeteroAr-OR" in reaction_key
    reactive = result.get("product_motifs_reactive") or []
    assert "HeteroAr-OR" in reactive
    assert "Ar-OR" not in reactive


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_key_keeps_ar_or_for_plain_aryl_chloride_coupling() -> None:
    rxn = "Oc1ccccc1.Clc1ccccc1>>Oc1ccccc1Oc1ccccc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "confirm_coupling_products": True,
            "motif_site_filter": "substituent",
        },
    )
    reaction_key = str(result.get("reaction_key") or "")
    assert "-> Ar-OR" in reaction_key
