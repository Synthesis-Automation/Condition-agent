import pytest

from chemtools.featurizers.formatters.reaction import _select_reactive_product_motifs
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def test_select_reactive_product_motifs_collapses_duplicate_group_representations() -> None:
    product_motifs = [
        {"id": "Ar-Alkyl", "atoms": {0, 1}},
        {"id": "Alkyl-NHCOR", "atoms": {3, 4, 5, 6, 7}, "priority": 5, "rank_score": 448.0},
        {"id": "Ar-CONHR", "atoms": {4, 5, 6, 7}, "priority": 7, "rank_score": 426.0},
    ]

    selected = _select_reactive_product_motifs(
        product_motifs,
        bond_key="break: C-O | form: C-N",
        formed_motifs=["Ar-Alkyl", "Alkyl-NHCOR", "Ar-CONHR"],
        reacted_motifs=["Ar-COOH", "Bn-NH2"],
        reaction_type="Unknown",
    )

    assert selected == ["Ar-CONHR"]


def test_select_reactive_product_motifs_uses_taxonomy_for_suzuki_products() -> None:
    product_motifs = [
        {"id": "Ar-Alkyl"},
        {"id": "Ar-Ar"},
        {"id": "Ar-CONHR"},
    ]

    selected = _select_reactive_product_motifs(
        product_motifs,
        bond_key="form: C(ar)-C; C(ar)-C(ar)",
        formed_motifs=["Ar-Alkyl", "Ar-Ar", "Ar-CONHR"],
        reacted_motifs=["Ar-Br", "Ar-B(OH)2"],
        reaction_type="Suzuki_miyaura",
    )

    assert selected == ["Ar-Alkyl", "Ar-Ar"]


def test_featurize_reaction_amide_product_key_omits_ar_alkyl() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    rxn = "O=C(O)c1cnccn1.NCc1ccc(F)cc1>>O=C(NCc1ccc(F)cc1)c1cnccn1"
    result = featurize_reaction(rxn)
    reactive = result.get("product_motifs_reactive") or []
    aggregates = result.get("aggregates") or {}

    assert "Ar-Alkyl" not in reactive
    assert any(
        motif.endswith(("-CONH2", "-CONHR", "-CONR2", "-NHCOR", "-NRCOR"))
        for motif in reactive
    )
    assert set(aggregates.get("formed_motifs_center") or []) == set(reactive)
    assert "Ar-Alkyl" in (aggregates.get("formed_motifs_context") or [])


def test_featurize_reaction_cli_like_options_detects_amide_center_motifs() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    rxn = "O=C(O)c1cnccn1.NCc1ccc(F)cc1>>O=C(NCc1ccc(F)cc1)c1cnccn1"
    result = featurize_reaction(
        rxn,
        options={
            "include_ar_h": False,
            "target_groups": None,
            "discovery_mode": True,
            "confirm_coupling_products": True,
            "motif_site_filter": "substituent",
            "detailed": True,
        },
    )
    reactive = result.get("product_motifs_reactive") or []

    assert result.get("reaction_type") == "Amide_formation"
    assert "Ar-Alkyl" not in reactive
    assert set(reactive) == {"HeteroAr-CONHR"}


def test_featurize_reaction_amide_uses_single_group_for_same_product_center() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    rxn = "COCCN.O=C(O)c1cc(I)cc(CNc2ccc(Cl)c(O)c2)c1>>COCCNC(=O)c1cc(I)cc(CNc2ccc(Cl)c(O)c2)c1"
    result = featurize_reaction(rxn, options={"detailed": True, "include_reaction_type": True})
    reaction_key = str(result.get("reaction_key") or "")
    reactive = result.get("product_motifs_reactive") or []

    assert result.get("reaction_type") == "Amide_formation"
    assert "Ar-CONHR" in reactive
    assert "Alkyl-NHCOR" not in reactive
    assert "-> Ar-CONHR" in reaction_key


def test_hydrazide_acylation_keeps_ar_alkyl_out_of_crk_reacted_list() -> None:
    if not rdkit_available():
        pytest.skip("rdkit not available")

    rxn = "NNC(=O)c1cccs1.O=C(O)Cc1cccc2c1OCCO2>>O=C(Cc1cccc2c1OCCO2)NNC(=O)c1cccs1"
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
    reaction_key = str(result.get("reaction_key") or "")
    summary = reaction_key.split(" | ", 1)[0]
    reacted = set((result.get("aggregates") or {}).get("reacted_motifs") or [])

    assert "Ar-Alkyl" not in reacted
    assert "|Ar-Alkyl|" not in summary

