import pytest

from chemtools.featurizers.formatters.reaction import _select_reactive_product_motifs
from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def test_select_reactive_product_motifs_prefers_amides_over_generic_aryl_alkyl() -> None:
    product_motifs = [
        {"id": "Ar-Alkyl"},
        {"id": "Alkyl-NHCOR"},
        {"id": "Ar-CONHR"},
    ]

    selected = _select_reactive_product_motifs(
        product_motifs,
        bond_key="break: C-O | form: C-N",
        formed_motifs=["Ar-Alkyl", "Alkyl-NHCOR", "Ar-CONHR"],
        reacted_motifs=["Ar-CO2H", "Bn-NH2"],
        reaction_type="Amide_formation",
    )

    assert selected == ["Alkyl-NHCOR", "Ar-CONHR"]


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
    assert set(reactive) == {"Bn-NHCOR", "HeteroAr-CONHR"}
