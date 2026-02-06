import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_crk_reactants_include_heteroar_h_for_minisci_like_reaction() -> None:
    rxn = "Clc1ccc2c(Cl)ccnc2c1.O=C(O)C1CCC1>>Clc1ccc2c(Cl)cc(C3CCC3)nc2c1"
    result = featurize_reaction(rxn, options={"detailed": True})
    reaction_key = result.get("reaction_key") or ""

    assert reaction_key
    summary = reaction_key.split(" | ", 1)[0]
    reactant_part = summary.split("->", 1)[0].strip().lstrip("|")
    reactants = set(token for token in reactant_part.split("|") if token and token != "[]")

    assert "HeteroAr-H" in reactants
    assert any(
        token == "R_acidic-H" or token.endswith("-CO2H")
        for token in reactants
    )
    assert result.get("reaction_type") == "Minisci_alkylation"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_minisci_detection_with_cli_like_options() -> None:
    rxn = "Clc1ccc2c(Cl)ccnc2c1.O=C(O)C1CCC1>>Clc1ccc2c(Cl)cc(C3CCC3)nc2c1"
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
    assert result.get("reaction_type") == "Minisci_alkylation"


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_minisci_acylation_detection_with_cli_like_options() -> None:
    rxn = "c1ccc2cnccc2c1.O=C(O)C(=O)c1cccc(Br)c1>>O=C(c1cccc(Br)c1)c1nccc2ccccc12"
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
    reaction_key = result.get("reaction_key") or ""
    assert result.get("reaction_type") == "Minisci_acylation"
    assert "|Acyl-CO2H|" in reaction_key
    assert "HeteroAr-H" in reaction_key.split("->", 1)[0]
