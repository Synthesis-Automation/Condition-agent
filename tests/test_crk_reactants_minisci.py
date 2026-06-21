import pytest

from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


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
        token == "R_acidic-H" or token.endswith("-COOH")
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
    assert "|Acyl-COOH|" in reaction_key
    assert "HeteroAr-H" in reaction_key.split("->", 1)[0]


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_minisci_acylation_with_thiophene_keto_acid() -> None:
    rxn = "c1ccc2cnccc2c1.O=C(O)C(=O)c1cccs1>>O=C(c1cccs1)c1nccc2ccccc12"
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
    assert "->" in reaction_key
    summary = reaction_key.split(" | ", 1)[0]
    product_part = summary.split("->", 1)[1].strip()
    assert product_part != "[]"
    assert "COR" in product_part
    assert "|Acyl-COOH|" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_minisci_alkylation_keeps_heteroar_h_in_crk_reactants() -> None:
    rxn = "CC(=O)O.Cc1ccc2cc(F)ccc2n1>>Cc1cc(C)c2cc(F)ccc2n1"
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
    summary = reaction_key.split(" | ", 1)[0]
    reactant_part = summary.split("->", 1)[0].strip().lstrip("|")
    reactants = {token for token in reactant_part.split("|") if token and token != "[]"}
    assert result.get("reaction_type") == "Minisci_alkylation"
    assert "HeteroAr-H" in reactants
    assert any(token.endswith("-COOH") for token in reactants)

