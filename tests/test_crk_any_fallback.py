import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_cs_coupling_uses_specific_thiol_over_any_fallback() -> None:
    rxn = "COc1ccc(Br)cc1.CCS>>COc1ccc(SCC)cc1"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "include_ar_h": False,
            "motif_site_filter": "substituent",
            "discovery_mode": False,
            "confirm_coupling_products": True,
        },
    )

    reaction_key = result.get("reaction_key") or ""
    assert reaction_key
    summary = reaction_key.split(" | ", 1)[0]
    reactant_part = summary.split("->", 1)[0].strip().lstrip("|")
    reactants = {token for token in reactant_part.split("|") if token and token != "[]"}

    assert "RCH2-SH" in reactants
    assert "Any-SH" not in reactants
    assert "|Any-SH|" not in reaction_key
