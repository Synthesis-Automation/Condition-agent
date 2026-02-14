import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


_RXN_MISSING_EQUIV = (
    "COc1nc(Cl)nc(Cl)n1.Nc1ccc(Cl)c(O)c1"
    ">>"
    "COc1nc(Cl)nc(Nc2ccc(Cl)c(Oc3nc(Cl)nc(OC)n3)c2)n1"
)


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_stoichiometry_precheck_adds_missing_repeated_reactant_by_default() -> None:
    result = featurize_reaction(
        _RXN_MISSING_EQUIV,
        options={"detailed": True, "confirm_coupling_products": True, "motif_site_filter": "substituent"},
    )

    reacted = set((result.get("aggregates") or {}).get("reacted_motifs") or [])
    assert "HeteroAr-Cl" in reacted

    detection = result.get("detection") or {}
    precheck = detection.get("reactant_precheck") or {}
    assert precheck.get("applied") is True
    assert precheck.get("added_copies") == 1


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_stoichiometry_precheck_can_be_disabled() -> None:
    result = featurize_reaction(
        _RXN_MISSING_EQUIV,
        options={
            "detailed": True,
            "confirm_coupling_products": True,
            "motif_site_filter": "substituent",
            "auto_repeat_reactants": False,
        },
    )

    reacted = set((result.get("aggregates") or {}).get("reacted_motifs") or [])
    assert "HeteroAr-Cl" not in reacted

    detection = result.get("detection") or {}
    assert "reactant_precheck" not in detection
