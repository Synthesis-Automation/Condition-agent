import pytest

from chemtools.featurizers.unified import featurize_molecule, featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def _motif_ids(payload: dict) -> set[str]:
    values: set[str] = set()
    for motif in payload.get("motifs") or []:
        if not isinstance(motif, dict):
            continue
        mid = motif.get("id") or motif.get("compound_id")
        if mid:
            values.add(str(mid))
    return values


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_documented_aryl_isonitrile_detected_without_discovery() -> None:
    result = featurize_molecule("[C-]#[N+]c1ccc(Cl)cc1", options={"detailed": True, "discovery_mode": False})
    motifs = _motif_ids(result)
    assert "Ar-NC" in motifs


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_documented_heteroaryl_isonitrile_detected_without_discovery() -> None:
    result = featurize_molecule("[C-]#[N+]c1ccsc1", options={"detailed": True, "discovery_mode": False})
    motifs = _motif_ids(result)
    assert "HeteroAr-NC" in motifs


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_crk_includes_aryl_isonitrile_without_discovery() -> None:
    rxn = "[C-]#[N+]c1ccc(Cl)cc1.Clc1ccc(SCc2ccccc2I)cc1>>Clc1ccc(Sc2cn(-c3ccc(Cl)cc3)c3ccccc23)cc1"
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
    assert "|Ar-I|Ar-NC ->" in reaction_key
