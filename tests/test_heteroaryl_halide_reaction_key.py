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
def test_trichloropyridine_uses_heteroaryl_halide_without_aryl_fallback() -> None:
    result = featurize_molecule("Clc1cc(Cl)nc(Cl)c1", options={"detailed": True, "discovery_mode": False})
    motifs = _motif_ids(result)
    assert "HeteroAr-Cl" in motifs
    assert "Ar-Cl" not in motifs


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_key_uses_single_heteroaryl_halide_class_for_chloropyridine_partner() -> None:
    rxn = (
        "COc1ccc(CO)cc1.Clc1cc(Cl)nc(Cl)c1"
        ">>"
        "COc1ccc(COc2cc(Cl)nc(OCc3ccc(OC)cc3)c2)cc1"
    )
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "confirm_coupling_products": True,
            "motif_site_filter": "substituent",
        },
    )
    reacted = set((result.get("aggregates") or {}).get("reacted_motifs") or [])
    assert "HeteroAr-Cl" in reacted
    assert "Ar-Cl" not in reacted
