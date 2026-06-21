import pytest

from chemtools.molecule import featurize_molecule
from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


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
def test_documented_aryl_iodonium_detected_without_discovery() -> None:
    reagent = "COc1ccc([I+]c2ccccc2NC(C)=O)cc1"
    result = featurize_molecule(reagent, options={"detailed": True, "discovery_mode": False})
    motifs = _motif_ids(result)
    assert "Ar-Iodonium" in motifs
    assert "Ar-I" not in motifs


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_neutral_aryl_iodide_still_maps_to_ar_i() -> None:
    result = featurize_molecule("Ic1ccccc1", options={"detailed": True, "discovery_mode": False})
    motifs = _motif_ids(result)
    assert "Ar-I" in motifs


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reaction_with_iodonium_reagent_exposes_iodonium_motif() -> None:
    rxn = (
        "Oc1ccccc1I.COC1=CC=C([I+]C2=C(NC(C)=O)C=CC=C2)C=C1.[F-][B+3]([F-])([F-])[F-]"
        ">>CC(=O)Nc1ccccc1Oc1ccccc1I"
    )
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "discovery_mode": False,
            "motif_site_filter": "substituent",
            "confirm_coupling_products": True,
        },
    )

    iodonium_reactants = [
        entry
        for entry in (result.get("reactants") or [])
        if isinstance(entry, dict) and "[I+]" in str(entry.get("smiles") or "")
    ]
    assert iodonium_reactants, "Expected a parsed iodonium reactant in normalized reaction input"

    motifs = _motif_ids(iodonium_reactants[0])
    assert "Ar-Iodonium" in motifs
    assert "Ar-I" not in motifs
