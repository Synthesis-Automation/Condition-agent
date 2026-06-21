import pytest

from chemtools.reaction import featurize_reaction
from chemtools.core.rdkit import rdkit_available


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_broad_definition_group_hydrazine_marked_reacted_on_acylation() -> None:
    rxn = "O=C(O)c1ccccc1F.NNc1ccccc1>>O=C(NNc1ccccc1)c1ccccc1F"
    result = featurize_reaction(rxn, options={"detailed": True, "include_reaction_type": True})
    reaction_key = str(result.get("reaction_key") or "")
    aggregates = result.get("aggregates") or {}

    reacted = set(aggregates.get("reacted_motifs") or [])
    formed = set(aggregates.get("formed_motifs") or [])
    spectator = set(aggregates.get("spectator_motifs") or [])

    assert "Ar-Hydrazine" in reacted
    assert "Ar-Hydrazine" not in spectator
    assert "Ar-Hydrazide" in formed
    assert "|Ar-COOH|Ar-Hydrazine -> Ar-Hydrazide" in reaction_key

