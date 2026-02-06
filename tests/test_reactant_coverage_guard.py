import pytest

from chemtools.featurizers.unified import featurize_reaction
from chemtools.util.rdkit_helpers import rdkit_available


def _motif_ids(entry: dict) -> set[str]:
    values: set[str] = set()
    for motif in entry.get("motifs") or []:
        if not isinstance(motif, dict):
            continue
        mid = motif.get("id") or motif.get("compound_id")
        if mid:
            values.add(str(mid))
    return values


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reactant_coverage_guard_adds_unclassified_token() -> None:
    rxn = "NotASmiles.CC(=O)O>>CC(=O)O"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "discovery_mode": True,
        },
    )
    reactants = result.get("reactants") or []
    assert reactants
    assert "Unclassified-Reactant" in _motif_ids(reactants[0])
    reaction_key = result.get("reaction_key") or ""
    assert "Unclassified-Reactant" in reaction_key


@pytest.mark.skipif(not rdkit_available(), reason="rdkit not available")
def test_reactant_coverage_guard_can_be_disabled() -> None:
    rxn = "NotASmiles.CC(=O)O>>CC(=O)O"
    result = featurize_reaction(
        rxn,
        options={
            "detailed": True,
            "discovery_mode": True,
            "reactant_coverage_guard": False,
        },
    )
    reactants = result.get("reactants") or []
    assert reactants
    assert "Unclassified-Reactant" not in _motif_ids(reactants[0])
