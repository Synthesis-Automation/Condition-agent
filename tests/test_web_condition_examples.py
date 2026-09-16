"""The displayed examples must retain complete, structure-supported chemistry."""

import json
from pathlib import Path

import pytest
from rdkit import Chem

from condition_recommender.reaction_completion import propose_reaction_completion
from reactive_taxonomy import featurize_reaction


EXAMPLES = json.loads(
    (
        Path(__file__).resolve().parents[1]
        / "web/reaction_recommender/src/conditions/examples.json"
    ).read_text(encoding="utf-8")
)


@pytest.mark.parametrize("example", EXAMPLES, ids=lambda item: item["label"])
def test_web_example_is_complete_and_structure_supported(example: dict) -> None:
    """Exercise the exact UI structures, including protected and heteroaryl sites."""
    reaction = example["reaction_smiles"]
    analysis = featurize_reaction(reaction)
    assert analysis.valid, analysis.error
    assert analysis.reaction_signature is not None
    assert analysis.reaction_core is not None
    assert not propose_reaction_completion(reaction).requirements
    reactants, _, product = reaction.split(">")
    assert len(reactants.split(".")) == 2
    molecule = Chem.MolFromSmiles(product)
    assert molecule is not None
    assert molecule.GetNumHeavyAtoms() >= 15
