"""Current structural observations can be reannotated without changing evidence."""

import json
from dataclasses import asdict

import pytest

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.observation_serialization import restore_reaction_observation
from reactive_taxonomy.reaction_api import reanalyze_reaction_observation


@pytest.mark.parametrize(
    "reaction",
    [
        "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]",
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "Brc1ccccc1.CO>>COc1ccccc1",
        "Brc1ccccc1.CS>>CSc1ccccc1",
        "CC(=O)O.CCO>>CC(=O)OCC",
        "[CH3:1][Br:1].CN>>CNC",
        "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
    ],
)
def test_json_observation_reannotation_equals_fresh_analysis(reaction):
    original = featurize_reaction(reaction)
    assert original.observation is not None
    payload = json.loads(json.dumps(asdict(original.observation)))
    restored = restore_reaction_observation(payload)
    actual = reanalyze_reaction_observation(restored)
    assert actual.to_dict() == original.to_dict()
    assert json.loads(json.dumps(asdict(restored))) == payload


def test_restoration_rejects_stale_schema_and_unknown_fields():
    original = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    payload = asdict(original.observation)
    with pytest.raises(ValueError, match="schema"):
        restore_reaction_observation({**payload, "schema_version": "stale"})
    with pytest.raises(ValueError, match="Unexpected"):
        restore_reaction_observation({**payload, "arbitrary_type": "executable"})
