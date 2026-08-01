"""Regressions for observation-derived reaction patterns."""

from dataclasses import asdict

import reactive_taxonomy.reaction_patterns as pattern_module
from reactive_taxonomy import featurize_reaction, observe_reaction
from reactive_taxonomy.reaction_patterns import match_reaction_patterns


def _pattern_ids(reaction: str) -> set[str]:
    observation = observe_reaction(reaction)
    return {pattern.pattern_id for pattern in match_reaction_patterns(observation)}


def test_generic_patterns_are_derived_from_observed_edits() -> None:
    assert "net_substitution" in _pattern_ids(
        "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    )
    assert "net_elimination" in _pattern_ids("CCBr>>C=C")
    assert "net_addition" in _pattern_ids("C=C>>CC")


def test_reductive_amination_is_a_synthesis_pattern() -> None:
    reaction = "O=Cc1ccncc1.Nc1ccccc1>>c1ccc(NCc2ccncc2)cc1"
    result = featurize_reaction(reaction)

    assert result.interpretation is not None
    patterns = {item.pattern_id: item for item in result.interpretation.pattern_matches}
    assert patterns["reductive_amination_like"].tier == "synthesis"
    assert result.interpretation.primary_pattern_id == "reductive_amination_like"


def test_decarboxylative_coupling_is_optional_not_part_of_core() -> None:
    reaction = (
        "O=[C:15](O)[CH2:2][CH3:1]."
        "[cH:3]1[cH:4][c:5]([CH3:6])[n:7][c:8]2[cH:9][cH:10]"
        "[c:11]([F:12])[cH:13][c:14]12"
        ">>[CH3:1][CH2:2][c:3]1[cH:4][c:5]([CH3:6])[n:7]"
        "[c:8]2[cH:9][cH:10][c:11]([F:12])[cH:13][c:14]12"
    )
    result = featurize_reaction(reaction)

    assert result.reaction_core is not None
    assert "abstraction" not in asdict(result.reaction_core)
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == (
        "decarboxylative_coupling_like"
    )
    assert result.reaction_label is not None
    assert result.reaction_label.source == "optional_pattern"
    assert result.reaction_label.concise == "decarboxylative C–C coupling"


def test_pattern_registry_cannot_change_structural_observation(monkeypatch) -> None:
    reaction = "CCBr>>C=C"
    before = featurize_reaction(reaction)
    monkeypatch.setattr(pattern_module, "load_reaction_pattern_definitions", lambda: ())
    after = featurize_reaction(reaction)

    assert asdict(before.observation) == asdict(after.observation)
    assert before.reaction_signature is not None
    assert after.reaction_signature is not None
    assert before.reaction_signature.signature_id == after.reaction_signature.signature_id
    assert after.interpretation is not None
    assert after.interpretation.pattern_matches == ()


def test_pattern_definitions_contain_no_structural_execution_data() -> None:
    forbidden = {
        "operator_id",
        "predicted_bond_changes",
        "reconstruction_rule_id",
        "slots",
    }
    for definition in pattern_module.load_reaction_pattern_definitions():
        assert not forbidden.intersection(definition)
