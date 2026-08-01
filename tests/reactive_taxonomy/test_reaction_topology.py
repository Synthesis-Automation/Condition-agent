"""Generic topology regressions independent of named reaction patterns."""

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction


def test_inter_and_intramolecular_scope_are_observed_from_edits() -> None:
    intermolecular = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")
    intramolecular = featurize_reaction("NCCCCBr>>C1CCCN1")
    assert intermolecular.reaction_topology is not None
    assert intramolecular.reaction_topology is not None
    assert intermolecular.reaction_topology.reaction_scope == "intermolecular"
    assert intramolecular.reaction_topology.reaction_scope == "intramolecular"
    assert intramolecular.reaction_topology.formed_ring_sizes == (5,)


def test_ring_closure_is_an_optional_pattern() -> None:
    result = featurize_reaction("NCCCCBr>>C1CCCN1")
    assert result.interpretation is not None
    assert {match.pattern_id for match in result.interpretation.pattern_matches} >= {
        "net_ring_closure",
        "net_substitution",
    }


def test_topology_serializes_inside_observation_and_signature() -> None:
    result = featurize_reaction("NCCCCBr>>C1CCCN1")
    assert result.observation is not None
    payload = asdict(result)
    assert payload["observation"]["topology"]["reaction_scope"] == "intramolecular"
    assert result.reaction_signature is not None
    assert result.reaction_signature.topology is not None
    assert result.reaction_signature.topology.reaction_scope == "intramolecular"
