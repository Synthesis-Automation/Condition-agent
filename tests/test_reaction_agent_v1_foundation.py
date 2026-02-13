"""Tests for reaction agent v1 foundation loop."""

from reaction_agent_v1 import ReactionAgentGateway


def test_foundation_gateway_classifies_known_reaction() -> None:
    gateway = ReactionAgentGateway(min_confidence=0.5)
    reaction = "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
    result = gateway.run(reaction)

    assert result.status == "completed"
    assert result.final_decision["reaction_type"] == "Amide_formation"
    assert float(result.final_decision["confidence"]) >= 0.9
    actions = [row.action for row in result.trace]
    assert actions[:3] == ["analyze", "validate", "finalize"]
    assert result.coverage_suggestions == []


def test_foundation_gateway_unknown_triggers_coverage_advice() -> None:
    gateway = ReactionAgentGateway(min_confidence=0.5)
    reaction = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = gateway.run(reaction)

    assert result.status == "completed"
    assert result.final_decision["reaction_type"] == "unknown"
    actions = [row.action for row in result.trace]
    assert "coverage" in actions
    assert len(result.coverage_suggestions) >= 1
    suggestion_ids = {row["suggestion_id"] for row in result.coverage_suggestions}
    assert "tool-candidate-retrieval-gap" in suggestion_ids
    assert "taxonomy-snar-cn-gap" in suggestion_ids
