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
    assert actions[:6] == [
        "reaction_diff",
        "validate",
        "confidence_calibration",
        "llm_rerank",
        "precedent_lookup",
        "finalize",
    ]
    assert result.tool_artifacts["fallback_candidates"]["status"] == "skipped"
    assert result.tool_artifacts["confidence_calibration"]["status"] == "not_implemented"
    assert result.tool_artifacts["llm_rerank"]["status"] == "not_implemented"
    assert result.tool_artifacts["precedent_lookup"]["status"] == "not_implemented"
    assert result.evidence["diff"]["principal_pair"]
    assert result.coverage_suggestions == []


def test_foundation_gateway_fallback_recovers_snar_candidate() -> None:
    gateway = ReactionAgentGateway(min_confidence=0.5)
    reaction = "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1"
    result = gateway.run(reaction)

    assert result.status == "completed"
    assert result.final_decision["reaction_type"] == "SNAr_CN"
    assert float(result.final_decision["confidence"]) >= 0.8
    actions = [row.action for row in result.trace]
    assert "reaction_diff" in actions
    assert "fallback_candidates" in actions
    assert "confidence_calibration" in actions
    assert "llm_rerank" in actions
    assert "precedent_lookup" in actions
    assert "coverage" not in actions
    assert result.tool_artifacts["fallback_candidates"]["status"] == "ok"
    assert result.tool_artifacts["fallback_candidates"]["candidate_reaction_types"][0] == "SNAr_CN"
    assert result.evidence["detection"]["reaction_key"]
    assert result.coverage_suggestions == []


def test_foundation_gateway_unknown_still_triggers_coverage_advice() -> None:
    gateway = ReactionAgentGateway(min_confidence=0.5)
    reaction = "CCO.CN>>CCN"
    result = gateway.run(reaction)

    assert result.status == "completed"
    assert result.final_decision["reaction_type"] == "unknown"
    actions = [row.action for row in result.trace]
    assert "fallback_candidates" in actions
    assert "coverage" in actions
    assert result.tool_artifacts["fallback_candidates"]["status"] == "ok"
    assert result.tool_artifacts["fallback_candidates"]["candidate_reaction_types"] == []
    assert len(result.coverage_suggestions) >= 1
