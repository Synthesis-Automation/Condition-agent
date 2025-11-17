"""
Smoke tests for the planner API stubs.

These tests ensure the deterministic spine (family detect → rules → DRFP → HTE)
can be invoked without raising, even when optional indexes are missing.
"""

from chem_assistant.planner.api import (
    ReactionInput,
    detect_family,
    auto_conditions,
    fetch_rule_candidates,
    find_similar_protocols,
    fetch_hte_stats,
    fuse_scores,
    plan_workflow,
    CandidateCondition,
    EvidenceScores,
)


def test_planner_chain_runs_for_cn_coupling() -> None:
    reaction = ReactionInput(
        reaction_smiles="Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    )

    family = detect_family(reaction)
    assert family.family == "cn_coupling"

    rule_candidates = fetch_rule_candidates(reaction, family)
    assert rule_candidates, "Rule candidates should be present for C-N coupling"
    assert rule_candidates[0].components

    protocol_candidates = find_similar_protocols(reaction, top_k=2)
    assert isinstance(protocol_candidates, list)

    hte_summary = fetch_hte_stats(reaction, family)
    if hte_summary:
        assert hte_summary.reaction_type in (None, "C_N_Coupling")

    fused = fuse_scores(rule_candidates, protocol_candidates, hte_summary)
    assert fused.ranked
    assert fused.ranked[0].candidate_id == rule_candidates[0].candidate_id


def test_plan_workflow_outputs_steps() -> None:
    reaction = ReactionInput(
        reaction_smiles="Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    )
    plan = plan_workflow(reaction)
    step_names = {step.name for step in plan.steps}
    assert {"rule_candidates", "drfp_protocols", "hte_stats", "ml_score"} <= step_names


def test_auto_conditions_pipeline_returns_protocols() -> None:
    reaction = ReactionInput(
        reaction_smiles="Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    )
    result = auto_conditions(reaction, top_k_protocols=2, max_protocols=2)

    assert result.family.family == "cn_coupling"
    assert result.fused.ranked, "Fused ranking should not be empty"
    assert result.protocols, "Protocols should be built from top candidates"
    assert result.protocols[0].additions, "Protocol should contain addition steps"


def test_fusion_prefers_ml_boosted_rule_candidate() -> None:
    rule_cand = CandidateCondition(
        candidate_id="rule::a",
        components={"ligand": "XPhos"},
        source="rule",
        raw_score=0.2,
        metadata={},
    )
    proto_cand = CandidateCondition(
        candidate_id="protocol::b",
        components={"solvent": "DMSO"},
        source="protocol",
        raw_score=0.6,
        metadata={},
    )
    ml_scores = [
        EvidenceScores(candidate_id="rule::a", scores={"ml": 0.9}),
        EvidenceScores(candidate_id="protocol::b", scores={"ml": 0.1}),
    ]
    fused = fuse_scores([rule_cand], [proto_cand], None, ml_scores=ml_scores)
    assert fused.ranked[0].candidate_id == "rule::a"
    assert fused.provenance["composite_scores"]["rule::a"] > fused.provenance["composite_scores"]["protocol::b"]
