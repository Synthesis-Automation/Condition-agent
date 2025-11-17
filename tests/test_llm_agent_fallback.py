"""Test the deterministic fallback summary for the LLM agent CLI."""

from chem_assistant.planner.llm_agent_cli import summarize_deterministic


def test_summarize_deterministic_returns_steps() -> None:
    reaction_smiles = "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    summary = summarize_deterministic(reaction_smiles, top_k_protocols=2)
    assert summary["family"]["family"] == "cn_coupling"
    assert summary["counts"]["rule_candidates"] >= 1
    assert summary["top_protocols"], "Expected at least one top protocol"
    assert summary["top_protocols"][0]["steps"], "Protocol summary should include steps"
