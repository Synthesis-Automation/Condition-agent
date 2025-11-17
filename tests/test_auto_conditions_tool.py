"""Smoke test for the LLM-callable auto-conditions tool."""

from chem_assistant.chemtools_wrapper import auto_conditions_llm_tool
from chem_assistant.chemtools_wrapper import (
    planner_detect_family_tool,
    planner_rule_candidates_tool,
    planner_protocol_candidates_tool,
    planner_hte_summary_tool,
)


def test_auto_conditions_llm_tool_returns_protocol() -> None:
    reaction_smiles = "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    result = auto_conditions_llm_tool.invoke(
        {
            "reaction_smiles": reaction_smiles,
            "top_k_protocols": 2,
            "max_protocols": 1,
        }
    )
    assert result["success"] is True
    assert result["family"]["family"] == "cn_coupling"
    protocols = result.get("protocols") or []
    assert protocols, "Expected at least one formatted protocol"
    assert protocols[0].get("additions"), "Protocol should contain addition steps"


def test_planner_tools_individual_calls() -> None:
    reaction_smiles = "Brc1ccccc1.N1CCOCC1>>Brc1ccccc1N1CCOCC1"
    detect = planner_detect_family_tool.invoke({"reaction_smiles": reaction_smiles})
    assert detect["success"] and detect["family"] == "cn_coupling"

    rules = planner_rule_candidates_tool.invoke({"reaction_smiles": reaction_smiles})
    assert rules["success"]
    assert rules["candidates"]

    protos = planner_protocol_candidates_tool.invoke({"reaction_smiles": reaction_smiles, "top_k": 1})
    assert protos["success"]

    hte = planner_hte_summary_tool.invoke({"reaction_smiles": reaction_smiles})
    assert hte["success"]
