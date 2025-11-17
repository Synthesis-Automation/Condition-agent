"""
LLM-assisted auto-conditions entrypoint with deterministic fallback.

If an LLM client (e.g., OpenAI) and API key are available, the agent will use
the `auto_conditions_llm_tool` to drive tool calls. If not, it falls back to
the deterministic auto-conditions pipeline and prints a concise summary.
"""

from __future__ import annotations

import argparse
import os
from typing import Any, Dict, List

from chem_assistant.planner import ReactionInput, auto_conditions
from chem_assistant.chemtools_wrapper import (
    auto_conditions_llm_tool,
    planner_detect_family_tool,
    planner_rule_candidates_tool,
    planner_protocol_candidates_tool,
    planner_hte_summary_tool,
    planner_score_candidates_tool,
    planner_fuse_tool,
)

try:  # Optional: only used when key + deps present
    from langchain_openai import ChatOpenAI
    from langchain.agents import create_agent
    LLM_AVAILABLE = True
except Exception:
    LLM_AVAILABLE = False
    ChatOpenAI = None  # type: ignore
    create_agent = None  # type: ignore


def summarize_deterministic(reaction_smiles: str, top_k_protocols: int = 5) -> Dict[str, Any]:
    """Run deterministic auto-conditions and produce a concise summary payload."""
    result = auto_conditions(
        ReactionInput(reaction_smiles=reaction_smiles),
        top_k_protocols=top_k_protocols,
        max_protocols=3,
        build_protocols=True,
    )
    summary: Dict[str, Any] = {
        "family": result.family.model_dump(),
        "counts": {
            "rule_candidates": len(result.rule_candidates),
            "protocol_candidates": len(result.protocol_candidates),
        },
        "hte": result.hte_summary.model_dump() if result.hte_summary else None,
        "top_protocols": [],
    }
    for proto in result.protocols:
        summary["top_protocols"].append(
            {
                "candidate_id": proto.candidate_id,
                "source": proto.source,
                "steps": proto.additions[:6],
            }
        )
    return summary


def run_cli() -> int:
    parser = argparse.ArgumentParser(
        description="LLM-assisted auto-conditions (with deterministic fallback)."
    )
    parser.add_argument(
        "--reaction",
        required=True,
        help="Reaction SMILES (reactants>>products).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of DRFP protocol candidates to retrieve.",
    )
    parser.add_argument(
        "--max-protocols",
        type=int,
        default=3,
        help="Maximum protocols to format.",
    )
    parser.add_argument(
        "--llm-model",
        default="gpt-4o-mini",
        help="LLM model for agent mode (requires API key and langchain_openai).",
    )
    args = parser.parse_args()

    use_llm = LLM_AVAILABLE and bool(os.getenv("OPENAI_API_KEY"))
    if not use_llm:
        summary = summarize_deterministic(args.reaction, args.top_k)
        print("[fallback] LLM unavailable; using deterministic pipeline.")
        print(f"Family: {summary['family'].get('family')} (confidence={summary['family'].get('confidence')})")
        print(f"Rule candidates: {summary['counts']['rule_candidates']}, protocol candidates: {summary['counts']['protocol_candidates']}")
        if summary["hte"]:
            print(f"HTE type: {summary['hte'].get('reaction_type')}")
        print("Top protocols:")
        for proto in summary["top_protocols"]:
            print(f"- {proto['candidate_id']} ({proto['source']})")
            for step in proto["steps"]:
                print(f"    • {step}")
        return 0

    # LLM path
    llm = ChatOpenAI(model=args.llm_model, temperature=0)
    tools = [
        auto_conditions_llm_tool,
        planner_detect_family_tool,
        planner_rule_candidates_tool,
        planner_protocol_candidates_tool,
        planner_hte_summary_tool,
        planner_score_candidates_tool,
        planner_fuse_tool,
    ]
    agent = create_agent(llm, tools)
    user_query = (
        "Recommend reaction conditions with rationale and provide up to 3 protocol steps. "
        "Use planner tools to detect family, gather rule/precedent/HTE signals, "
        "score/fuse candidates, and fall back to auto_conditions_llm_tool if needed. "
        f"Reaction: {args.reaction}"
    )
    response = agent.invoke({"messages": [user_query]})
    print(response["messages"][-1].content)
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())
