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
        default="gpt-4o",
        help="LLM model for agent mode (requires API key and langchain_openai).",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.0,
        help="LLM temperature (0.0=deterministic, higher=more creative).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show agent reasoning steps and tool calls.",
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

    # LLM path with enhanced system prompt
    llm = ChatOpenAI(model=args.llm_model, temperature=args.temperature)
    tools = [
        auto_conditions_llm_tool,
        planner_detect_family_tool,
        planner_rule_candidates_tool,
        planner_protocol_candidates_tool,
        planner_hte_summary_tool,
        planner_score_candidates_tool,
        planner_fuse_tool,
    ]
    
    # Enhanced system prompt for better chemistry reasoning
    system_prompt = """You are an expert chemistry assistant specializing in organic synthesis and reaction condition optimization.

Your goal is to recommend reliable, practical reaction conditions based on:
1. Rule-based knowledge from validated reaction databases
2. Similar precedent reactions (DRFP similarity search)
3. High-throughput experimentation (HTE) statistics
4. Synthesis best practices and safety considerations

When analyzing a reaction:
- First detect the reaction family to understand the transformation
- Gather evidence from multiple sources (rules, precedents, HTE data)
- Compare and synthesize information critically
- Prioritize conditions with strong experimental support
- Explain your reasoning clearly, citing evidence sources
- Consider practical factors: scalability, cost, safety, availability
- Provide specific, actionable recommendations with quantities and conditions
- Highlight any uncertainties or alternative approaches

Always use the available tools systematically to build a comprehensive recommendation."""

    agent = create_agent(llm, tools, system_prompt=system_prompt)
    
    # Enhanced user query with clearer instructions
    user_query = f"""Analyze this reaction and recommend optimal conditions:

Reaction SMILES: {args.reaction}

Please:
1. Detect the reaction family and confidence
2. Find rule-based conditions if available
3. Search for similar precedent reactions (top {args.top_k})
4. Check HTE database for statistical insights
5. Score and rank all candidates
6. Provide your top {args.max_protocols} recommendations with:
   - Complete reaction setup (reagents, catalyst, ligand, base, solvent)
   - Specific quantities and conditions (temperature, time, concentration)
   - Success rate or yield expectations from evidence
   - Brief rationale explaining why these conditions are recommended
   - Any important considerations (safety, handling, alternatives)

Be specific and practical - a chemist should be able to run this reaction from your recommendations."""

    if args.verbose:
        print("🤖 Agent thinking...\n")
    
    response = agent.invoke({"messages": [user_query]})
    
    if args.verbose:
        print("\n" + "="*80)
        print("AGENT RESPONSE:")
        print("="*80)
    
    print(response["messages"][-1].content)
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())
