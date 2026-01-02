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
    # Core analysis tools
    detect_reaction_family_tool,
    classify_reactant_tool,
    get_functional_groups_tool,
    analyze_bond_changes_tool,
    
    # Recommendation tools from all databases
    recommend_conditions_tool,          # ML-based precedent search
    rule_based_conditions_tool,         # Rule database
    unified_recommender_tool,           # Unified dataset/protocol/HTE search
    search_precedents_tool,             # Precedent database
    protocol_recommendation_tool,       # Protocol database
    hte_recommend_tool,                 # HTE database (66K experiments)
    hte_analytics_tool,                 # HTE analytics
    
    # Supporting tools
    list_supported_cores_tool,          # Catalyst enumeration
    reaction_dataset_analytics_tool,    # Dataset statistics
    
    # Planner tools (compact wrappers)
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

    # LLM path with comprehensive tool set
    llm = ChatOpenAI(model=args.llm_model, temperature=args.temperature)
    tools = [
        # Phase 1: Structural Analysis
        detect_reaction_family_tool,
        classify_reactant_tool,
        get_functional_groups_tool,
        analyze_bond_changes_tool,
        
        # Phase 2: Multi-Database Search
        recommend_conditions_tool,          # ML precedent database
        rule_based_conditions_tool,         # Rule database
        unified_recommender_tool,           # Unified dataset/protocol/HTE search
        search_precedents_tool,             # Precedent reactions
        protocol_recommendation_tool,       # Full protocols
        hte_recommend_tool,                 # HTE experiments (66K+)
        hte_analytics_tool,                 # HTE statistics
        
        # Phase 3: Analysis & Synthesis
        list_supported_cores_tool,
        reaction_dataset_analytics_tool,
        
        # Fallback: All-in-one tool
        auto_conditions_llm_tool,
    ]
    
    # Comprehensive system prompt for thorough analysis
    system_prompt = """You are an expert chemistry assistant specializing in organic synthesis and reaction condition optimization.

Your mission is to provide COMPREHENSIVE, evidence-based recommendations by systematically gathering information from ALL available sources.

## Your Analytical Process (FOLLOW THIS ORDER):

### Phase 1: STRUCTURAL ANALYSIS (Required)
1. Use detect_reaction_family_tool to identify the reaction type
2. Use classify_reactant_tool to analyze EACH reactant's chemical class
3. Use get_functional_groups_tool to identify key functional groups in reactants
4. Use analyze_bond_changes_tool to understand the transformation

### Phase 2: MULTI-DATABASE SEARCH (Query ALL sources)
You must gather evidence from EVERY available database:

A. Rule Database:
   - Use rule_based_conditions_tool for validated rule-based conditions
   
B. ML Precedent Database:
   - Use recommend_conditions_tool for ML-ranked precedent reactions
   - Use search_precedents_tool for detailed precedent exploration
   
C. Protocol Database:
   - Use protocol_recommendation_tool for complete experimental protocols
   - Use unified_recommender_tool for unified dataset/protocol/HTE search
   
D. HTE Database (66,308 experiments):
   - Use hte_recommend_tool for HTE-based recommendations with catalyst filtering
   - Use hte_analytics_tool to get statistical insights (top catalysts, success rates)

E. Dataset Analytics:
   - Use reaction_dataset_analytics_tool to understand frequency/yield patterns
   - Use list_supported_cores_tool to enumerate validated catalyst options

### Phase 3: SYNTHESIS & RANKING
After collecting evidence from ALL sources:
1. Compare recommendations across databases
2. Identify consensus (conditions appearing in multiple sources)
3. Note contradictions or alternatives
4. Rank based on: evidence strength, success rates, practical feasibility
5. Consider: catalyst availability, cost, safety, scalability

### Phase 4: FINAL RECOMMENDATIONS
Provide 2-3 top recommendations with:
- Complete reaction setup (specific reagents, quantities, conditions)
- Evidence sources (which databases support this)
- Expected outcomes (yield, success rate if available)
- Rationale explaining why this is recommended
- Practical notes (safety, handling, alternatives)

## Critical Rules:
- NEVER skip Phase 1 (structural analysis)
- ALWAYS query at least 3 different database sources
- CITE which tools provided each piece of evidence
- If databases disagree, discuss the differences
- Prioritize conditions with multiple independent sources of support
- Be specific: no vague ranges - give actionable numbers

## Output Format:
Structure your response clearly:
1. Reaction Analysis Summary
2. Evidence from Each Database (what each source says)
3. Consensus & Contradictions
4. Top 3 Ranked Recommendations (with complete details)
5. Practical Considerations

Remember: Your goal is COMPREHENSIVE analysis, not quick answers. Use ALL available tools systematically."""

    agent = create_agent(llm, tools, system_prompt=system_prompt)
    
    # Comprehensive user query demanding thorough analysis
    user_query = f"""Perform a COMPREHENSIVE analysis of this reaction using ALL available tools and databases:

Reaction SMILES: {args.reaction}

REQUIRED ANALYSIS STEPS:

## Phase 1: Structural Analysis (USE ALL 4 TOOLS)
1. Detect reaction family (detect_reaction_family_tool)
2. Classify each reactant type (classify_reactant_tool for each substrate)
3. Identify functional groups (get_functional_groups_tool for key reactants)
4. Analyze bond changes (analyze_bond_changes_tool)

## Phase 2: Multi-Database Search (QUERY ALL SOURCES)

### A. Rule Database
- Get rule-based conditions (rule_based_conditions_tool)

### B. ML Precedent Databases
- Get ML recommendations (recommend_conditions_tool with top_k={args.top_k})
- Search precedents (search_precedents_tool with k={args.top_k})
- Try unified search (unified_recommender_tool)

### C. Protocol Database
- Search experimental protocols (protocol_recommendation_tool)

### D. HTE Database (66K+ experiments)
- Get HTE recommendations (hte_recommend_tool with top_k={args.top_k})
- Get HTE analytics (hte_analytics_tool for catalyst/metal statistics)

### E. Supporting Analysis
- List available catalysts (list_supported_cores_tool)
- Get dataset statistics (reaction_dataset_analytics_tool)

## Phase 3: Evidence Synthesis
After collecting data from ALL sources above:
1. Create a table comparing what each database recommends
2. Identify CONSENSUS conditions (appearing in multiple sources)
3. Note CONFLICTS or alternatives where databases disagree
4. Rank based on evidence strength and practical feasibility

## Phase 4: Final Recommendations
Provide your top {args.max_protocols} recommendations, each with:

### For Each Recommendation:
- **Complete Setup**: Catalyst (with loading %), ligand (with %), base (with eq), solvent, additives
- **Conditions**: Temperature (°C), time (h), concentration (M), atmosphere, scale considerations
- **Evidence**: Which databases support this (cite tool names and key stats)
- **Expected Outcome**: Yield % and/or success rate from the evidence
- **Rationale**: Why this is recommended (mechanistic logic + empirical support)
- **Practical Notes**: Cost, availability, safety, handling requirements, troubleshooting tips

### Additional Requirements:
- Compare your recommendations: explain pros/cons of each
- If databases strongly disagree, discuss why and which to trust
- Provide a decision tree: "Use conditions 1 if X, use conditions 2 if Y"
- Include alternative approaches if primary recommendations fail

## Critical Instructions:
✓ DO use EVERY tool in Phases 1 and 2 - comprehensive data gathering is mandatory
✓ DO cite specific tool outputs as evidence (e.g., "HTE database shows 67% success rate")
✓ DO give specific numbers (not ranges like "2-3 eq" but "2.5 equivalents")
✓ DO explain contradictions between sources
✗ DON'T skip any database queries
✗ DON'T provide vague conditions
✗ DON'T guess - use the tools to get actual data

Begin your comprehensive analysis now. A chemist should be able to confidently execute your recommended protocol based on your thorough evidence."""

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
