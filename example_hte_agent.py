#!/usr/bin/env python3
"""
Example: Using HTE Tools with Chem Assistant LangChain Agent

This script demonstrates how to use the integrated HTE recommendation and
analytics tools through a LangChain agent for natural language chemistry queries.

Requirements:
    pip install langchain langgraph langchain-openai
    
Environment:
    export OPENAI_API_KEY=your_key_here  # Or set in .env
"""

import sys
import os
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

print("=" * 80)
print("HTE Tools + Chem Assistant Demo")
print("=" * 80)

# Check for OpenAI API key
if not os.getenv("OPENAI_API_KEY"):
    print("\n⚠️  Warning: OPENAI_API_KEY not set")
    print("Set it with: export OPENAI_API_KEY=your_key_here")
    print("\nSkipping LangChain agent demo. Showing direct tool usage instead.\n")
    USE_AGENT = False
else:
    USE_AGENT = True

# Import tools
from chem_assistant.chemtools_wrapper import (
    CHEMTOOLS_TOOLS,
    hte_recommend_tool,
    hte_analytics_tool,
    HTE_AVAILABLE
)

print(f"\n✓ Loaded {len(CHEMTOOLS_TOOLS)} chemistry tools")
print(f"✓ HTE tools available: {HTE_AVAILABLE}\n")

if not HTE_AVAILABLE:
    print("❌ HTE tools not available. Install with:")
    print("   pip install chemtools[hte]")
    sys.exit(1)

# ============================================================================
# Part 1: Direct Tool Usage (No LLM)
# ============================================================================

print("=" * 80)
print("Part 1: Direct Tool Usage")
print("=" * 80)

print("\n1. Suzuki Coupling Recommendation:")
print("   Query: Recommend Pd catalysts for bromobenzene + phenylboronic acid")
print("-" * 80)

result = hte_recommend_tool.invoke({
    "reactant_a_smiles": "Brc1ccccc1",        # Bromobenzene
    "reactant_b_smiles": "c1ccc(B(O)O)cc1",  # Phenylboronic acid
    "top_k": 3,
    "catalyst_filter": "Pd"
})

if result["success"]:
    print(f"   Predicted reaction: {result['predicted_reaction_type']}")
    print(f"   Matching experiments: {result['matching_experiments']}")
    print(f"\n   Top 3 recommendations:")
    
    for i, rec in enumerate(result["recommendations"], 1):
        print(f"\n   {i}. {rec['catalyst']} + {rec['base']} in {rec['solvent']}")
        print(f"      Success: {rec['success_rate']}%, Yield: {rec['avg_yield']}%")
        print(f"      Confidence: {rec['confidence_score']} ({rec['num_experiments']} experiments)")
else:
    print(f"   Error: {result['error']}")

print("\n\n2. Database Analytics - Top Pd Catalysts for Suzuki:")
print("-" * 80)

result = hte_analytics_tool.invoke({
    "query_type": "catalysts",
    "reaction_type": "Suzuki",
    "top_n": 5
})

if result["success"]:
    print(f"   Found {result['total_results']} catalysts\n")
    
    for i, cat in enumerate(result["results"], 1):
        print(f"   {i}. {cat['catalyst']} ({cat['metal']})")
        print(f"      {cat['experiments']} experiments, {cat['success_rate']}% success")
else:
    print(f"   Error: {result['error']}")

print("\n\n3. Metal Usage Comparison:")
print("-" * 80)

result = hte_analytics_tool.invoke({
    "query_type": "metals",
    "top_n": 5
})

if result["success"]:
    print(f"   Total experiments: {result['total_experiments']:,}\n")
    
    for metal_data in result["metals"][:3]:
        print(f"   {metal_data['metal']}: {metal_data['experiments']:,} exp ({metal_data['percentage']}%)")
        if "top_reactions" in metal_data:
            top_rxn = metal_data["top_reactions"][0]
            print(f"      Top reaction: {top_rxn['reaction']} ({top_rxn['count']} exp)")
else:
    print(f"   Error: {result['error']}")

# ============================================================================
# Part 2: LangChain Agent Usage (With LLM)
# ============================================================================

if USE_AGENT:
    print("\n\n" + "=" * 80)
    print("Part 2: LangChain Agent Usage (Natural Language Queries)")
    print("=" * 80)
    
    try:
        from langchain_openai import ChatOpenAI
        from langgraph.prebuilt import create_react_agent
        
        print("\n🤖 Initializing LangChain agent with GPT-4...")
        llm = ChatOpenAI(model="gpt-4", temperature=0)
        agent = create_react_agent(llm, CHEMTOOLS_TOOLS)
        
        print("✓ Agent created with 25 chemistry tools\n")
        
        # Query 1: Simple recommendation
        print("=" * 80)
        print("Query 1: Natural Language Recommendation")
        print("=" * 80)
        print('User: "I need to run a Suzuki coupling with bromobenzene and')
        print('       phenylboronic acid. What palladium catalyst should I use?"')
        print("-" * 80)
        
        response = agent.invoke({
            "messages": [
                "I need to run a Suzuki coupling with bromobenzene and phenylboronic acid. "
                "What palladium catalyst should I use? Give me the top recommendation with success rate."
            ]
        })
        
        print("\nAgent Response:")
        print(response["messages"][-1].content)
        
        # Query 2: Database exploration
        print("\n\n" + "=" * 80)
        print("Query 2: Database Exploration")
        print("=" * 80)
        print('User: "What are the most commonly tested palladium catalysts')
        print('       for Suzuki reactions? Show me the top 3."')
        print("-" * 80)
        
        response = agent.invoke({
            "messages": [
                "What are the most commonly tested palladium catalysts for Suzuki reactions? "
                "Show me the top 3 with their success rates and number of experiments."
            ]
        })
        
        print("\nAgent Response:")
        print(response["messages"][-1].content)
        
        # Query 3: Catalyst comparison
        print("\n\n" + "=" * 80)
        print("Query 3: Catalyst Metal Comparison")
        print("=" * 80)
        print('User: "Compare palladium vs copper for C-N coupling.')
        print('       Which one performs better?"')
        print("-" * 80)
        
        response = agent.invoke({
            "messages": [
                "Compare palladium vs copper catalysts for C-N coupling reactions. "
                "Which metal performs better in terms of success rate? Give me statistics."
            ]
        })
        
        print("\nAgent Response:")
        print(response["messages"][-1].content)
        
        # Query 4: Complex analytical query
        print("\n\n" + "=" * 80)
        print("Query 4: Complex Analytical Query")
        print("=" * 80)
        print('User: "What aryl halide and amine combinations have been')
        print('       tested most extensively for C-N coupling?"')
        print("-" * 80)
        
        response = agent.invoke({
            "messages": [
                "What aryl halide and amine combinations have been tested most extensively "
                "for C-N coupling reactions? Show me the top 3 pairs with their statistics."
            ]
        })
        
        print("\nAgent Response:")
        print(response["messages"][-1].content)
        
    except ImportError as e:
        print(f"\n❌ LangChain not available: {e}")
        print("Install with: pip install langchain langgraph langchain-openai")
    except Exception as e:
        print(f"\n❌ Error creating agent: {e}")

else:
    print("\n\n" + "=" * 80)
    print("Part 2: LangChain Agent Demo Skipped")
    print("=" * 80)
    print("\nTo enable agent demo, set OPENAI_API_KEY environment variable:")
    print("   export OPENAI_API_KEY=your_key_here")
    print("\nOr create a .env file with:")
    print("   OPENAI_API_KEY=your_key_here")

# ============================================================================
# Summary
# ============================================================================

print("\n\n" + "=" * 80)
print("Summary")
print("=" * 80)

print(f"""
✓ HTE recommendation tool: Access to 66,308 experiments
✓ HTE analytics tool: 5 query types for database exploration
✓ Integration: 25 total tools available to LangChain agents

Key Capabilities:
  • Fast condition recommendations (<100ms)
  • Catalyst filtering (Pd, Cu, Ni)
  • Statistical confidence with success rates
  • Database analytics (pairs, catalysts, metals)
  • Natural language queries via LLM agent

Next Steps:
  1. Try your own SMILES strings with hte_recommend_tool
  2. Explore the database with hte_analytics_tool
  3. Build custom agents with CHEMTOOLS_TOOLS
  4. Check docs/HTE_ANALYTICS.md for complete API reference

Documentation:
  • API Reference: docs/HTE_ANALYTICS.md
  • Integration Guide: HTE_CHEM_ASSISTANT_INTEGRATION.md
  • Quick Reference: HTE_ANALYTICS_QUICKREF.md
  • Test Script: test_hte_integration.py
""")

print("=" * 80)
