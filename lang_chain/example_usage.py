"""
Simple example demonstrating ChemTools LangChain integration.

This script shows how to use the ChemTools agent for common chemistry tasks
without needing the interactive CLI.

Usage:
    python lang_chain/example_usage.py
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))


def check_api_key():
    """Check if API key is configured."""
    provider = os.getenv("LLM_PROVIDER", "openai")
    
    if provider == "openai":
        if not os.getenv("OPENAI_API_KEY"):
            print("❌ Error: OPENAI_API_KEY not set")
            print("\nSet your API key:")
            print('  $env:OPENAI_API_KEY = "sk-your-key-here"  # PowerShell')
            print('  export OPENAI_API_KEY="sk-your-key-here"  # Bash')
            return False
    elif provider == "aliyun":
        if not os.getenv("ALIYUN_API_KEY"):
            print("❌ Error: ALIYUN_API_KEY not set")
            return False
    
    return True


def example_1_normalize_smiles():
    """Example 1: Normalize SMILES strings."""
    print("\n" + "="*70)
    print("Example 1: Normalize SMILES")
    print("="*70)
    
    from lang_chain.chemtools_agent import quick_query
    
    query = "Normalize this SMILES: c1ccccc1"
    print(f"\nQuery: {query}")
    print("\n🤖 Agent is working...", end="", flush=True)
    response = quick_query(query)
    print("\r" + " " * 30 + "\r", end="")  # Clear the working message
    print("Agent Response:")
    print(response)


def example_2_detect_reaction():
    """Example 2: Detect reaction family."""
    print("\n" + "="*70)
    print("Example 2: Detect Reaction Family")
    print("="*70)
    
    from lang_chain.chemtools_agent import quick_query
    
    query = "What reaction family is this: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    print(f"\nQuery: {query}")
    print("\n🤖 Agent is working...", end="", flush=True)
    response = quick_query(query)
    print("\r" + " " * 30 + "\r", end="")  # Clear the working message
    print("Agent Response:")
    print(response)


def example_3_recommend_conditions():
    """Example 3: Recommend reaction conditions."""
    print("\n" + "="*70)
    print("Example 3: Recommend Conditions")
    print("="*70)
    
    from lang_chain.chemtools_agent import quick_query
    
    query = "Recommend optimal conditions for this Buchwald coupling: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    print(f"\nQuery: {query}")
    print("\n🤖 Agent is working...", end="", flush=True)
    response = quick_query(query)
    print("\r" + " " * 30 + "\r", end="")  # Clear the working message
    print("Agent Response:")
    print(response)


def example_4_reagent_lookup():
    """Example 4: Look up reagent information."""
    print("\n" + "="*70)
    print("Example 4: Reagent Lookup")
    print("="*70)
    
    from lang_chain.chemtools_agent import quick_query
    
    query = "What is the CAS number and role of Cs2CO3?"
    print(f"\nQuery: {query}")
    print("\n🤖 Agent is working...", end="", flush=True)
    response = quick_query(query)
    print("\r" + " " * 30 + "\r", end="")  # Clear the working message
    print("Agent Response:")
    print(response)


def example_5_conversational():
    """Example 5: Multi-turn conversation."""
    print("\n" + "="*70)
    print("Example 5: Conversational Agent")
    print("="*70)
    
    from lang_chain.chemtools_agent import ChemToolsAgent
    
    agent = ChemToolsAgent()
    history = []
    
    queries = [
        "What functional groups are in CCO?",
        "Is it an alcohol?",
        "What about c1ccccc1O?",
    ]
    
    for query in queries:
        print(f"\nYou: {query}")
        response, history = agent.chat(query, history)
        print(f"Agent: {response}")


def example_6_custom_agent():
    """Example 6: Custom agent configuration."""
    print("\n" + "="*70)
    print("Example 6: Custom Agent Configuration")
    print("="*70)
    
    from lang_chain.chemtools_agent import ChemToolsAgent
    
    # Create agent with custom settings
    agent = ChemToolsAgent(
        model="gpt-4o",
        temperature=0.1,
        verbose=True  # Show tool execution
    )
    
    query = "Classify this reactant: Brc1ccccc1"
    print(f"\nQuery: {query}")
    print("\nAgent Response (with verbose output):")
    response = agent.run(query)
    print(response)


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("ChemTools LangChain Integration - Examples")
    print("="*70)
    
    # Check API key
    if not check_api_key():
        sys.exit(1)
    
    print("\nRunning examples...")
    print("Note: Each example will call the LLM API (may take a few seconds)")
    
    try:
        # Run examples
        example_1_normalize_smiles()
        example_2_detect_reaction()
        example_3_recommend_conditions()
        example_4_reagent_lookup()
        example_5_conversational()
        # example_6_custom_agent()  # Uncomment for verbose output
        
        print("\n" + "="*70)
        print("✅ All examples completed successfully!")
        print("="*70)
        print("\nNext steps:")
        print("  - Try the interactive CLI: python -m lang_chain.chemtools_cli")
        print("  - Read the docs: lang_chain/README.md")
        print("  - Build your own agent with custom tools")
        
    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        print("\nTroubleshooting:")
        print("  - Ensure API key is set correctly")
        print("  - Check internet connection")
        print("  - Verify LangChain packages are installed: pip install -r lang_chain/requirements.txt")
        sys.exit(1)


if __name__ == "__main__":
    main()
