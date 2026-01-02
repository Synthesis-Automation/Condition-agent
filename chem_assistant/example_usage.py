"""
Simple examples demonstrating ChemTools featurization/analysis integration.
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))


def check_api_key() -> bool:
    """Check if API key is configured."""
    provider = os.getenv("LLM_PROVIDER", "openai")

    if provider == "openai":
        if not os.getenv("OPENAI_API_KEY"):
            print("Error: OPENAI_API_KEY not set")
            print("Set your API key:")
            print('  $env:OPENAI_API_KEY = "sk-your-key-here"  # PowerShell')
            print('  export OPENAI_API_KEY="sk-your-key-here"  # Bash')
            return False
    elif provider == "aliyun":
        if not os.getenv("ALIYUN_API_KEY"):
            print("Error: ALIYUN_API_KEY not set")
            return False

    return True


def example_1_normalize_smiles() -> None:
    """Example 1: Normalize SMILES strings."""
    print("\n" + "=" * 70)
    print("Example 1: Normalize SMILES")
    print("=" * 70)

    from chem_assistant.chemtools_agent import quick_query

    query = "Normalize this SMILES: c1ccccc1"
    print(f"\nQuery: {query}")
    response = quick_query(query)
    print("Agent Response:")
    print(response)


def example_2_analyze_reaction() -> None:
    """Example 2: Analyze reaction (taxonomy + family)."""
    print("\n" + "=" * 70)
    print("Example 2: Analyze Reaction")
    print("=" * 70)

    from chem_assistant.chemtools_agent import quick_query

    query = "Analyze reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
    print(f"\nQuery: {query}")
    response = quick_query(query)
    print("Agent Response:")
    print(response)


def example_3_featurize_molecule() -> None:
    """Example 3: Unified molecule featurization."""
    print("\n" + "=" * 70)
    print("Example 3: Unified Molecule Featurization")
    print("=" * 70)

    from chem_assistant.chemtools_agent import quick_query

    query = "Featurize molecule: c1ccccc1O"
    print(f"\nQuery: {query}")
    response = quick_query(query)
    print("Agent Response:")
    print(response)


def example_4_featurize_reaction_pair() -> None:
    """Example 4: Reaction pair featurization."""
    print("\n" + "=" * 70)
    print("Example 4: Reaction Pair Featurization")
    print("=" * 70)

    from chem_assistant.chemtools_agent import quick_query

    query = "Featurize pair: electrophile=Brc1ccccc1 nucleophile=Nc1ccccc1"
    print(f"\nQuery: {query}")
    response = quick_query(query)
    print("Agent Response:")
    print(response)


def example_5_conversational() -> None:
    """Example 5: Multi-turn conversation."""
    print("\n" + "=" * 70)
    print("Example 5: Conversational Agent")
    print("=" * 70)

    from chem_assistant.chemtools_agent import ChemToolsAgent

    agent = ChemToolsAgent()
    history = []

    queries = [
        "What functional groups are in CCO?",
        "Featurize molecule: c1ccccc1O",
        "Analyze reaction: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    ]

    for query in queries:
        print(f"\nYou: {query}")
        response, history = agent.chat(query, history)
        print(f"Agent: {response}")


def main() -> None:
    """Run all examples."""
    print("\n" + "=" * 70)
    print("ChemTools Featurization Examples")
    print("=" * 70)

    if not check_api_key():
        sys.exit(1)

    print("\nRunning examples...")
    print("Note: Each example will call the LLM API.")

    try:
        example_1_normalize_smiles()
        example_2_analyze_reaction()
        example_3_featurize_molecule()
        example_4_featurize_reaction_pair()
        example_5_conversational()

        print("\n" + "=" * 70)
        print("All examples completed successfully")
        print("=" * 70)
        print("\nNext steps:")
        print("  - Try the CLI: python -m chem_assistant.chemtools_cli")
        print("  - Read the docs: chem_assistant/README.md")
    except Exception as exc:
        print(f"\nError running examples: {exc}")
        print("\nTroubleshooting:")
        print("  - Ensure API key is set correctly")
        print("  - Check internet connection")
        print("  - Verify dependencies: pip install -r requirements.txt")
        sys.exit(1)


if __name__ == "__main__":
    main()
