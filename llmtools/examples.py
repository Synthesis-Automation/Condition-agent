"""
LLM Tools Examples
==================

Demonstrates various llmtools capabilities:
1. Basic LLM client usage
2. Chemistry agent recommendations
3. Mechanism explanations
4. Protocol generation
5. Troubleshooting assistance

Run:
    python llmtools/examples.py
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from llmtools import LLMClient, ChemistryAgent


def example_1_basic_client():
    """Example 1: Basic LLM client usage."""
    print("\n" + "="*70)
    print("Example 1: Basic LLM Client")
    print("="*70)
    
    # Create client (reads API key from environment)
    client = LLMClient(provider="openai", model="gpt-4o-mini")
    
    # Simple query
    response = client.chat("Suggest a good solvent for Suzuki coupling reactions.")
    
    print("\nPrompt: Suggest a good solvent for Suzuki coupling reactions.")
    print(f"\nResponse from {response.model}:")
    print("-" * 70)
    print(response.content)
    print("-" * 70)
    print(f"\nTokens: {response.total_tokens} | Latency: {response.latency_ms:.0f}ms")


def example_2_condition_recommendation():
    """Example 2: Recommend reaction conditions with precedents."""
    print("\n" + "="*70)
    print("Example 2: Condition Recommendation (LLM + ChemTools)")
    print("="*70)
    
    client = LLMClient(provider="openai", model="gpt-4o")
    agent = ChemistryAgent(client, use_chemtools=True, verbose=True)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "Buchwald_CN"
    
    print(f"\nReaction: {reaction}")
    print(f"Type: {reaction_type}")
    print("\nGenerating recommendations...")
    
    result = agent.suggest_conditions(
        reaction=reaction,
        reaction_type=reaction_type,
        context="Prefer mild conditions suitable for scale-up"
    )
    
    print("\n" + "-" * 70)
    print(result["llm_response"])
    print("-" * 70)
    print(f"\nPrecedents used: {result['precedents_used']}")
    print(f"Model: {result['model']}")
    print(f"Tokens: {result['tokens']}")


def example_3_mechanism_explanation():
    """Example 3: Explain reaction mechanism."""
    print("\n" + "="*70)
    print("Example 3: Mechanism Explanation")
    print("="*70)
    
    client = LLMClient(provider="openai", model="gpt-4o-mini")
    agent = ChemistryAgent(client)
    
    reaction = "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
    reaction_type = "Suzuki"
    
    print(f"\nReaction: {reaction}")
    print(f"Type: {reaction_type}")
    print("\nGenerating mechanism explanation...")
    
    result = agent.explain_mechanism(
        reaction=reaction,
        reaction_type=reaction_type,
        reagents="Pd(PPh3)4, K2CO3",
        detail_level="Concise overview with key steps"
    )
    
    print("\n" + "-" * 70)
    print(result["mechanism_explanation"])
    print("-" * 70)


def example_4_troubleshooting():
    """Example 4: Troubleshoot a problematic reaction."""
    print("\n" + "="*70)
    print("Example 4: Reaction Troubleshooting")
    print("="*70)
    
    client = LLMClient(provider="openai", model="gpt-4o")
    agent = ChemistryAgent(client)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "Buchwald_CN"
    problem = "Getting only 15% yield after 24 hours at 100°C"
    
    print(f"\nReaction: {reaction}")
    print(f"Type: {reaction_type}")
    print(f"Problem: {problem}")
    print("\nGenerating troubleshooting advice...")
    
    result = agent.troubleshoot_reaction(
        reaction=reaction,
        reaction_type=reaction_type,
        problem=problem,
        current_conditions="Pd(PPh3)4 (5 mol%), K2CO3 (2 eq), toluene, 100°C",
        observed_result="15% conversion by HPLC, mostly starting material"
    )
    
    print("\n" + "-" * 70)
    print(result["troubleshooting_advice"])
    print("-" * 70)


def example_5_protocol_generation():
    """Example 5: Generate experimental protocol."""
    print("\n" + "="*70)
    print("Example 5: Protocol Generation")
    print("="*70)
    
    client = LLMClient(provider="openai", model="gpt-4o-mini")
    agent = ChemistryAgent(client)
    
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
    reaction_type = "Buchwald_CN"
    
    conditions = {
        "catalyst": "Pd(OAc)2 (2 mol%)",
        "ligand": "BINAP (4 mol%)",
        "base": "NaOtBu (1.5 eq)",
        "solvent": "Toluene (0.5 M)",
        "temperature": "100°C",
        "time": "16 h"
    }
    
    print(f"\nReaction: {reaction}")
    print(f"Type: {reaction_type}")
    print(f"Scale: 5 mmol")
    print("\nGenerating experimental protocol...")
    
    result = agent.generate_protocol(
        reaction=reaction,
        reaction_type=reaction_type,
        conditions=conditions,
        scale="5 mmol"
    )
    
    print("\n" + "-" * 70)
    print(result["protocol"])
    print("-" * 70)


def check_api_keys():
    """Check if required API keys are set."""
    openai_key = os.getenv("OPENAI_API_KEY")
    aliyun_key = os.getenv("ALIYUN_API_KEY")
    
    if not openai_key and not aliyun_key:
        print("\n⚠️  Warning: No API keys found!")
        print("\nPlease set at least one using PowerShell:")
        print("  [System.Environment]::SetEnvironmentVariable('OPENAI_API_KEY', 'sk-...', 'User')")
        print("\nThen restart PowerShell and try again.")
        print("\nSee docs/WINDOWS_ENV_SETUP.md for detailed instructions.")
        return False
    
    if openai_key:
        print(f"\n✓ OPENAI_API_KEY found: {openai_key[:10]}...")
    if aliyun_key:
        print(f"✓ ALIYUN_API_KEY found: {aliyun_key[:10]}...")
    
    return True


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("LLM Tools Examples")
    print("="*70)
    
    if not check_api_keys():
        print("\nSkipping examples due to missing API keys.")
        return 1
    
    try:
        # Run examples
        example_1_basic_client()
        
        # Uncomment to run more examples (costs tokens!)
        # example_2_condition_recommendation()
        # example_3_mechanism_explanation()
        # example_4_troubleshooting()
        # example_5_protocol_generation()
        
        print("\n" + "="*70)
        print("Examples completed!")
        print("="*70)
        print("\nTo run more examples, uncomment them in llmtools/examples.py")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
