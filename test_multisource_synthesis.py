"""
Test Multi-Source LLM Synthesis

Tests the synthesize_recommendations_llm() function with diverse reactions
to validate the multi-source synthesis approach.

Usage:
    python test_multisource_synthesis.py

Author: Condition-Agent Team
Date: October 12, 2025
"""

import sys
import json
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from llmtools.clients import LLMClient
from llmtools.recommendation_llm import synthesize_recommendations_llm


# Test cases with simulated results from different sources
TEST_CASES = [
    {
        "name": "Suzuki Coupling - High Consensus",
        "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "ml_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(PPh3)4",
                    "solvent": "toluene",
                    "temperature": "80°C",
                    "base": "K3PO4",
                    "similarity": 0.92
                },
                {
                    "catalyst": "Pd(PPh3)4",
                    "solvent": "DMF",
                    "temperature": "100°C",
                    "base": "K2CO3",
                    "similarity": 0.88
                }
            ]
        },
        "rule_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(PPh3)4",
                    "solvent": "toluene",
                    "temperature": "80°C",
                    "base": "K3PO4"
                }
            ]
        },
        "protocol_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(PPh3)4",
                    "solvent": "toluene",
                    "temperature": "80-90°C",
                    "base": "K3PO4",
                    "time": "12-16h"
                }
            ]
        },
        "constraints": None,
        "expected_confidence": "high"
    },
    {
        "name": "Suzuki Coupling - Nitro Substrate (Medium Consensus)",
        "reaction_smiles": "Brc1ccc([N+](=O)[O-])cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc([N+](=O)[O-])cc2)cc1",
        "ml_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(PPh3)4",
                    "solvent": "toluene",
                    "temperature": "80°C",
                    "base": "K3PO4",
                    "similarity": 0.85
                }
            ]
        },
        "rule_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(dppf)Cl2",
                    "solvent": "DMF",
                    "temperature": "100°C",
                    "base": "K3PO4"
                }
            ]
        },
        "protocol_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(dppf)Cl2",
                    "solvent": "DMF",
                    "temperature": "100°C",
                    "base": "K3PO4",
                    "time": "18-24h",
                    "note": "Electron-poor substrates require electron-rich ligands"
                }
            ]
        },
        "constraints": {"scale": "multigram", "cost": "low"},
        "expected_confidence": "medium"
    },
    {
        "name": "Buchwald-Hartwig - Low Consensus",
        "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "ml_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd2(dba)3",
                    "ligand": "BINAP",
                    "solvent": "toluene",
                    "temperature": "100°C",
                    "base": "NaOtBu",
                    "similarity": 0.78
                }
            ]
        },
        "rule_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(OAc)2",
                    "ligand": "XPhos",
                    "solvent": "dioxane",
                    "temperature": "90°C",
                    "base": "Cs2CO3"
                }
            ]
        },
        "protocol_results": {
            "recommended_conditions": [
                {
                    "catalyst": "Pd(dba)2",
                    "ligand": "SPhos",
                    "solvent": "toluene",
                    "temperature": "80°C",
                    "base": "NaOtBu"
                }
            ]
        },
        "constraints": None,
        "expected_confidence": "low"
    }
]


def run_test(test_case: dict, llm_client: LLMClient):
    """Run a single test case."""
    print(f"\n{'='*80}")
    print(f"TEST: {test_case['name']}")
    print(f"{'='*80}")
    print(f"Reaction: {test_case['reaction_smiles'][:60]}...")
    print(f"Constraints: {test_case['constraints']}")
    print(f"Expected confidence: {test_case['expected_confidence']}")
    
    # Run synthesis
    result = synthesize_recommendations_llm(
        reaction_smiles=test_case['reaction_smiles'],
        ml_results=test_case['ml_results'],
        rule_results=test_case['rule_results'],
        protocol_results=test_case['protocol_results'],
        constraints=test_case['constraints'],
        llm_client=llm_client
    )
    
    # Check status
    if result['status'] != 'success':
        print(f"❌ FAILED: {result.get('error')}")
        if 'raw_response' in result:
            print(f"\nRaw LLM response:\n{result['raw_response'][:500]}...")
        return False
    
    synthesis = result['synthesis']
    
    # Display results
    print(f"\n✅ SUCCESS")
    print(f"\nConfidence Level: {synthesis['confidence_level']}")
    print(f"Confidence Reasoning: {synthesis['confidence_reasoning'][:100]}...")
    
    print(f"\n--- Consensus Analysis ---")
    for param, analysis in synthesis['consensus_analysis'].items():
        print(f"{param}: {analysis['agreement']} - {analysis.get('consensus_value', 'N/A')}")
    
    print(f"\n--- Recommended Condition ---")
    rec = synthesis['recommended_condition']
    print(f"Catalyst: {rec.get('catalyst')}")
    print(f"Ligand: {rec.get('ligand')}")
    print(f"Solvent: {rec.get('solvent')}")
    print(f"Temperature: {rec.get('temperature')}")
    print(f"Base: {rec.get('base')}")
    print(f"\nRationale: {rec.get('rationale', 'N/A')[:150]}...")
    
    print(f"\n--- Backup Conditions ---")
    for i, backup in enumerate(synthesis.get('backup_conditions', []), 1):
        print(f"{i}. {backup.get('catalyst', 'N/A')} in {backup.get('solvent', 'N/A')}")
        print(f"   When: {backup.get('when_to_use', 'N/A')[:100]}...")
    
    print(f"\n--- Warnings ---")
    for warning in synthesis.get('warnings', []):
        print(f"⚠️  {warning}")
    
    print(f"\n--- Source Contributions ---")
    comp = synthesis.get('source_comparison', {})
    print(f"ML: {comp.get('ml_contribution', 'N/A')[:80]}...")
    print(f"Rule: {comp.get('rule_contribution', 'N/A')[:80]}...")
    print(f"Protocol: {comp.get('protocol_contribution', 'N/A')[:80]}...")
    
    # Validate confidence matches expectation
    expected = test_case['expected_confidence']
    actual = synthesis['confidence_level']
    if expected == actual:
        print(f"\n✓ Confidence level matches expected: {expected}")
    else:
        print(f"\n⚠️  Confidence mismatch - Expected: {expected}, Got: {actual}")
    
    # Display LLM metadata
    meta = result['llm_metadata']
    print(f"\nLLM: {meta['model']} | Tokens: {meta['tokens']} | Latency: {meta['latency_ms']}ms")
    
    return True


def main():
    """Run all tests."""
    print("="*80)
    print("MULTI-SOURCE LLM SYNTHESIS TEST SUITE")
    print("="*80)
    
    # Initialize LLM client
    try:
        llm_client = LLMClient(provider="aliyun", model="deepseek-v3.2-exp")
        print(f"✓ LLM client initialized: {llm_client.model}")
    except Exception as e:
        print(f"❌ Failed to initialize LLM client: {e}")
        return
    
    # Run tests
    passed = 0
    failed = 0
    
    for test_case in TEST_CASES:
        try:
            success = run_test(test_case, llm_client)
            if success:
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"\n❌ TEST EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    # Summary
    print(f"\n{'='*80}")
    print(f"TEST SUMMARY")
    print(f"{'='*80}")
    print(f"Passed: {passed}/{len(TEST_CASES)}")
    print(f"Failed: {failed}/{len(TEST_CASES)}")
    
    if failed == 0:
        print(f"\n🎉 ALL TESTS PASSED!")
    else:
        print(f"\n⚠️  Some tests failed. Review output above.")


if __name__ == "__main__":
    main()
