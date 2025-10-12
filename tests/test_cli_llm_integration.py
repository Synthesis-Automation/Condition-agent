"""
Quick test script to verify LLM synthesis CLI integration.
Tests that all arguments are properly wired up.
"""

import sys
sys.path.insert(0, '.')

from scripts.local_recommendation_cli import HAS_LLM

print("=" * 60)
print("CLI LLM Integration Test")
print("=" * 60)

# Check if LLM support is available
print(f"\n1. LLM Support Available: {HAS_LLM}")

if not HAS_LLM:
    print("   ⚠️  LLM modules not installed. Install with:")
    print("      pip install openai python-dotenv")
    print("\n   Skipping functional tests...")
    sys.exit(0)

# Check imports
try:
    from llmtools.clients import LLMClient
    from llmtools.recommendation_llm import synthesize_recommendations_llm
    print("   ✅ LLM imports successful")
except ImportError as e:
    print(f"   ❌ Import error: {e}")
    sys.exit(1)

# Check utils import
try:
    from scripts.recommendation_cli_utils import summarize_llm_synthesis
    print("   ✅ summarize_llm_synthesis imported")
except ImportError as e:
    print(f"   ❌ Import error: {e}")
    sys.exit(1)

# Check local_llm_synthesis function exists
try:
    from scripts.local_recommendation_cli import local_llm_synthesis
    print("   ✅ local_llm_synthesis function available")
except ImportError as e:
    print(f"   ❌ Import error: {e}")
    sys.exit(1)

# Verify function signature
import inspect
sig = inspect.signature(local_llm_synthesis)
params = list(sig.parameters.keys())
print(f"\n2. Function Parameters ({len(params)}):")
for param in params:
    print(f"   - {param}")

expected_params = [
    'reaction', 'ml_result', 'rule_result', 'protocol_result',
    'constraints', 'llm_provider', 'llm_model', 'prompt_version'
]

missing = set(expected_params) - set(params)
if missing:
    print(f"\n   ❌ Missing parameters: {missing}")
    sys.exit(1)
else:
    print("\n   ✅ All expected parameters present")

# Test summarize_llm_synthesis with mock data
print("\n3. Testing summarize_llm_synthesis():")

mock_result = {
    "status": "success",
    "synthesis": {
        "confidence_level": "high",
        "confidence_reasoning": "All sources agree on Pd(PPh3)4",
        "recommended_condition": {
            "catalyst": "Pd(PPh3)4",
            "ligand": "PPh3",
            "solvent": "THF",
            "temperature": "80°C",
            "base": "K2CO3",
            "rationale": "High consensus across all sources"
        },
        "backup_conditions": [
            {
                "catalyst": "Pd(dppf)Cl2",
                "when_to_use": "If Pd(PPh3)4 unavailable"
            }
        ],
        "warnings": ["Test warning"],
        "consensus_analysis": {
            "catalyst": {"agreement": "high"},
            "solvent": {"agreement": "high"}
        }
    },
    "sources_used": {
        "ml_precedents": 1,
        "rule_matches": 1,
        "protocol_procedures": 1
    },
    "llm_metadata": {
        "model": "deepseek-v3.2-exp",
        "tokens": 1480,
        "latency_ms": 43600,
        "processing_time_ms": 45823
    }
}

try:
    summarize_llm_synthesis(mock_result)
    print("\n   ✅ Summary function works correctly")
except Exception as e:
    print(f"\n   ❌ Error: {e}")
    sys.exit(1)

print("\n" + "=" * 60)
print("✅ All tests passed! CLI LLM integration is ready.")
print("=" * 60)
print("\nTo test with a real reaction, run:")
print('  python scripts/local_recommendation_cli.py \\')
print('    --strategy llm \\')
print('    --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"')
print()
