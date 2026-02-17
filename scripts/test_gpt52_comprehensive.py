"""
Test GPT-5.2 comprehensive mode separately.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools.quick_reaction_glance import quick_reaction_glance
from llmtools.clients import LLMClient

# Your reaction
rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

print("Testing GPT-5 models comprehensive mode...")
print("="*80)
print()

# Try different GPT-5 models
models_to_test = ["gpt-5.2", "gpt-5-mini", "gpt-4o"]

for model_name in models_to_test:
    print(f"\n--- Testing {model_name} ---")
    try:
        client = LLMClient(provider="openai", model=model_name)
        print(f"✓ {model_name} client created")

        print("Running comprehensive analysis...")
        result = quick_reaction_glance(
            rxn_smiles,
            client,
            prompt_style="comprehensive",
            thorough=True
        )

        if result.get('success'):
            print("✓ Analysis successful!")
            print()
            print(f"Summary: {result.get('summary', 'N/A')}")
            print()

            # Protecting groups
            pg = result.get('protecting_groups', {})
            if pg.get('removed'):
                print(f"✓ Detected PG removal: {pg['removed']}")
            if pg.get('added'):
                print(f"✓ Detected PG addition: {pg['added']}")

            # All changes
            changes = result.get('all_changes', [])
            if changes:
                print()
                print("All changes detected:")
                for change in changes:
                    print(f"  • {change}")

            # Reaction types
            print()
            print(f"Reaction types: {result.get('reaction_types', [])}")
            print(f"Complexity: {result.get('complexity', 'N/A')}")
            print(f"Confidence: {result.get('confidence', 0.0):.2f}")

            if 'metadata' in result:
                print()
                print(f"Model: {result['metadata'].get('model', 'N/A')}")
                print(f"Latency: {result['metadata'].get('latency_ms', 0):.0f}ms")
                print(f"Tokens: {result['metadata'].get('total_tokens', 0)}")

            # If successful, no need to test other models
            print(f"\n✓ {model_name} works! Using this model.")
            break

        else:
            print(f"✗ Failed: {result.get('error', 'Unknown error')}")
            if 'raw_response' in result:
                print()
                print("Raw response (first 500 chars):")
                print(result['raw_response'][:500] if result['raw_response'] else "(empty)")

    except Exception as e:
        print(f"✗ Error: {e}")
        continue  # Try next model
