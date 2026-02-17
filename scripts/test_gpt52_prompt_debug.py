"""
Test GPT-5.2 with the actual comprehensive prompt to see what's wrong.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from llmtools.clients import LLMClient
from chemtools.quick_reaction_glance import _get_comprehensive_prompt

rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

reactants, products = rxn_smiles.split(">>")

prompt = _get_comprehensive_prompt(reactants, products)

print("Testing GPT-5.2 with comprehensive chemistry prompt...")
print("="*80)
print(f"\nPrompt length: {len(prompt)} characters")
print(f"Approx tokens: {len(prompt) // 4}")  # Rough estimate
print(f"\nFirst 500 chars of prompt:")
print(prompt[:500])
print("...")
print(f"\nLast 200 chars of prompt:")
print(prompt[-200:])

client = LLMClient(provider="openai", model="gpt-5.2")

print("\n\nCalling GPT-5.2 (Test 1: WITH reasoning_effort)...")
try:
    response = client.chat(
        prompt=prompt,
        reasoning_effort="low",
        max_tokens=1000
    )
    print(f"✓ Success!")
    print(f"Model: {response.model}")
    print(f"Content length: {len(response.content)} chars")
    print(f"Content: {response.content[:500] if response.content else '(empty)'}")
except Exception as e:
    print(f"✗ Failed: {e}")

print("\n\nCalling GPT-5.2 (Test 2: WITHOUT reasoning_effort)...")
try:
    response2 = client.chat(
        prompt=prompt,
        temperature=0.0,
        max_tokens=1000
    )
    print(f"✓ Success!")
    print(f"Model: {response2.model}")
    print(f"Content length: {len(response2.content)} chars")
    print(f"Content: {response2.content[:500] if response2.content else '(empty)'}")
except Exception as e:
    print(f"✗ Failed: {e}")
