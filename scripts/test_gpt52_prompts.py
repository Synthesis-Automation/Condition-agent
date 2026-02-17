"""
Quick test to see if GPT-5.2 works with different prompt styles.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from chemtools.quick_reaction_glance import quick_reaction_glance
from llmtools.clients import LLMClient

rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"

client = LLMClient(provider="openai", model="gpt-5.2")

print("Testing GPT-5.2 with different prompt styles...")
print("="*80)

# Test 1: Structured prompt (simpler)
print("\n1. Testing 'structured' prompt style...")
result1 = quick_reaction_glance(rxn_smiles, client, prompt_style="structured", thorough=False)
print(f"Success: {result1.get('success')}")
if not result1.get('success'):
    print(f"Error: {result1.get('error')}")

# Test 2: Comprehensive prompt
print("\n2. Testing 'comprehensive' prompt style (thorough=True)...")
result2 = quick_reaction_glance(rxn_smiles, client, prompt_style="comprehensive", thorough=True)
print(f"Success: {result2.get('success')}")
if not result2.get('success'):
    print(f"Error: {result2.get('error')}")

# Test 3: Concise prompt
print("\n3. Testing 'concise' prompt style...")
result3 = quick_reaction_glance(rxn_smiles, client, prompt_style="concise", thorough=False)
print(f"Success: {result3.get('success')}")
if not result3.get('success'):
    print(f"Error: {result3.get('error')}")
