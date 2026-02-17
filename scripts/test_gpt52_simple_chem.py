"""
Test GPT-5.2 with chemistry prompt but natural language output.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from llmtools.clients import LLMClient

rxn_smiles = "CC1(C)OB(c2cnn(CCOC3CCCC O3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"
reactants, products = rxn_smiles.split(">>")

# Simple natural language prompt
prompt1 = f"""Analyze this organic chemistry reaction:

Reactants: {reactants}
Products: {products}

What protecting groups were removed or added?"""

# JSON prompt without complexity
prompt2 = f"""Analyze this reaction:

Reactants: {reactants}
Products: {products}

Return as JSON: {{"protecting_groups_removed": ["list"], "reaction_type": "type"}}"""

client = LLMClient(provider="openai", model="gpt-5.2")

print("Test 1: Natural language prompt...")
response1 = client.chat(prompt=prompt1, reasoning_effort="low", max_tokens=500)
print(f"Content length: {len(response1.content)}")
print(f"Content: {response1.content[:200]}")

print("\n\nTest 2: Simple JSON prompt...")
response2 = client.chat(prompt=prompt2, reasoning_effort="low", max_tokens=500)
print(f"Content length: {len(response2.content)}")
print(f"Content: {response2.content[:200]}")
