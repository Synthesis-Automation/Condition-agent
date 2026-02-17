"""
Test if GPT-5.2 refuses SMILES strings specifically.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from llmtools.clients import LLMClient

client = LLMClient(provider="openai", model="gpt-5.2")

# Test 1: With SMILES
print("Test 1: With SMILES strings...")
prompt1 = """Analyze this: CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C"""
response1 = client.chat(prompt=prompt1, reasoning_effort="low", max_tokens=100)
print(f"Content length: {len(response1.content)}")
print(f"Content: {response1.content}")

# Test 2: Without SMILES (plain chemistry)
print("\n\nTest 2: Plain English chemistry question...")
prompt2 = """What is a THP protecting group and how is it removed?"""
response2 = client.chat(prompt=prompt2, reasoning_effort="low", max_tokens=100)
print(f"Content length: {len(response2.content)}")
print(f"Content: {response2.content}")

# Test 3: Single SHORT SMILES
print("\n\nTest 3: Very short SMILES...")
prompt3 = """What functional groups are in this molecule: CCO"""
response3 = client.chat(prompt=prompt3, reasoning_effort="low", max_tokens=100)
print(f"Content length: {len(response3.content)}")
print(f"Content: {response3.content}")
