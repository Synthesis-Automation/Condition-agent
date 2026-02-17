"""
Direct test of GPT-5.2 with minimal prompt.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from llmtools.clients import LLMClient

print("Testing GPT-5.2 with simple prompt...")
print("="*80)

client = LLMClient(provider="openai", model="gpt-5.2")

# Test 1: Simple question with reasoning_effort
print("\n1. With reasoning_effort='low':")
try:
    response = client.chat(
        prompt="What is 2+2?",
        reasoning_effort="low",
        max_tokens=100
    )
    print(f"Success!")
    print(f"Model: {response.model}")
    print(f"Content: {response.content}")
    print(f"Tokens: {response.total_tokens}")
except Exception as e:
    print(f"Failed: {e}")

# Test 2: Without reasoning_effort (should fail for GPT-5)
print("\n2. With temperature=0.0 (should fail for GPT-5):")
try:
    response = client.chat(
        prompt="What is 2+2?",
        temperature=0.0,
        max_tokens=100
    )
    print(f"Success!")
    print(f"Model: {response.model}")
    print(f"Content: {response.content}")
except Exception as e:
    print(f"Failed: {e}")

# Test 3: JSON output with reasoning_effort
print("\n3. JSON output with reasoning_effort='low':")
try:
    response = client.chat(
        prompt='Return ONLY valid JSON: {"result": "4"}',
        reasoning_effort="low",
        max_tokens=100
    )
    print(f"Success!")
    print(f"Model: {response.model}")
    print(f"Content: {response.content}")
except Exception as e:
    print(f"Failed: {e}")
