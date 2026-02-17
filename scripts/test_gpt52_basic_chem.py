"""
Test basic chemistry with GPT-5.2.
"""

import sys
sys.path.insert(0, 'c:/Git-softwares/Condition-agent')

from llmtools.clients import LLMClient

client = LLMClient(provider="openai", model="gpt-5.2")

tests = [
    "Describe the Suzuki coupling reaction.",
    "What is THP?",
    "Explain nucleophilic substitution.",
    "What are protecting groups in chemistry?",
    "List common organic reactions.",
]

for i, prompt in enumerate(tests, 1):
    print(f"\nTest {i}: {prompt[:50]}...")
    response = client.chat(prompt=prompt, reasoning_effort="low", max_tokens=200)
    print(f"Length: {len(response.content)} | Content: {response.content[:100] if response.content else '(empty)'}")
